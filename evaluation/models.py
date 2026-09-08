"""
All regression model architectures. 
Each model has 
    - fit (X, Y) where X: (n_compounds, d)  Y: (n_compounds, n_genes)
    - predict(X) -> (n_compounds, n_genes)

Y: delta, mean expression under the perturbation - the control mean

adapted from the LLM-Pert repo (simple_models.py, GP/multioutputGP.py, GP/kernels.py)   
note: compounds keyed by InChIKey
"""

from __future__ import annotations
from dataclasses import dataclass, field
import dataclasses
from functools import partial
import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize
from sklearn.decomposition import PCA
from sklearn.linear_model import RidgeCV 
from sklearn.model_selection import KFold
from sklearn.neighbors import NearestNeighbors
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import make_pipeline
from sklearn.cross_decomposition import PLSRegression


# Baseline Models
class ControlBaseline:
    """ predicts no change = zeros in the delta space"""

    def fit(self, X, Y):
        self.n_genes = Y.shape[1]
        return self

    def predict(self, X):
        preds = np.zeros((len(X), self.n_genes), dtype = np.float32)
        return preds


class MeanBaseline:
    """ predicts the average training response (embeddings not taken into account) """

    def fit(self, X, Y):
        self.mean_ = Y.mean(axis = 0)
        return self

    def predict(self, X):
        preds = np.tile(self.mean_, (len(X), 1))
        return preds


@dataclass
class RandomBaseline:
    """ predicts a randomly chosen training compound's response """
    seed: int = 10

    def fit(self, X, Y):
        self.Y_ = Y
        return self

    def predict(self, X):
        rng = np.random.default_rng(self.seed)
        preds = self.Y_[rng.integers(0, len(self.Y_), size = len(X))]
        return preds

# KNN
@dataclass
class KNN:
    # tested hyperparams, selected!
    k: int = 20
    metric: str = "cosine"

    def fit(self, X, Y):
        self.mu_ = X.mean(axis = 0)              
        self.Y_ = Y
        self.k_ = min(self.k, len(X))
        self.nn_ = NearestNeighbors(
            n_neighbors = self.k_, metric = self.metric
        ).fit(X - self.mu_)
        return self
    
    def predict(self, X):
        _, idx = self.nn_.kneighbors(X - self.mu_, n_neighbors = self.k_)
        preds = self.Y_[idx].mean(axis = 1)
        return preds

# Ridge 
@dataclass
class Ridge:
    """ multi-output ridge. alpha by exact leave-one-out on the training fold."""

    alphas: tuple = tuple(np.logspace(-2, 8, 41))
    with_std: bool = True

    def fit(self, X, Y):
        self.scaler_ = StandardScaler(with_std = self.with_std).fit(X)
        self.model_ = RidgeCV(alphas = self.alphas, 
                              alpha_per_target = True).fit(self.scaler_.transform(X), Y)
        self.alpha_ = self.model_.alpha_       
        return self

    def predict(self, X):
        return self.model_.predict(self.scaler_.transform(X))


# MLPs
@dataclass
class MLP:
    """ sklearn MLPRegressor w/ early stopping """

    # tested hyperparams, selected!
    hidden_layer_sizes: tuple = (2048, 1024)
    alpha: float = 0.00001
                     
    lr: float = 1e-3
    max_iter: int = 1500
    validation_fraction: float = 0.2          
    with_std: bool = True                    
    seed: int = 17

    def fit(self, X, Y):
        self.scaler_ = StandardScaler(with_std = self.with_std).fit(X)
        self.model_ = MLPRegressor(
            hidden_layer_sizes = self.hidden_layer_sizes,
            alpha = self.alpha,
            learning_rate_init = self.lr,
            max_iter = self.max_iter,
            early_stopping = True,
            validation_fraction = self.validation_fraction,
            n_iter_no_change = 20,
            random_state=self.seed,
        ).fit(self.scaler_.transform(X), Y)
        return self

    def predict(self, X):
        preds = self.model_.predict(self.scaler_.transform(X))
        return preds


# Gaussian Porcess Model from the LLMPert Repo
@dataclass
class GaussianProcess:
    kernel: str = "rbf"
    reduce_dim: int | None = None # tested w hyperparam config run 
    
    prior_mean: str = "train_mean"
    feature_dims: tuple | None = None
    with_std: bool = True
    noise_frac: float = 0.1
    jitter: float = 1e-8
    verbose: bool = True
    
    
    # model functions 
    def fit(self, X, Y):
        """
        - setup X, Y, prior mean offset, hyperparameters
        - compute the kernel matrix K, add noise, Cholesky factorization, alpha = K^-1 Yc
        - optimise the hyperparameter vector by min. the neg. log marginal likelihood
        """
        X = self.fit_reduce(X)
        
        if self.prior_mean == "train_mean":
            self.offset_ = Y.mean(axis = 0)
        else:
            self.offset_ = 0.0
        Yc = Y - self.offset_
 
        theta0, bounds = self.initial_hyperparams(X, Yc)
        result = minimize(self.negative_log_marginal_likelihood, theta0,
                          args = (X, Yc), 
                          method = "L-BFGS-B", 
                          bounds = bounds)
        self.theta = result.x

        if self.verbose and np.allclose(result.x, theta0, atol = 1e-6):
            print("[gp] hyperparameters didn't changed from the initialised values")
            
        training_kernel = self.calc_K(X, X)
        noise_variance = self.gp_hyperparams()[2]
        I = np.eye(len(X))
        K = training_kernel + noise_variance * I
        self.chol = cho_factor(K, lower = True)
        self.alpha = cho_solve(self.chol, Yc)
        
        self.X = X
        self.Yc = Yc
        return self
 
    def predict(self, X):
        """ calculate the GP posterior mean for new compounds """
        kernel = self.calc_K(self.reduce(X), self.X)
        y_pred = kernel @ self.alpha + self.offset_
        return y_pred
 
    def predict_with_weights(self, X):
        """ same as predict, but also returns the weight matrix, adapted from LLMPert """
        kernel = self.calc_K(self.reduce(X), self.X)
        weights = kernel @ cho_solve(self.chol, np.eye(len(self.X)))
        y_pred = weights @ self.Yc + self.offset_
        return y_pred, weights
 
    @property
    def block_lengthscales(self) -> np.ndarray:
        return self.gp_hyperparams()[0]
 
    @property
    def block_amplitudes(self) -> np.ndarray:
        return self.gp_hyperparams()[1]
 
    # hyperparameter related functions, might discard ?? not sure if this is gonna improve performance
    def gp_hyperparams(self, theta = None):
        """ hyperaprameters: (lengthscales, amplitudes, noise) """
        theta = self.theta if theta is None else theta
        m = len(self.dims_)
        
        lengthscales = np.exp(theta[:m])
        amplitudes = np.exp(theta[m:2 * m])
        noise = np.exp(theta[-1]) + self.jitter
        return lengthscales, amplitudes, noise
 
    def initial_hyperparams(self, X, Yc):
        """ initial hyperparameters for the GP, derived from data  """
        blocks = [X[:, sl] for sl in self.slices()]
        m = len(blocks)
 
        lengthscales = np.array([self.median_distance(b) for b in blocks])
        amplitude = max(float(Yc.var(axis = 0).mean()), 1e-8) / m
        self.block_lengthscales_init_ = lengthscales
 
        theta0 = np.log(np.concatenate([lengthscales,
                                        np.full(m, amplitude),
                                        [self.noise_frac * amplitude]]))
        
        lengthscale_bounds = [(t - 4.0, t + 4.0) for t in theta0[:m]]
        amplitude_bounds = [(np.log(1e-10), np.log(1e4))] * m
        noise_bounds = [(np.log(1e-12), np.log(1e2))]
        bounds = lengthscale_bounds + amplitude_bounds + noise_bounds
        return theta0, bounds
 
    @staticmethod
    def median_distance(B) -> float:
        """ median-heuristic lengthscale: sqrt(median pairwise sq distance / 2) """
        if len(B) < 2:
            return 1.0
        pairwise_sq = ((B ** 2).sum(1)[:, None] + (B ** 2).sum(1)[None, :] - 2 * B @ B.T)
        median = np.median(pairwise_sq[np.triu_indices(len(B), k = 1)])
        ret = float(np.sqrt(max(median, 1e-12) / 2.0))
        return ret
 

    # kernel functions
    def calc_K(self, A, B, theta = None) -> np.ndarray:
        """
        computes the kernel matrix b/w two sets of reduced inputs.
        adapted from kernel_mat in AdditiveMultiOutputGP
        """
        lengthscales, amplitudes, _ = self.gp_hyperparams(theta)
        K = np.zeros((len(A), len(B)))
        for sl, ls, amp in zip(self.slices(), lengthscales, amplitudes):
            K += self.block_kernel(A[:, sl], B[:, sl], ls, amp)
        return K
 
    def block_kernel(self, A, B, lengthscale, amplitude) -> np.ndarray:
        """
        computes one block's contribution: RBF or linear, scaled by amplitude
        adapted from RBF in kernels.py in LLM-Pert
        """
        if self.kernel == "linear":
            K_one = amplitude * (A @ B.T)
            return K_one

        A = A / lengthscale
        B = B / lengthscale
        s_a = (A ** 2).sum(1)[:, None]
        s_b = (B ** 2).sum(1)[None, :]
        
        sq = s_a + s_b - 2 * A @ B.T
        K_one = amplitude * np.exp(-0.5 * np.maximum(sq, 0))
        return K_one


    def slices(self):
        """ block locations, used by calc_K to split A and B into blocks """
        start = 0
        for width in self.dims_:
            yield slice(start, start + width)
            start += width
 
    def negative_log_marginal_likelihood(self, theta, X, Yc) -> float:
        """ obj for the optimiser """
        K = self.calc_K(X, X, theta) + self.gp_hyperparams(theta)[2] * np.eye(len(X))
        try:
            factor = cho_factor(K, lower = True)
        except np.linalg.LinAlgError:
            return 1e12
        quadratic = np.einsum("ng,ng->", Yc, cho_solve(factor, Yc))
        log_k = 2.0 * np.log(np.diag(factor[0])).sum()
        
        term1 = 0.5 * quadratic
        term2 = 0.5 * Yc.shape[1] * log_k
        return float(term1 + term2) 
 
    # dimensionality reduction helpers 
    def fit_reduce(self, X) -> np.ndarray:
        """
        fits + applies a StandardScaler + PCA per block on a training fold
        sets dims - the block widths after reduction, which the kernel indexes with.
        """
        widths = list(self.feature_dims) if self.feature_dims else [X.shape[1]]
 
        if self.reduce_dim is None:
            self.reducer_ = None
            self.dims_ = widths
            self.retained_variance_ = None
            return X
 
        blocks, start = [], 0
        for i, width in enumerate(widths):
            n_components = min(self.reduce_dim, width, len(X) - 1)
            blocks.append((
                f"block{i}",
                make_pipeline(StandardScaler(with_std = self.with_std),
                              PCA(n_components = n_components, 
                                random_state = 17)),
                list(range(start, start + width)),
            ))
            start += width
 
        self.reducer_ = ColumnTransformer(blocks)
        reduced = self.reducer_.fit_transform(X)
 
        fitted = [pipe.named_steps["pca"] for _, pipe, _ in self.reducer_.transformers_]
        self.dims_ = [pca.n_components_ for pca in fitted]
        self.retained_variance_ = [float(pca.explained_variance_ratio_.sum())
                                   for pca in fitted]
        return reduced
 
    def reduce(self, X) -> np.ndarray:
        return X if self.reducer_ is None else self.reducer_.transform(X)
    
    @property
    def block_importance_(self) -> np.ndarray:
        # Normalised share of fitted kernel amplitude per block.
        amplitudes = self.block_amplitudes
        return amplitudes / amplitudes.sum()
 



def _patch_sklearn_for_mbpls() -> None:
    import functools
    import sklearn.utils as _sk_utils
    import sklearn.utils.validation as _sk_validation

    original = _sk_validation.check_array
    if getattr(original, "_mbpls_patched", False):
        return

    @functools.wraps(original)
    def check_array(*args, **kwargs):
        if "force_all_finite" in kwargs:
            kwargs["ensure_all_finite"] = kwargs.pop("force_all_finite")
        return original(*args, **kwargs)

    check_array._mbpls_patched = True
    _sk_validation.check_array = check_array
    _sk_utils.check_array = check_array


_patch_sklearn_for_mbpls()
from mbpls.mbpls import MBPLS


@dataclass
class SeqOrthPLS:
    feature_dims: tuple
    n_components: tuple = (2, 2)          
    block_order: tuple | None = None     
    with_std: bool = True

    def __post_init__(self):
        if len(self.n_components) != len(self.feature_dims):
            raise ValueError(
                f"n_components has {len(self.n_components)} entries but "
                f"feature_dims describes {len(self.feature_dims)} blocks"
            )

    def _split(self, X):
        blocks, start = [], 0
        for width in self.feature_dims:
            blocks.append(X[:, start:start + width])
            start += width
        return blocks

    def fit(self, X, Y):
        blocks_raw = self._split(X)
        self.order_ = list(self.block_order) if self.block_order is not None else list(range(len(blocks_raw)))

        self.scalers_ = [StandardScaler(with_std=self.with_std).fit(b) for b in blocks_raw]
        blocks = [s.transform(b) for s, b in zip(self.scalers_, blocks_raw)]

        self.pls_models_ = {}
        self.ortho_weights_ = {}
        self.var_reduction_ = {}                  
        Y_resid = Y.astype(float).copy()
        T_cum = None

        for b in self.order_:
            Xb = blocks[b]
            if T_cum is None:
                Xb_orth = Xb
                self.ortho_weights_[b] = None
            else:
                W = np.linalg.pinv(T_cum) @ Xb
                Xb_orth = Xb - T_cum @ W
                self.ortho_weights_[b] = W

            nc = min(self.n_components[b], Xb_orth.shape[1], len(X) - 1)
            pls = PLSRegression(n_components=nc)
            pls.fit(Xb_orth, Y_resid)
            self.pls_models_[b] = pls

            var_before = Y_resid.var(axis=0).sum()      
            Y_resid = Y_resid - pls.predict(Xb_orth)
            var_after = Y_resid.var(axis=0).sum()      
            self.var_reduction_[b] = max(var_before - var_after, 0.0) 

            Tb = pls.x_scores_
            T_cum = Tb if T_cum is None else np.hstack([T_cum, Tb])

        return self

    @property
    def block_importance_(self) -> np.ndarray:
        """ normalized per-block share of Y variance explained, indexed by the feature dims for the table later"""
        reductions = np.array([self.var_reduction_[b] for b in range(len(self.feature_dims))])
        total = reductions.sum()
        return reductions / total if total > 1e-12 else reductions

    def predict(self, X):
        blocks_raw = self._split(X)
        blocks = [s.transform(b) for s, b in zip(self.scalers_, blocks_raw)]

        total_pred = None
        T_cum = None
        for b in self.order_:
            Xb = blocks[b]
            W = self.ortho_weights_[b]
            Xb_orth = Xb if W is None else Xb - T_cum @ W

            pls = self.pls_models_[b]
            pred_b = pls.predict(Xb_orth)
            total_pred = pred_b if total_pred is None else total_pred + pred_b

            Tb = pls.transform(Xb_orth)
            T_cum = Tb if T_cum is None else np.hstack([T_cum, Tb])

        return total_pred


@dataclass
class MultiBlockPLS:
    """ wrapper around mbpls.MBPLS. n_components"""
    feature_dims: tuple
    n_components: int = 3
    with_std: bool = True

    def _split(self, X):
        blocks, start = [], 0
        for width in self.feature_dims:
            blocks.append(X[:, start:start + width])
            start += width
        return blocks

    def fit(self, X, Y):
        blocks = self._split(X)
        nc = min(self.n_components, len(X) - 1, min(self.feature_dims))
        self.model_ = MBPLS(n_components=nc, standardize=self.with_std)
        self.model_.fit(X=blocks, Y=Y)
        return self

    def predict(self, X):
        blocks = self._split(X)
        preds = self.model_.predict(X=blocks)
        return np.asarray(preds)

    @property
    def block_importance_(self) -> np.ndarray:
        weights = np.abs(np.asarray(self.model_.A_corrected_))
        if weights.ndim > 1:
            weights = weights.sum(axis=1)
        return weights / weights.sum()
    
    @property
    def explained_variance_(self) -> dict:
        return {
            "y": float(np.asarray(self.model_.explained_var_y_).sum()),
            "x_blocks": np.asarray(self.model_.explained_var_xblocks_).tolist(),
        }

# model dictionary 
MODELS = {
    "control": ControlBaseline,
    "mean": MeanBaseline,
    "random": RandomBaseline,
    "ridge": Ridge,
    "knn": KNN,
    "mlp": MLP,
    "gp": partial(GaussianProcess, prior_mean = "train_mean"),
    
    # fusion embeddings 
    "gp_additive": partial(GaussianProcess, prior_mean = "train_mean", reduce_dim = 50),

}


def register_grid(prefix: str, cls, grid: dict) -> list[str]:
    """Register one model variant per hyperparameter combination.
    -> adds hidden layer size + alpha value to MODELS
       and returns their names.
    """
    from itertools import product
 
    def short(key, value):
        head = "".join(w[0] for w in key.split("_"))
        if value is None:
            value = "full"
        if isinstance(value, tuple):
            value = "x".join(str(v) for v in value)
        elif isinstance(value, float):
            value = f"{value:g}"                       # 10.0 -> "10", 0.001 -> "0.001"
        return f"{head}{value}".replace(".", "p").replace("-", "m")
 
    names = []
    keys = list(grid)
    for combo in product(*(grid[k] for k in keys)):
        kwargs = dict(zip(keys, combo))
        name = "_".join([prefix] + [short(k, v) for k, v in kwargs.items()])
        MODELS[name] = partial(cls, **kwargs)
        names.append(name)
    return names
 

# # hyperparameter testing for MLP, KNN, GP! report in appendix? 
# MLP_GRID = register_grid("mlp", MLP, {
#     "hidden_layer_sizes": [(64,), (528, 256,), (1024,)],
#     "alpha": [0, 0.001, 0.1, 1.0, 10.0],
# })

# MLP_GRID = register_grid("mlp", MLP, {
#     "hidden_layer_sizes": [(528, 256), (1024, 512), (512, 256, 128)],
#     "alpha": [0.0, 0.001, 0.0001, 10.0],
# })

# MLP_GRID = register_grid("mlp", MLP, {
#     "hidden_layer_sizes": [(1024, 512), (2048, 1024), (1024, 512, 256)],
#     "alpha": [0.0, 0.0001, 0.00001],
# })

# MLP_GRID = register_grid("mlp", MLP, {
#     "hidden_layer_sizes": [(2048, 1024), (2048, 1024, 512)],
#     "alpha": [0.0, 0.00001],
# })

# KNN_GRID = register_grid("knn", KNN, {
#     "k": [1, 3, 5, 10, 20, 40],
#     "metric": ["cosine", "euclidean"],
# })

# KNN_GRID = register_grid("knn", KNN, {
#     "k": [40, 60, 80, 100, 120, 149],
#     "metric": ["cosine"],
# })

# GP_GRID = register_grid("gp", GaussianProcess, {
#     "reduce_dim": [10, 25, 50, 100, 148, None],
#     "prior_mean": ["train_mean"],
# })
 
# 2-block fusion configs: sand_x_semantic, sand_x_combined_semantic(*), str_x_semantic, sand_x_rdkit
SOPLS_GRID_2BLOCK = register_grid("so_pls2", SeqOrthPLS, {
    "n_components": [(1, 1), (2, 2), (3, 3), (5, 5), (2, 5), (5, 2), (8, 8)],
})

# 3-block: sand_x_combined_semantic
SOPLS_GRID_3BLOCK = register_grid("so_pls3", SeqOrthPLS, {
    "n_components": [(2, 2, 2), (3, 3, 3), (5, 5, 5), (2, 3, 5)],
})

MBPLS_GRID = register_grid("mb_pls", MultiBlockPLS, {
    "n_components": [1, 2, 3, 5, 8, 12],
})

def build(name: str, **kwargs):
    if name not in MODELS:
        raise ValueError(f"unknown model {name!r}, have {sorted(MODELS)}")
    factory = MODELS[name]
    cls = factory.func if isinstance(factory, partial) else factory
    accepted = ({f.name for f in dataclasses.fields(cls)}
                if dataclasses.is_dataclass(cls) else set())
    return factory(**{k: v for k, v in kwargs.items() if k in accepted})