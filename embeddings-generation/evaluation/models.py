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
    k: int = 5
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

    hidden_layer_sizes: tuple = (256,) #TODO: ask if there's a better layer size to test
    alpha: float = 10   #TODO: test a few different ones , tried out 1                      
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
    prior_mean: str = "train_mean"
    feature_dims: tuple | None = None
    reduce_dim: int | None = 50 # TODO: ask if this should be increased
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
 


# model dictionary 
MODELS = {
    "control": ControlBaseline,
    "mean": MeanBaseline,
    "random": RandomBaseline,
    "ridge": Ridge,
    "knn": KNN,
    "mlp": MLP,
    
    "gp": partial(GaussianProcess, prior_mean = "train_mean"),
    # "gp_zero": partial(GaussianProcess, prior_mean = "zero"),
    # "gp_linear": partial(GaussianProcess, kernel = "linear"),
}

def build(name: str, **kwargs):
    if name not in MODELS:
        raise ValueError(f"unknown model {name!r}, have {sorted(MODELS)}")
    factory = MODELS[name]
    cls = factory.func if isinstance(factory, partial) else factory
    accepted = ({f.name for f in dataclasses.fields(cls)}
                if dataclasses.is_dataclass(cls) else set())
    return factory(**{k: v for k, v in kwargs.items() if k in accepted})