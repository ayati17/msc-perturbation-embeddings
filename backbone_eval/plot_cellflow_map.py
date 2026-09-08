from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
import umap

from backbone_eval import config as cfg
from backbone_eval.data import load_adata, compound_folds, split


# ============================================================
# SETTINGS
# ============================================================

# ecfp_raw/counts sand/z  text_semantic/f2llm/curated  
SPEC = "fusion/pca_concat50/sand_x_combined_semantic"
FOLD = 0
CONDITION_ID = None   # If None, the first condition in conditions.parquet is used.
N_NEIGHBORS = 15
MIN_DIST = 0.9
RANDOM_STATE = 42

# Output
OUTPUT = "cellflow_umap_example.png"

# Only visualize conditions where the number of predicted
# cells is reasonably comparable to the number of observed cells.
MAX_PRED_TO_OBS_RATIO = 10

# Optional: require at least this many observed cells
MIN_OBSERVED_CELLS = 200


# ============================================================
# LOAD DATA
# ============================================================

print("Loading SciPlex data...")
adata = load_adata()

fold_map = compound_folds(adata)

print(f"Loading fold {FOLD}...")
_, test = split(adata, fold_map[FOLD])


# ============================================================
# LOAD SAVED CELLFLOW PREDICTION
# ============================================================

pred_dir = (
    cfg.RESULTS
    / "predictions"
    / SPEC
    / f"fold_{FOLD}"
)

conditions_path = pred_dir / "conditions.parquet"

conditions = pd.read_parquet(conditions_path)

print(f"\nFound {len(conditions)} conditions.")


# ============================================================
# FIND A CONDITION WITH A GOOD OBSERVED/PREDICTED RATIO
# ============================================================

if CONDITION_ID is not None:

    # --------------------------------------------------------
    # User explicitly requested a condition
    # --------------------------------------------------------

    matches = conditions[
        conditions["condition_id"] == CONDITION_ID
    ]

    if len(matches) == 0:
        raise ValueError(
            f"Could not find condition: {CONDITION_ID}"
        )

    candidate_conditions = matches

else:

    # --------------------------------------------------------
    # Otherwise, cycle through conditions until we find one
    # with a visually reasonable predicted/observed ratio.
    # --------------------------------------------------------

    candidate_conditions = conditions


condition_row = None
predicted = None
observed_condition = None

for idx, candidate in candidate_conditions.iterrows():

    perturbation = candidate["perturbation"]
    log_dose = candidate["log_dose"]
    cell_line = candidate["cell_line"]

    # --------------------------------------------------------
    # Find observed cells for this condition
    # --------------------------------------------------------

    is_this_condition = (
        (test.obs["perturbation"] == perturbation)
        & (test.obs["log_dose"] == log_dose)
        & (test.obs["cell_line"] == cell_line)
        & (~test.obs["is_control"])
    ).values

    n_observed = int(is_this_condition.sum())

    # Skip conditions with too few observed cells
    if n_observed < MIN_OBSERVED_CELLS:
        print(
            f"[skip] {candidate['condition_id']}: "
            f"only {n_observed} observed cells"
        )
        continue

    # --------------------------------------------------------
    # Load CellFlow prediction
    # --------------------------------------------------------

    pred_path = pred_dir / candidate["prediction_file"]

    if not pred_path.exists():
        print(
            f"[skip] {candidate['condition_id']}: "
            f"prediction file not found"
        )
        continue

    candidate_predicted = np.load(
        pred_path,
        allow_pickle=False,
    )

    n_predicted = len(candidate_predicted)

    # --------------------------------------------------------
    # Calculate ratio
    # --------------------------------------------------------

    ratio = n_predicted / n_observed

    print(
        f"[check] {candidate['condition_id']}: "
        f"observed={n_observed:,} | "
        f"predicted={n_predicted:,} | "
        f"ratio={ratio:.2f}x"
    )

    # --------------------------------------------------------
    # Accept the first visually reasonable condition
    # --------------------------------------------------------

    if ratio <= MAX_PRED_TO_OBS_RATIO:

        condition_row = candidate
        predicted = candidate_predicted
        observed_condition = test[is_this_condition].copy()

        print()
        print(
            f"[selected] {candidate['condition_id']}"
        )
        print(
            f"[selected] observed={n_observed:,} | "
            f"predicted={n_predicted:,} | "
            f"ratio={ratio:.2f}x"
        )

        break


# ============================================================
# MAKE SURE WE FOUND SOMETHING
# ============================================================

if condition_row is None:

    raise RuntimeError(
        "\nNo condition satisfied the visualization criteria.\n"
        f"Required:\n"
        f"  observed cells >= {MIN_OBSERVED_CELLS}\n"
        f"  predicted / observed <= "
        f"{MAX_PRED_TO_OBS_RATIO:.1f}x\n\n"
        "Try increasing MAX_PRED_TO_OBS_RATIO."
    )


# ============================================================
# PRINT SELECTED CONDITION
# ============================================================

condition_id = condition_row["condition_id"]

print("\nSelected condition:")
print(condition_row)

print(
    f"\nPrediction shape: {predicted.shape}"
)


# ============================================================
# FIND OBSERVED / CONTROL CELLS FOR THIS CONDITION
# ============================================================

perturbation = condition_row["perturbation"]
log_dose = condition_row["log_dose"]
cell_line = condition_row["cell_line"]

is_control = test.obs["is_control"].values

is_this_condition = (
    (test.obs["perturbation"] == perturbation)
    & (test.obs["log_dose"] == log_dose)
    & (test.obs["cell_line"] == cell_line)
    & (~test.obs["is_control"])
).values

# observed_condition was selected above
controls = test[test.obs["is_control"]].copy()

print(
    f"Observed perturbed cells: {len(observed_condition)}"
)
print(
    f"Control cells: {len(controls)}"
)
print(
    f"Predicted cells: {len(predicted)}"
)
print(
    f"Predicted / observed ratio: "
    f"{len(predicted) / len(observed_condition):.2f}x"
)

# SAMPLING
# rng = np.random.default_rng(17)

# n = len(observed_condition)

# if len(predicted) > n:
#     idx = rng.choice(
#         len(predicted),
#         size=n,
#         replace=False,
#     )
#     predicted_plot = predicted[idx]
# else:
#     predicted_plot = predicted


# ============================================================
# GET PCA COORDINATES
# ============================================================

# These are already in the fold-specific PCA space used by CellFlow.
X_test_pca = test.obsm["X_pca"]

X_condition_pca = observed_condition.obsm["X_pca"]

X_control_pca = controls.obsm["X_pca"]


# ============================================================
# BUILD ONE SHARED UMAP
# ============================================================

# IMPORTANT:
# Actual cells + CellFlow predictions are embedded together.
# This means every panel uses EXACTLY the same UMAP coordinates.

X_all = np.vstack([
    X_test_pca,
    predicted,
])

n_test = len(X_test_pca)

umap_model = umap.UMAP(
    n_neighbors=N_NEIGHBORS,
    min_dist=MIN_DIST,
    metric="euclidean",
    random_state=RANDOM_STATE,
)

X_umap = umap_model.fit_transform(X_all)

# Split back into actual and predicted coordinates
test_umap = X_umap[:n_test]
pred_umap = X_umap[n_test:]


# ============================================================
# COORDINATES FOR DIFFERENT GROUPS
# ============================================================

# Control cells
control_mask = test.obs["is_control"].values

# True perturbed cells for this particular condition
condition_mask = is_this_condition

# All other cells
other_mask = ~(control_mask | condition_mask)

# Dose values
dose_values = test.obs["log_dose"].astype(float).values


# ============================================================
# FIGURE
# ============================================================

fig, axes = plt.subplots(
    1,
    4,
    figsize=(18, 4.5),
)


# ------------------------------------------------------------
# Shared plotting parameters
# ------------------------------------------------------------

POINT_SIZE = 5

xmin = X_umap[:, 0].min()
xmax = X_umap[:, 0].max()
ymin = X_umap[:, 1].min()
ymax = X_umap[:, 1].max()

xpad = 0.03 * (xmax - xmin)
ypad = 0.03 * (ymax - ymin)

xlim = (xmin - xpad, xmax + xpad)
ylim = (ymin - ypad, ymax + ypad)


# ============================================================
# PANEL A — CONTROL VS OTHER
# ============================================================

ax = axes[0]

# Everything else in gray
ax.scatter(
    test_umap[~control_mask, 0],
    test_umap[~control_mask, 1],
    s=POINT_SIZE,
    c="lightgray",
    alpha=0.5,
    linewidths=0,
)

# Controls in black
ax.scatter(
    test_umap[control_mask, 0],
    test_umap[control_mask, 1],
    s=POINT_SIZE,
    c="black",
    alpha=0.9,
    linewidths=0,
)

ax.set_title("Control")
ax.set_xlim(xlim)
ax.set_ylim(ylim)


legend_elements = [
    Line2D(
        [0], [0],
        marker="o",
        color="w",
        label="control",
        markerfacecolor="black",
        markersize=7,
    ),
    Line2D(
        [0], [0],
        marker="o",
        color="w",
        label="other",
        markerfacecolor="lightgray",
        markersize=7,
    ),
]

ax.legend(
    handles=legend_elements,
    title="Perturbation",
    loc="lower right",
    frameon=False,
)


# ============================================================
# PANEL B — TRUE PERTURBED
# ============================================================

ax = axes[1]

# Everything else
ax.scatter(
    test_umap[other_mask, 0],
    test_umap[other_mask, 1],
    s=POINT_SIZE,
    c="lightgray",
    alpha=0.5,
    linewidths=0,
)

# True cells for this drug/dose/cell line
ax.scatter(
    test_umap[condition_mask, 0],
    test_umap[condition_mask, 1],
    s=POINT_SIZE,
    c="#2f8fb3",
    alpha=0.85,
    linewidths=0,
)

ax.set_title("True perturbed")
ax.set_xlim(xlim)
ax.set_ylim(ylim)

legend_elements = [
    Line2D(
        [0], [0],
        marker="o",
        color="w",
        label="True perturbed",
        markerfacecolor="#2f8fb3",
        markersize=7,
    ),
    Line2D(
        [0], [0],
        marker="o",
        color="w",
        label="other",
        markerfacecolor="lightgray",
        markersize=7,
    ),
]

ax.legend(
    handles=legend_elements,
    title="Perturbation",
    loc="lower right",
    frameon=False,
)


# ============================================================
# PANEL C — CELLFLOW PREDICTIONS
# ============================================================

ax = axes[2]

# Actual cells in gray
ax.scatter(
    test_umap[:, 0],
    test_umap[:, 1],
    s=POINT_SIZE,
    c="lightgray",
    alpha=0.35,
    linewidths=0,
)

# CellFlow predictions
ax.scatter(
    pred_umap[:, 0],
    pred_umap[:, 1],
    s=POINT_SIZE,
    c="#c02a9b",
    alpha=0.85,
    linewidths=0,
)

ax.set_title("CellFlow predictions")
ax.set_xlim(xlim)
ax.set_ylim(ylim)

legend_elements = [
    Line2D(
        [0], [0],
        marker="o",
        color="w",
        label="CellFlow predictions",
        markerfacecolor="#c02a9b",
        markersize=7,
    ),
    Line2D(
        [0], [0],
        marker="o",
        color="w",
        label="other",
        markerfacecolor="lightgray",
        markersize=7,
    ),
]

ax.legend(
    handles=legend_elements,
    title="Perturbation",
    loc="lower right",
    frameon=False,
)


# # ============================================================
# # PANEL D — DOSE
# # ============================================================

# ax = axes[3]

# # Background
# ax.scatter(
#     test_umap[~condition_mask, 0],
#     test_umap[~condition_mask, 1],
#     s=POINT_SIZE,
#     c="lightgray",
#     alpha=0.3,
#     linewidths=0,
# )

# # This treatment, colored by dose
# dose_plot = ax.scatter(
#     test_umap[condition_mask, 0],
#     test_umap[condition_mask, 1],
#     s=POINT_SIZE,
#     c=dose_values[condition_mask],
#     cmap="viridis",
#     alpha=0.9,
#     linewidths=0,
# )

# ax.set_title("Treatment dose")
# ax.set_xlim(xlim)
# ax.set_ylim(ylim)

# cbar = fig.colorbar(
#     dose_plot,
#     ax=ax,
#     fraction=0.046,
#     pad=0.04,
# )

# cbar.set_label(r"$\log_{10}$ dose in M")


# ============================================================
# FORMATTING
# ============================================================

for ax in axes:
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)


fig.suptitle(
    f"{perturbation} treatment on "
    f"{cell_line} cell line",
    fontsize=13,
    y=1.02,
)

plt.tight_layout()

plt.savefig(
    OUTPUT,
    dpi=300,
    bbox_inches="tight",
)

plt.show()

print(f"\nSaved figure to: {OUTPUT}")