# PENSA-md-extensions

Extensions to [PENSA](https://github.com/drorlab/pensa) adding CLoNe and TICC as new clustering methods for MD ensemble comparison.

> This is a fork of drorlab/pensa. All original functionality is preserved. The changes described here are additive.

---

## New clustering methods

PENSA originally provides **k-means** and **r-space** clustering. This fork adds:

| Method | Type | k selection | Temporal sensitivity |
|--------|------|-------------|----------------------|
| **CLoNe** | Density-based | Automatic | No |
| **TICC** | MRF / temporal | Manual | Yes |

### CLoNe (Clustering by Local Density)
Source: [LBM-EPFL/CLoNe](https://github.com/LBM-EPFL/CLoNe)

- Detects clusters based on local density — no need to specify `k` in advance
- Returns a `pdc` (density-distance cutoff) sensitivity plot to guide parameter choice
- Frames that do not fit any cluster are labeled as outliers (`-1`)
- Useful when the number of conformational states is unknown

### TICC (Toeplitz Inverse Covariance-based Clustering)
Source: [davidhallac/TICC](https://github.com/davidhallac/TICC)

- Clusters time-series data by learning a Markov Random Field (MRF) per cluster
- Sensitive to the **temporal order** of frames — identifies when the receptor transitions between states
- Automatically applies PCA before clustering to handle high-dimensional MD features
- Useful when conformational transitions are of interest, not just final state populations

---

## New files

```
pensa/clusters/
├── clone.py                  # CLoNe clustering wrapper
├── cluster_selection.py      # WSS, PDC sensitivity, TICC BIC functions
└── ticc/
    ├── TICC_solver.py        # TICC solver (added fit_array() for numpy input)
    ├── __init__.py
    └── src/
        ├── TICC_helper.py
        ├── admm_solver.py
        └── __init__.py
```

**Modified files:**
- `pensa/clusters/clustering.py` — added `obtain_clusters`, `obtain_combined_clusters`, `obtain_mult_combined_clusters` with CLoNe/TICC support; uses `ClusterResult` and `CombinedClusterResult` namedtuples
- `pensa/clusters/trajectory.py` — fixed bounds check for TICC (which drops last `window_size - 1` frames)
- `pensa/clusters/__init__.py` — updated exports
- `scripts/calculate_combined_clusters.py` — added `--algorithm` argument, algorithm name in all output paths
- `tutorial/4-clustering.sh` — added k-means, CLoNe, and TICC tutorial blocks

---

## Usage

### Via script (recommended)

#### k-means — two-step workflow

k-means requires choosing `k` before the main run. The recommended workflow is:

**Step 1 — generate the WSS elbow plot (no output written)**
```bash
python scripts/calculate_combined_clusters.py \
    --wss --no-write \
    --ref_file_a traj/condition-a_receptor.gro \
    --trj_file_a traj/condition-a_receptor.xtc \
    --ref_file_b traj/condition-b_receptor.gro \
    --trj_file_b traj/condition-b_receptor.xtc \
    --label_a apo --label_b bu72 \
    --out_plots plots/receptor \
    --out_results results/receptor \
    --out_frames_a clusters/receptor_apo \
    --out_frames_b clusters/receptor_bu72 \
    --start_frame 2000 \
    --algorithm kmeans \
    --max_num_clusters 12
```

Inspect `plots/receptor_kmeans_wss_bb-torsions.pdf` and `plots/receptor_kmeans_wss_sc-torsions.pdf`.
Choose the `k` at the elbow of each curve.

**Step 2 — run clustering with chosen k**
```bash
python scripts/calculate_combined_clusters.py \
    --no-wss --write \
    --ref_file_a traj/condition-a_receptor.gro \
    --trj_file_a traj/condition-a_receptor.xtc \
    --ref_file_b traj/condition-b_receptor.gro \
    --trj_file_b traj/condition-b_receptor.xtc \
    --label_a apo --label_b bu72 \
    --out_plots plots/receptor \
    --out_results results/receptor \
    --out_frames_a clusters/receptor_apo \
    --out_frames_b clusters/receptor_bu72 \
    --start_frame 2000 \
    --algorithm kmeans \
    --write_num_clusters 2   # <-- k chosen from WSS plot
```

> Both steps can be combined in a single call (`--wss --write --write_num_clusters 2`) if you already know k from prior experience with the system.

#### CLoNe — two-step workflow

CLoNe selects k automatically, but requires choosing the `pdc` (percentile density cutoff). Same two-step pattern:

**Step 1 — PDC sensitivity sweep**
```bash
python scripts/calculate_combined_clusters.py \
    --no-write --pdc-sensitivity \
    --ref_file_a traj/condition-a_receptor.gro \
    --trj_file_a traj/condition-a_receptor.xtc \
    --ref_file_b traj/condition-b_receptor.gro \
    --trj_file_b traj/condition-b_receptor.xtc \
    --label_a apo --label_b bu72 \
    --out_plots plots/receptor \
    --out_results results/receptor \
    --out_frames_a clusters/receptor_apo \
    --out_frames_b clusters/receptor_bu72 \
    --start_frame 2000 \
    --algorithm clone \
    --n_resize 4 \
    --filt 0.1
```

Inspect `plots/receptor_clone_pdc_sensitivity_bb-torsions.pdf`. Choose the lowest `pdc` value where the number of clusters is stable.

**Step 2 — run clustering with chosen pdc**
```bash
python scripts/calculate_combined_clusters.py \
    --no-pdc-sensitivity --write \
    --ref_file_a traj/condition-a_receptor.gro \
    --trj_file_a traj/condition-a_receptor.xtc \
    --ref_file_b traj/condition-b_receptor.gro \
    --trj_file_b traj/condition-b_receptor.xtc \
    --label_a apo --label_b bu72 \
    --out_plots plots/receptor \
    --out_results results/receptor \
    --out_frames_a clusters/receptor_apo \
    --out_frames_b clusters/receptor_bu72 \
    --start_frame 2000 \
    --algorithm clone \
    --pdc 4.0 \   # <-- pdc chosen from sensitivity plot
    --n_resize 4 \
    --filt 0.1

#### TICC

TICC requires specifying `k` manually. Use `--ticc-bic` to run a BIC sweep over a range of k values and inspect the plot, or rely on prior knowledge of the system.

For bb-torsions with the MOR tutorial dataset, the following parameters are validated (see [TICC parameter notes](#ticc-parameter-notes) below):

```bash
python scripts/calculate_combined_clusters.py \
    --write --no-wss \
    --ref_file_a traj/condition-a_receptor.gro \
    --trj_file_a traj/condition-a_receptor.xtc \
    --ref_file_b traj/condition-b_receptor.gro \
    --trj_file_b traj/condition-b_receptor.xtc \
    --label_a apo --label_b bu72 \
    --out_plots plots/receptor \
    --out_results results/receptor \
    --out_frames_a clusters/receptor_apo \
    --out_frames_b clusters/receptor_bu72 \
    --start_frame 2000 \
    --algorithm ticc \
    --ticc_clusters 2 \
    --window_size 5 \
    --lambda_parameter 0.05 \
    --beta 50 \
    --pca_variance_threshold 0.80 \
    --ticc_max_iter 100
```

### Output naming

All outputs include the algorithm name to avoid overwriting across runs:

| File type | Pattern |
|-----------|---------|
| Cluster plot | `{out_plots}_{algo}_combined-clusters_{ftype}.pdf` |
| WSS elbow plot | `{out_plots}_{algo}_wss_{ftype}.pdf` |
| PDC sensitivity plot | `{out_plots}_{algo}_pdc_sensitivity_{ftype}.pdf` |
| Cluster indices CSV | `{out_results}_{algo}_combined-cluster-indices_{ftype}.csv` |
| Cluster trajectory | `{out_frames}_{algo}_{ftype}_c{N}.xtc` |

### Via Python API

```python
from pensa.clusters import obtain_combined_clusters, ClusterResult, CombinedClusterResult

# CLoNe
result = obtain_combined_clusters(
    data_a, data_b,
    label_a='apo', label_b='bu72',
    algorithm='clone',
    pdc=4.0, n_resize=4, filt=0.1,
    saveas='plots/combined-clusters.pdf'
)
print(result.cidx)   # cluster index per frame
print(result.cond)   # condition per frame (0=A, 1=B)

# TICC (PCA applied automatically when n_features > 50)
result = obtain_combined_clusters(
    data_a, data_b,
    label_a='apo', label_b='bu72',
    algorithm='ticc',
    num_clusters=2,
    window_size=5, lambda_parameter=0.05, beta=50,
    pca_variance_threshold=0.80
)
```

---

## TICC parameter notes

TICC is sensitive to its parameters. Key considerations:

- **`pca_variance_threshold`**: PCA is applied automatically when `n_features > 50`. For bb-torsions (574 features), `0.80` gives ~43 components (workable). `0.90` gives ~97 components and may produce singular covariance matrices.
- **`beta`**: switch penalty between clusters. Values above ~100 tend to collapse all frames into one cluster.
- **`window_size`**: TICC drops the last `window_size - 1` frames from labeling — this is handled automatically.
- **sc-torsions**: with the standard MOR dataset, 80% PCA variance requires ~231 components → stacked feature dimension of 1155 (231 × window_size=5). This makes covariance estimation expensive. Whether this is feasible depends on available memory and compute.

Validated parameters for MOR bb-torsions: `window_size=5, lambda=0.05, beta=50, pca_variance_threshold=0.80`

---

## Dependencies

In addition to PENSA's original dependencies:

```
scikit-learn   # for PCA pre-reduction in TICC
cvxpy          # required by TICC admm_solver
```

Install with:
```bash
pip install scikit-learn cvxpy
```

---

## Compatibility notes

- NumPy >= 2.0: `np.NINF` was removed — fixed in `clone.py` (uses `-np.inf`)
- TICC `fit_array()` method added to `TICC_solver.py` to accept numpy arrays directly (original only reads from CSV files)
