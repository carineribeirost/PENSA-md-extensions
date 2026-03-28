import numpy as np
import deeptime
import matplotlib.pyplot as plt
from collections import namedtuple

from pensa.clusters.clone import CLoNe
from pensa.clusters.ticc.TICC_solver import TICC


# --- Named return types ---

ClusterResult = namedtuple(
    'ClusterResult',
    ['cidx', 'wss', 'centroids']
)

CombinedClusterResult = namedtuple(
    'CombinedClusterResult',
    ['cidx', 'cond', 'oidx', 'wss', 'centroids']
)


# --- Internal helpers ---

def _run_kmeans(data, num_clusters, max_iter):
    clusters = deeptime.clustering.KMeans(
        n_clusters=num_clusters, max_iter=max_iter).fit(data)
    clusters = clusters.fetch_model()
    cidx = clusters.transform(data)
    return cidx, clusters.inertia, clusters.cluster_centers


def _run_rspace(data, min_dist):
    clusters = deeptime.clustering.RegularSpace(dmin=min_dist).fit(data)
    clusters = clusters.fetch_model()
    cidx = clusters.transform(data)
    return cidx, clusters.inertia, clusters.cluster_centers


def _run_clone(data, pdc=4, n_resize=4, filt=0.1, verbose=False):
    """
    Run CLoNe density-based clustering.

    Parameters
    ----------
    data : float array, shape (n_frames, n_features)
    pdc : float
        Percentile distance cutoff.
    n_resize : int
        Resize parameter.
    filt : float
        Outlier fraction threshold.
    verbose : bool

    Returns
    -------
    cidx : int array
        Cluster index per frame (-1 = outlier).
    wss : float
        Within-cluster sum of squares (np.nan — CLoNe does not compute WSS).
    centroids : float array
        Cluster centres.
    """
    model = CLoNe(pdc=pdc, n_resize=n_resize, filt=filt, verbose=verbose)
    model.fit(data)
    cidx = model.labels_
    centroids = model.centers
    # WSS is not defined for density-based clustering; return nan
    wss = float('nan')
    return cidx, wss, centroids


def _run_ticc(data, num_clusters=5, window_size=10, lambda_parameter=11e-2,
              beta=400, max_iter=1000, compute_bic=False,
              pca_variance_threshold=0.90, pca_n_components=None):
    """
    Run TICC temporal clustering.

    Parameters
    ----------
    data : float array, shape (n_frames, n_features)
    num_clusters : int
    window_size : int
    lambda_parameter : float
        Sparsity parameter.
    beta : float
        Temporal consistency (switch penalty).
    max_iter : int
    compute_bic : bool
    pca_variance_threshold : float or None
        Reduce data by PCA until this fraction of variance is explained.
        Applied automatically when n_features > 50. Default: 0.90.
        Set to None to skip PCA entirely.
    pca_n_components : int or None
        If set, overrides pca_variance_threshold and uses exactly this many
        PCA components. Default: None.

    Returns
    -------
    cidx : int array
        Cluster index per frame.
    wss : float
        np.nan — TICC does not compute WSS.
    centroids : None
        TICC does not expose cluster centres.
    bic : float or None
        BIC value, returned only if compute_bic=True.
    """
    n_features = data.shape[1]
    do_pca = (pca_n_components is not None) or \
             (pca_variance_threshold is not None and n_features > 50)

    if do_pca:
        from sklearn.decomposition import PCA as SklearnPCA
        if pca_n_components is not None:
            # explicit override
            pca = SklearnPCA(n_components=pca_n_components)
            data = pca.fit_transform(data)
        else:
            # auto: find number of components for variance threshold
            pca_full = SklearnPCA()
            pca_full.fit(data)
            cumvar = np.cumsum(pca_full.explained_variance_ratio_)
            n_comp = int(np.searchsorted(cumvar, pca_variance_threshold) + 1)
            pca = SklearnPCA(n_components=n_comp)
            data = pca.fit_transform(data)
        explained = pca.explained_variance_ratio_.sum()
        print(f"  [ticc] PCA: {data.shape[1]} components "
              f"({100*explained:.1f}% variance explained)")

    model = TICC(
        window_size=window_size,
        number_of_clusters=num_clusters,
        lambda_parameter=lambda_parameter,
        beta=beta,
        maxIters=max_iter,
        compute_BIC=compute_bic,
    )
    result = model.fit_array(data)
    if compute_bic:
        clustered_points, _, bic = result
    else:
        clustered_points, _ = result
        bic = None
    cidx = np.array(clustered_points, dtype=int)
    wss = float('nan')
    centroids = None
    if compute_bic:
        return cidx, wss, centroids, bic
    return cidx, wss, centroids


# --- Public API ---


def obtain_clusters(data, algorithm='kmeans',
                    num_clusters=2, min_dist=12, max_iter=100,
                    pdc=4, n_resize=4, filt=0.1,
                    window_size=10, lambda_parameter=11e-2, beta=400,
                    compute_bic=False,
                    pca_variance_threshold=0.90, pca_n_components=None,
                    plot=True, saveas=None):
    """
    Clusters the provided data.

    Parameters
    ----------
    data : float array
        Trajectory data. Format: [frames, frame_data]
    algorithm : str
        Clustering algorithm. Options: kmeans, rspace, clone, ticc.
        Default: kmeans
    num_clusters : int, optional
        Number of clusters (kmeans, ticc). Default: 2.
    min_dist : float, optional
        Minimum distance for regspace clustering. Default: 12.
    max_iter : int, optional
        Maximum number of iterations (kmeans, ticc). Default: 100.
    pdc : float, optional
        Percentile distance cutoff (CLoNe). Default: 4.
    n_resize : int, optional
        Resize parameter (CLoNe). Default: 4.
    filt : float, optional
        Outlier threshold (CLoNe). Default: 0.1.
    window_size : int, optional
        Sliding window size (TICC). Default: 10.
    lambda_parameter : float, optional
        Sparsity parameter (TICC). Default: 11e-2.
    beta : float, optional
        Temporal consistency parameter (TICC). Default: 400.
    compute_bic : bool, optional
        If True, TICC also returns BIC. Default: False.
    pca_n_components : int or None, optional
        PCA pre-reduction before TICC. Strongly recommended for high-dimensional
        features (e.g. bb-torsions). Default: None.
    plot : bool, optional
        Create a bar chart of cluster populations. Default: True.
    saveas : str, optional
        File name to save the plot.

    Returns
    -------
    ClusterResult
        .cidx       : int array — cluster index per frame (-1 = outlier for CLoNe)
        .wss        : float — within-sum-of-squares (nan for CLoNe/TICC)
        .centroids  : float array or None — cluster centres
    """
    if algorithm == 'kmeans':
        cidx, wss, centroids = _run_kmeans(data, num_clusters, max_iter)
    elif algorithm == 'rspace':
        cidx, wss, centroids = _run_rspace(data, min_dist)
    elif algorithm == 'clone':
        cidx, wss, centroids = _run_clone(data, pdc=pdc, n_resize=n_resize, filt=filt)
    elif algorithm == 'ticc':
        out = _run_ticc(data, num_clusters=num_clusters, window_size=window_size,
                        lambda_parameter=lambda_parameter, beta=beta, max_iter=max_iter,
                        compute_bic=compute_bic,
                        pca_variance_threshold=pca_variance_threshold,
                        pca_n_components=pca_n_components)
        if compute_bic:
            cidx, wss, centroids, _ = out
        else:
            cidx, wss, centroids = out
    else:
        raise ValueError(f"Unknown algorithm '{algorithm}'. "
                         f"Choose from: kmeans, rspace, clone, ticc.")

    # Plot cluster populations
    if plot:
        fig, ax = plt.subplots(1, 1, figsize=[4, 3], dpi=300)
        valid = cidx[cidx >= 0]
        c, nc = np.unique(valid, return_counts=True)
        ax.bar(c, nc)
        ax.set_xlabel('cluster')
        ax.set_ylabel('population')
        fig.tight_layout()
        if saveas is not None:
            fig.savefig(saveas, dpi=300)

    return ClusterResult(cidx=cidx, wss=wss, centroids=centroids)


def obtain_combined_clusters(data_a, data_b, label_a='Sim A', label_b='Sim B', start_frame=0,
                              algorithm='kmeans', num_clusters=2, min_dist=12, max_iter=100,
                              pdc=4, n_resize=4, filt=0.1,
                              window_size=10, lambda_parameter=11e-2, beta=400,
                              compute_bic=False,
                              pca_variance_threshold=0.90, pca_n_components=None,
                              plot=True, saveas=None):
    """
    Clusters a combination of two data sets.

    Parameters
    ----------
    data_a : float array
        Trajectory data [frames, frame_data]
    data_b : float array
        Trajectory data [frames, frame_data]
    label_a : str, optional
        Label for the plot. Default: Sim A.
    label_b : str, optional
        Label for the plot. Default: Sim B.
    start_frame : int
        Frame from which the clustering data starts. Default: 0.
    algorithm : str
        Clustering algorithm. Options: kmeans, rspace, clone, ticc.
        Default: kmeans
    num_clusters : int, optional
        Number of clusters (kmeans, ticc). Default: 2.
    min_dist : float, optional
        Minimum distance for regspace clustering. Default: 12.
    max_iter : int, optional
        Maximum iterations (kmeans, ticc). Default: 100.
    pdc : float, optional
        Percentile distance cutoff (CLoNe). Default: 4.
    n_resize : int, optional
        Resize parameter (CLoNe). Default: 4.
    filt : float, optional
        Outlier threshold (CLoNe). Default: 0.1.
    window_size : int, optional
        Sliding window size (TICC). Default: 10.
    lambda_parameter : float, optional
        Sparsity parameter (TICC). Default: 11e-2.
    beta : float, optional
        Temporal consistency (TICC). Default: 400.
    compute_bic : bool, optional
        If True, TICC also returns BIC. Default: False.
    plot : bool, optional
        Create a bar chart of cluster populations. Default: True.
    saveas : str, optional
        File name to save the plot.

    Returns
    -------
    CombinedClusterResult
        .cidx       : int array — cluster index per frame
        .cond       : int array — 0 = data_a, 1 = data_b
        .oidx       : int array — original frame index per point
        .wss        : float — WSS (nan for CLoNe/TICC)
        .centroids  : float array or None
    """
    data = np.concatenate([data_a, data_b], 0)
    cond = np.concatenate([np.zeros(len(data_a)), np.ones(len(data_b))])
    oidx = np.concatenate(
        [np.arange(len(data_a)) + start_frame,
         np.arange(len(data_b)) + start_frame])

    if algorithm == 'kmeans':
        cidx, wss, centroids = _run_kmeans(data, num_clusters, max_iter)
    elif algorithm == 'rspace':
        cidx, wss, centroids = _run_rspace(data, min_dist)
    elif algorithm == 'clone':
        cidx, wss, centroids = _run_clone(data, pdc=pdc, n_resize=n_resize, filt=filt)
    elif algorithm == 'ticc':
        out = _run_ticc(data, num_clusters=num_clusters, window_size=window_size,
                        lambda_parameter=lambda_parameter, beta=beta, max_iter=max_iter,
                        compute_bic=compute_bic,
                        pca_variance_threshold=pca_variance_threshold,
                        pca_n_components=pca_n_components)
        if compute_bic:
            cidx, wss, centroids, _ = out
        else:
            cidx, wss, centroids = out
        # TICC only labels training frames (drops last window_size frames).
        # Trim cond/oidx to match.
        cond = cond[:len(cidx)]
        oidx = oidx[:len(cidx)]
    else:
        raise ValueError(f"Unknown algorithm '{algorithm}'. "
                         f"Choose from: kmeans, rspace, clone, ticc.")

    # Plot cluster populations split by condition
    if plot:
        fig, ax = plt.subplots(1, 1, figsize=[4, 3], dpi=300)
        c, _ = np.unique(cidx[cidx >= 0], return_counts=True)
        ca, nca = np.unique(cidx[cond == 0], return_counts=True)
        cb, ncb = np.unique(cidx[cond == 1], return_counts=True)
        # Filter out outliers (-1) from bars
        ca_mask = ca >= 0
        cb_mask = cb >= 0
        ax.bar(ca[ca_mask] - 0.15, nca[ca_mask], 0.3, label=label_a)
        ax.bar(cb[cb_mask] + 0.15, ncb[cb_mask], 0.3, label=label_b)
        ax.legend()
        ax.set_xticks(c)
        ax.set_xlabel('cluster')
        ax.set_ylabel('population')
        fig.tight_layout()
        if saveas is not None:
            fig.savefig(saveas, dpi=300)

    return CombinedClusterResult(cidx=cidx, cond=cond, oidx=oidx, wss=wss, centroids=centroids)


def obtain_mult_combined_clusters(data, start_frame=0, algorithm='kmeans',
                                  num_clusters=2, min_dist=12, max_iter=100,
                                  pdc=4, n_resize=4, filt=0.1,
                                  window_size=10, lambda_parameter=11e-2, beta=400,
                                  compute_bic=False,
                                  pca_variance_threshold=0.90, pca_n_components=None,
                                  plot=True, saveas=None, labels=None, colors=None):
    """
    Clusters a combination of multiple data sets.

    Parameters
    ----------
    data : list of float arrays
        Trajectory data [frames, frame_data]
    start_frame : int
        Frame from which the clustering data starts. Default: 0.
    algorithm : str
        Clustering algorithm. Options: kmeans, rspace, clone, ticc.
        Default: kmeans
    num_clusters : int, optional
        Number of clusters (kmeans, ticc). Default: 2.
    min_dist : float, optional
        Minimum distance for regspace clustering. Default: 12.
    max_iter : int, optional
        Maximum iterations (kmeans, ticc). Default: 100.
    pdc : float, optional
        Percentile distance cutoff (CLoNe). Default: 4.
    n_resize : int, optional
        Resize parameter (CLoNe). Default: 4.
    filt : float, optional
        Outlier threshold (CLoNe). Default: 0.1.
    window_size : int, optional
        Sliding window size (TICC). Default: 10.
    lambda_parameter : float, optional
        Sparsity parameter (TICC). Default: 11e-2.
    beta : float, optional
        Temporal consistency (TICC). Default: 400.
    compute_bic : bool, optional
        If True, TICC also returns BIC. Default: False.
    plot : bool, optional
        Create a bar chart of cluster populations. Default: True.
    saveas : str, optional
        File name to save the plot.
    labels : list of str, optional
        Labels for the plot.
    colors : list of str, optional
        Colors for the plot.

    Returns
    -------
    CombinedClusterResult
        .cidx       : int array
        .cond       : int array — trajectory index per frame
        .oidx       : int array — original frame index
        .wss        : float
        .centroids  : float array or None
    """
    if labels is not None:
        assert len(labels) == len(data)
    else:
        labels = [None for _ in range(len(data))]
    if colors is not None:
        assert len(colors) == len(data)
    else:
        colors = ['C%i' % num for num in range(len(data))]

    num_frames = [len(d) for d in data]
    num_traj = len(data)
    combined = np.concatenate(data, 0)
    cond = np.concatenate([i * np.ones(num_frames[i], dtype=int) for i in range(num_traj)])
    oidx = np.concatenate([np.arange(num_frames[i]) + start_frame for i in range(num_traj)])

    if algorithm == 'kmeans':
        cidx, wss, centroids = _run_kmeans(combined, num_clusters, max_iter)
    elif algorithm == 'rspace':
        cidx, wss, centroids = _run_rspace(combined, min_dist)
    elif algorithm == 'clone':
        cidx, wss, centroids = _run_clone(combined, pdc=pdc, n_resize=n_resize, filt=filt)
    elif algorithm == 'ticc':
        out = _run_ticc(combined, num_clusters=num_clusters, window_size=window_size,
                        lambda_parameter=lambda_parameter, beta=beta, max_iter=max_iter,
                        compute_bic=compute_bic,
                        pca_variance_threshold=pca_variance_threshold,
                        pca_n_components=pca_n_components)
        if compute_bic:
            cidx, wss, centroids, _ = out
        else:
            cidx, wss, centroids = out
        # TICC only labels training frames (drops last window_size frames).
        # Trim cond/oidx to match.
        cond = cond[:len(cidx)]
        oidx = oidx[:len(cidx)]
    else:
        raise ValueError(f"Unknown algorithm '{algorithm}'. "
                         f"Choose from: kmeans, rspace, clone, ticc.")

    if plot:
        fig, ax = plt.subplots(1, 1, figsize=[4, 3], dpi=300)
        c, _ = np.unique(cidx[cidx >= 0], return_counts=True)
        for n in range(num_traj):
            ca, nca = np.unique(cidx[cond == n], return_counts=True)
            ca_mask = ca >= 0
            bar_pos = ca[ca_mask] - 0.2 + n * 0.4 / num_traj
            bar_width = 0.4 / num_traj
            ax.bar(bar_pos, nca[ca_mask], bar_width, label=labels[n], color=colors[n])
        ax.legend()
        ax.set_xticks(c)
        ax.set_xlabel('cluster')
        ax.set_ylabel('population')
        fig.tight_layout()
        if saveas is not None:
            fig.savefig(saveas, dpi=300)

    return CombinedClusterResult(cidx=cidx, cond=cond, oidx=oidx, wss=wss, centroids=centroids)


def find_closest_frames(data, points):
    """
    Finds the frames in a timeseries that are closest to given points.

    Parameters
    ----------
    data : float array
        Trajectory data [frames, frame_data]
    points : list of float arrays
        Points to find the closest frames for.

    Returns
    -------
    frames : list of int
    distances : list of float
    """
    frames, distances = [], []
    for p in points:
        diff = data - p
        dist = np.sqrt(np.sum(np.square(diff), axis=1))
        frames.append(np.argmin(dist))
        distances.append(np.min(dist))
    return frames, distances
