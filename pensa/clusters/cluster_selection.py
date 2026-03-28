import numpy as np
import matplotlib.pyplot as plt

from pensa.clusters import obtain_clusters, obtain_combined_clusters


# --- WSS (k-means / rspace) ---


def wss_over_number_of_clusters(data, algorithm='kmeans',
                                max_iter=100, num_repeats=5,
                                max_num_clusters=12,
                                plot_file=None):
    """
    Calculates the within-sum-of-squares (WSS) for different numbers of clusters,
    averaged over several iterations.

    Parameters
    ----------
    data : float array
        Trajectory data [frames, frame_data]
    algorithm : str
        Options: kmeans, rspace. Default: kmeans
    max_iter : int, optional
        Default: 100.
    num_repeats : int, optional
        Default: 5.
    max_num_clusters : int, optional
        Default: 12.
    plot_file : str, optional
        File path to save the elbow plot.

    Returns
    -------
    all_wss : list of float
    std_wss : list of float
    """
    all_wss, std_wss = [], []
    for nc in range(1, max_num_clusters):
        rep_wss = []
        for _ in range(num_repeats):
            cc = obtain_clusters(data, algorithm=algorithm, max_iter=max_iter,
                                 num_clusters=nc, plot=False)
            rep_wss.append(cc.wss)
        all_wss.append(np.mean(rep_wss))
        std_wss.append(np.std(rep_wss))

    fig, ax = plt.subplots(1, 1, figsize=[4, 3], dpi=300)
    ax.errorbar(np.arange(len(all_wss)) + 1, np.array(all_wss),
                yerr=np.array(std_wss) / np.sqrt(num_repeats))
    ax.set_xlabel('number of clusters')
    ax.set_ylabel('total WSS')
    fig.tight_layout()
    if plot_file:
        fig.savefig(plot_file)

    return all_wss, std_wss


def wss_over_number_of_combined_clusters(data_a, data_b,
                                         label_a='Sim A', label_b='Sim B',
                                         start_frame=0,
                                         algorithm='kmeans', max_iter=100,
                                         num_repeats=5, max_num_clusters=12,
                                         plot_file=None):
    """
    Calculates the WSS for different numbers of clusters on combined data,
    averaged over several iterations.

    Parameters
    ----------
    data_a, data_b : float array
        Trajectory data [frames, frame_data]
    label_a, label_b : str, optional
    start_frame : int, optional
    algorithm : str
        Options: kmeans, rspace. Default: kmeans
    max_iter : int, optional
    num_repeats : int, optional
    max_num_clusters : int, optional
    plot_file : str, optional

    Returns
    -------
    all_wss : list of float
    std_wss : list of float
    """
    all_wss, std_wss = [], []
    for nc in range(1, max_num_clusters):
        rep_wss = []
        for _ in range(num_repeats):
            cc = obtain_combined_clusters(data_a, data_b, label_a, label_b, start_frame,
                                          algorithm=algorithm, max_iter=max_iter,
                                          num_clusters=nc, plot=False)
            rep_wss.append(cc.wss)
        all_wss.append(np.mean(rep_wss))
        std_wss.append(np.std(rep_wss))

    fig, ax = plt.subplots(1, 1, figsize=[4, 3], dpi=300)
    ax.errorbar(np.arange(len(all_wss)) + 1, np.array(all_wss),
                yerr=np.array(std_wss) / np.sqrt(num_repeats))
    ax.set_xlabel('number of clusters')
    ax.set_ylabel('total WSS')
    fig.tight_layout()
    if plot_file:
        fig.savefig(plot_file)

    return all_wss, std_wss


# --- CLoNe pdc sensitivity ---


def pdc_sensitivity_over_clusters(data, pdc_values=None,
                                   n_resize=4, filt=0.1,
                                   plot_file=None):
    """
    Sweeps pdc values for CLoNe on a single dataset and plots
    the number of clusters detected at each pdc.

    Parameters
    ----------
    data : float array
        Trajectory data [frames, frame_data]
    pdc_values : array-like, optional
        pdc values to sweep. Default: 1..20.
    n_resize : int, optional
        CLoNe resize parameter. Default: 4.
    filt : float, optional
        CLoNe outlier threshold. Default: 0.1.
    plot_file : str, optional

    Returns
    -------
    pdc_values : list
    n_clusters : list of int
    """
    if pdc_values is None:
        pdc_values = list(range(1, 21))

    n_clusters = []
    for pdc in pdc_values:
        cc = obtain_clusters(data, algorithm='clone', pdc=pdc,
                             n_resize=n_resize, filt=filt, plot=False)
        k = int(np.sum(np.unique(cc.cidx) >= 0))
        n_clusters.append(k)
        print(f"  pdc={pdc:5.1f}  ->  {k} cluster(s)")

    fig, ax = plt.subplots(1, 1, figsize=[5, 3], dpi=300)
    ax.plot(pdc_values, n_clusters, 'o-')
    ax.set_xlabel('pdc')
    ax.set_ylabel('number of clusters')
    ax.set_title('CLoNe pdc sensitivity')
    fig.tight_layout()
    if plot_file:
        fig.savefig(plot_file)

    return pdc_values, n_clusters


def pdc_sensitivity_for_combined_clusters(data_a, data_b,
                                           label_a='Sim A', label_b='Sim B',
                                           start_frame=0,
                                           pdc_values=None,
                                           n_resize=4, filt=0.1,
                                           plot_file=None):
    """
    Sweeps pdc values for CLoNe on combined data and plots
    the number of clusters detected at each pdc.

    Use this to choose pdc before the main CLoNe run — analogous
    to the WSS elbow plot for k-means.

    Parameters
    ----------
    data_a, data_b : float array
        Trajectory data [frames, frame_data]
    label_a, label_b : str, optional
    start_frame : int, optional
    pdc_values : array-like, optional
        Default: 1..20.
    n_resize : int, optional
    filt : float, optional
    plot_file : str, optional

    Returns
    -------
    pdc_values : list
    n_clusters : list of int
    """
    if pdc_values is None:
        pdc_values = list(range(1, 21))

    n_clusters = []
    for pdc in pdc_values:
        cc = obtain_combined_clusters(data_a, data_b, label_a, label_b, start_frame,
                                      algorithm='clone', pdc=pdc,
                                      n_resize=n_resize, filt=filt, plot=False)
        k = int(np.sum(np.unique(cc.cidx) >= 0))
        n_clusters.append(k)
        print(f"  pdc={pdc:5.1f}  ->  {k} cluster(s)")

    fig, ax = plt.subplots(1, 1, figsize=[5, 3], dpi=300)
    ax.plot(pdc_values, n_clusters, 'o-')
    ax.set_xlabel('pdc')
    ax.set_ylabel('number of clusters')
    ax.set_title('CLoNe pdc sensitivity')
    fig.tight_layout()
    if plot_file:
        fig.savefig(plot_file)

    return pdc_values, n_clusters


# --- TICC BIC ---


def ticc_bic_over_clusters(data, min_clusters=2, max_clusters=10,
                            window_size=10, lambda_parameter=11e-2,
                            beta=400, max_iter=1000,
                            pca_variance_threshold=0.90, pca_n_components=None,
                            plot_file=None):
    """
    Runs TICC for a range of cluster numbers and plots the BIC curve.

    Use to choose the number of clusters for TICC on a single dataset.

    Parameters
    ----------
    data : float array
        Trajectory data [frames, frame_data]
    min_clusters : int, optional
        Default: 2.
    max_clusters : int, optional
        Default: 10.
    window_size : int, optional
        Default: 10.
    lambda_parameter : float, optional
        Default: 11e-2.
    beta : float, optional
        Default: 400.
    max_iter : int, optional
        Default: 1000.
    plot_file : str, optional

    Returns
    -------
    cluster_range : list of int
    bic_values : list of float
    """
    from pensa.clusters.ticc.TICC_solver import TICC

    cluster_range = list(range(min_clusters, max_clusters + 1))
    bic_values = []
    n_features = data.shape[1]
    do_pca = (pca_n_components is not None) or \
             (pca_variance_threshold is not None and n_features > 50)
    if do_pca:
        from sklearn.decomposition import PCA as SklearnPCA
        if pca_n_components is not None:
            pca = SklearnPCA(n_components=pca_n_components)
            data = pca.fit_transform(data)
        else:
            pca_full = SklearnPCA()
            pca_full.fit(data)
            cumvar = np.cumsum(pca_full.explained_variance_ratio_)
            n_comp = int(np.searchsorted(cumvar, pca_variance_threshold) + 1)
            pca = SklearnPCA(n_components=n_comp)
            data = pca.fit_transform(data)
        print(f"  [ticc] PCA: {data.shape[1]} components "
              f"({100*pca.explained_variance_ratio_.sum():.1f}% variance)")

    for k in cluster_range:
        model = TICC(window_size=window_size, number_of_clusters=k,
                     lambda_parameter=lambda_parameter, beta=beta,
                     maxIters=max_iter, compute_BIC=True)
        _, _, bic = model.fit_array(data)
        bic_values.append(bic)
        print(f"  k={k}  BIC={bic:.4f}")

    fig, ax = plt.subplots(1, 1, figsize=[5, 3], dpi=300)
    ax.plot(cluster_range, bic_values, 'o-')
    ax.set_xlabel('number of clusters')
    ax.set_ylabel('BIC')
    ax.set_title('TICC BIC')
    fig.tight_layout()
    if plot_file:
        fig.savefig(plot_file)

    return cluster_range, bic_values


def ticc_bic_over_combined_clusters(data_a, data_b,
                                     label_a='Sim A', label_b='Sim B',
                                     start_frame=0,
                                     min_clusters=2, max_clusters=10,
                                     window_size=10, lambda_parameter=11e-2,
                                     beta=400, max_iter=1000,
                                     pca_variance_threshold=0.90, pca_n_components=None,
                                     plot_file=None):
    """
    Runs TICC for a range of cluster numbers on combined data and plots
    the BIC curve.

    Use to choose the number of clusters for TICC on two ensembles.

    Parameters
    ----------
    data_a, data_b : float array
        Trajectory data [frames, frame_data]
    label_a, label_b : str, optional
    start_frame : int, optional
    min_clusters : int, optional
    max_clusters : int, optional
    window_size : int, optional
    lambda_parameter : float, optional
    beta : float, optional
    max_iter : int, optional
    plot_file : str, optional

    Returns
    -------
    cluster_range : list of int
    bic_values : list of float
    """
    from pensa.clusters.ticc.TICC_solver import TICC

    data = np.concatenate([data_a, data_b], 0)
    n_features = data.shape[1]
    do_pca = (pca_n_components is not None) or \
             (pca_variance_threshold is not None and n_features > 50)
    if do_pca:
        from sklearn.decomposition import PCA as SklearnPCA
        if pca_n_components is not None:
            pca = SklearnPCA(n_components=pca_n_components)
            data = pca.fit_transform(data)
        else:
            pca_full = SklearnPCA()
            pca_full.fit(data)
            cumvar = np.cumsum(pca_full.explained_variance_ratio_)
            n_comp = int(np.searchsorted(cumvar, pca_variance_threshold) + 1)
            pca = SklearnPCA(n_components=n_comp)
            data = pca.fit_transform(data)
        print(f"  [ticc] PCA: {data.shape[1]} components "
              f"({100*pca.explained_variance_ratio_.sum():.1f}% variance)")

    cluster_range = list(range(min_clusters, max_clusters + 1))
    bic_values = []
    for k in cluster_range:
        model = TICC(window_size=window_size, number_of_clusters=k,
                     lambda_parameter=lambda_parameter, beta=beta,
                     maxIters=max_iter, compute_BIC=True)
        _, _, bic = model.fit_array(data)
        bic_values.append(bic)
        print(f"  k={k}  BIC={bic:.4f}")

    fig, ax = plt.subplots(1, 1, figsize=[5, 3], dpi=300)
    ax.plot(cluster_range, bic_values, 'o-')
    ax.set_xlabel('number of clusters')
    ax.set_ylabel('BIC')
    ax.set_title('TICC BIC (combined)')
    fig.tight_layout()
    if plot_file:
        fig.savefig(plot_file)

    return cluster_range, bic_values
