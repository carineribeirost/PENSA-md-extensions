import argparse
import numpy as np

from pensa.features import \
    read_structure_features
from pensa.clusters import \
    obtain_combined_clusters, \
    write_cluster_traj, \
    wss_over_number_of_combined_clusters, \
    pdc_sensitivity_for_combined_clusters, \
    ticc_bic_over_combined_clusters


# ----------------------- #
# --- Workflow helpers --- #
# ----------------------- #

def workflow_clusters(args, data_a, data_b, ftype):
    """
    Run parameter selection and main clustering for one feature type.

    Parameters
    ----------
    args : argparse.Namespace
    data_a, data_b : dict
        Feature dictionaries from read_structure_features.
    ftype : str
        Feature type key (e.g. 'bb-torsions', 'sc-torsions', 'bb-distances').

    Returns
    -------
    cc : CombinedClusterResult
    """
    algo = args.algorithm
    base_plot = f"{args.out_plots}_{algo}"

    # WSS elbow for k-means / rspace
    if args.wss and algo in ('kmeans', 'rspace'):
        wss_over_number_of_combined_clusters(
            data_a[ftype], data_b[ftype],
            label_a=args.label_a, label_b=args.label_b,
            start_frame=args.start_frame,
            algorithm=algo,
            max_iter=100, num_repeats=5,
            max_num_clusters=args.max_num_clusters,
            plot_file=f"{base_plot}_wss_{ftype}.pdf"
        )

    # CLoNe pdc sensitivity sweep
    if args.pdc_sensitivity and algo == 'clone':
        pdc_sensitivity_for_combined_clusters(
            data_a[ftype], data_b[ftype],
            label_a=args.label_a, label_b=args.label_b,
            start_frame=args.start_frame,
            n_resize=args.n_resize, filt=args.filt,
            plot_file=f"{base_plot}_pdc_sensitivity_{ftype}.pdf"
        )

    # TICC BIC sweep
    if args.ticc_bic and algo == 'ticc':
        ticc_bic_over_combined_clusters(
            data_a[ftype], data_b[ftype],
            label_a=args.label_a, label_b=args.label_b,
            start_frame=args.start_frame,
            min_clusters=args.ticc_min_clusters,
            max_clusters=args.ticc_max_clusters,
            window_size=args.window_size,
            lambda_parameter=args.lambda_parameter,
            beta=args.beta, max_iter=args.ticc_max_iter,
            pca_variance_threshold=args.pca_variance_threshold,
            pca_n_components=args.pca_n_components,
            plot_file=f"{base_plot}_bic_{ftype}.pdf"
        )

    # Main combined clustering
    if algo == 'kmeans':
        cc = obtain_combined_clusters(
            data_a[ftype], data_b[ftype],
            args.label_a, args.label_b,
            args.start_frame,
            algorithm='kmeans', max_iter=100,
            num_clusters=args.write_num_clusters,
            saveas=f"{base_plot}_combined-clusters_{ftype}.pdf"
        )
    elif algo == 'rspace':
        cc = obtain_combined_clusters(
            data_a[ftype], data_b[ftype],
            args.label_a, args.label_b,
            args.start_frame,
            algorithm='rspace',
            saveas=f"{base_plot}_combined-clusters_{ftype}.pdf"
        )
    elif algo == 'clone':
        cc = obtain_combined_clusters(
            data_a[ftype], data_b[ftype],
            args.label_a, args.label_b,
            args.start_frame,
            algorithm='clone',
            pdc=args.pdc, n_resize=args.n_resize, filt=args.filt,
            saveas=f"{base_plot}_combined-clusters_{ftype}.pdf"
        )
    elif algo == 'ticc':
        cc = obtain_combined_clusters(
            data_a[ftype], data_b[ftype],
            args.label_a, args.label_b,
            args.start_frame,
            algorithm='ticc',
            num_clusters=args.ticc_clusters,
            window_size=args.window_size,
            lambda_parameter=args.lambda_parameter,
            beta=args.beta, max_iter=args.ticc_max_iter,
            pca_variance_threshold=args.pca_variance_threshold,
            pca_n_components=args.pca_n_components,
            saveas=f"{base_plot}_combined-clusters_{ftype}.pdf"
        )

    return cc


# -------------#
# --- MAIN --- #
# -------------#

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    # --- Trajectory files ---
    parser.add_argument("--ref_file_a", type=str,
                        default='traj/rhodopsin_arrbound_receptor.gro')
    parser.add_argument("--trj_file_a", type=str,
                        default='traj/rhodopsin_arrbound_receptor.xtc')
    parser.add_argument("--ref_file_b", type=str,
                        default='traj/rhodopsin_gibound_receptor.gro')
    parser.add_argument("--trj_file_b", type=str,
                        default='traj/rhodopsin_gibound_receptor.xtc')
    parser.add_argument("--label_a", type=str, default='Sim A')
    parser.add_argument("--label_b", type=str, default='Sim B')

    # --- Output paths ---
    parser.add_argument("--out_plots", type=str,
                        default='plots/rhodopsin_receptor')
    parser.add_argument("--out_results", type=str,
                        default='results/rhodopsin_receptor')
    parser.add_argument("--out_frames_a", type=str,
                        default='clusters/rhodopsin_arrbound_receptor')
    parser.add_argument("--out_frames_b", type=str,
                        default='clusters/rhodopsin_gibound_receptor')

    # --- Shared options ---
    parser.add_argument("--start_frame", type=int, default=0)
    parser.add_argument("--algorithm", type=str, default='kmeans',
                        choices=['kmeans', 'rspace', 'clone', 'ticc'],
                        help='Clustering algorithm.')

    # --- k-means / rspace ---
    parser.add_argument("--max_num_clusters", type=int, default=12,
                        help='Max clusters for WSS sweep (kmeans/rspace).')
    parser.add_argument("--write_num_clusters", type=int, default=2,
                        help='Number of clusters for the main kmeans/rspace run.')

    # --- CLoNe ---
    parser.add_argument("--pdc", type=float, default=4.0,
                        help='Percentile distance cutoff for CLoNe.')
    parser.add_argument("--n_resize", type=int, default=4,
                        help='Resize parameter for CLoNe.')
    parser.add_argument("--filt", type=float, default=0.1,
                        help='Outlier threshold for CLoNe.')

    # --- TICC ---
    parser.add_argument("--ticc_clusters", type=int, default=5,
                        help='Number of clusters for TICC.')
    parser.add_argument("--window_size", type=int, default=10,
                        help='Sliding window size for TICC.')
    parser.add_argument("--lambda_parameter", type=float, default=11e-2,
                        help='Sparsity parameter for TICC.')
    parser.add_argument("--beta", type=float, default=400,
                        help='Temporal consistency parameter for TICC.')
    parser.add_argument("--ticc_max_iter", type=int, default=1000,
                        help='Max iterations for TICC.')
    parser.add_argument("--ticc_min_clusters", type=int, default=2,
                        help='Min clusters for TICC BIC sweep.')
    parser.add_argument("--ticc_max_clusters", type=int, default=10,
                        help='Max clusters for TICC BIC sweep.')
    parser.add_argument("--pca_variance_threshold", type=float, default=0.90,
                        help='Fraction of variance for automatic PCA pre-reduction '
                             'before TICC (applied when n_features > 50). Default: 0.90.')
    parser.add_argument("--pca_n_components", type=int, default=None,
                        help='Override: use exactly this many PCA components before TICC '
                             '(takes priority over --pca_variance_threshold).')

    # --- Flags ---
    parser.add_argument('--write', dest='write', action='store_true',
                        help='Write cluster trajectory files.')
    parser.add_argument('--no-write', dest='write', action='store_false')
    parser.add_argument('--wss', dest='wss', action='store_true',
                        help='Run WSS elbow sweep (kmeans/rspace).')
    parser.add_argument('--no-wss', dest='wss', action='store_false')
    parser.add_argument('--pdc-sensitivity', dest='pdc_sensitivity',
                        action='store_true',
                        help='Run CLoNe pdc sensitivity sweep.')
    parser.add_argument('--no-pdc-sensitivity', dest='pdc_sensitivity',
                        action='store_false')
    parser.add_argument('--ticc-bic', dest='ticc_bic', action='store_true',
                        help='Run TICC BIC sweep.')
    parser.add_argument('--no-ticc-bic', dest='ticc_bic', action='store_false')

    parser.set_defaults(write=True, wss=True, pdc_sensitivity=False, ticc_bic=False)
    args = parser.parse_args()

    # -- FEATURES --

    feat_a, data_a = read_structure_features(
        args.ref_file_a, args.trj_file_a,
        args.start_frame, cossin=True
    )
    feat_b, data_b = read_structure_features(
        args.ref_file_b, args.trj_file_b,
        args.start_frame, cossin=True
    )
    print('Feature dimensions from', args.trj_file_a)
    for k in data_a.keys():
        print(k, data_a[k].shape)
    print('Feature dimensions from', args.trj_file_b)
    for k in data_b.keys():
        print(k, data_b[k].shape)

    # -- BACKBONE TORSIONS --

    print('\nBACKBONE TORSIONS')
    cc = workflow_clusters(args, data_a, data_b, 'bb-torsions')
    cidx, cond, oidx = cc.cidx, cc.cond, cc.oidx
    np.savetxt(
        f"{args.out_results}_{args.algorithm}_combined-cluster-indices_bb-torsions.csv",
        np.array([cidx, cond, oidx], dtype=int).T,
        delimiter=', ', fmt='%i',
        header='Cluster, Condition, Index within condition'
    )
    if args.write:
        write_cluster_traj(
            cidx[cond == 0], args.ref_file_a, args.trj_file_a,
            f"{args.out_frames_a}_{args.algorithm}_bb-torsions", args.start_frame
        )
        write_cluster_traj(
            cidx[cond == 1], args.ref_file_b, args.trj_file_b,
            f"{args.out_frames_b}_{args.algorithm}_bb-torsions", args.start_frame
        )

    # -- SIDECHAIN TORSIONS --

    print('\nSIDECHAIN TORSIONS')
    cc = workflow_clusters(args, data_a, data_b, 'sc-torsions')
    cidx, cond, oidx = cc.cidx, cc.cond, cc.oidx
    np.savetxt(
        f"{args.out_results}_{args.algorithm}_combined-cluster-indices_sc-torsions.csv",
        np.array([cidx, cond, oidx], dtype=int).T,
        delimiter=', ', fmt='%i',
        header='Cluster, Condition, Index within condition'
    )
    if args.write:
        write_cluster_traj(
            cidx[cond == 0], args.ref_file_a, args.trj_file_a,
            f"{args.out_frames_a}_{args.algorithm}_sc-torsions", args.start_frame
        )
        write_cluster_traj(
            cidx[cond == 1], args.ref_file_b, args.trj_file_b,
            f"{args.out_frames_b}_{args.algorithm}_sc-torsions", args.start_frame
        )

    # -- BACKBONE C-ALPHA DISTANCES --
    # NOTE: commented out — 41k features impractical for clustering.
    # For distance analysis use compare_feature_distributions.py instead.

    # print('\nBACKBONE C-ALPHA DISTANCES')
    # cc = workflow_clusters(args, data_a, data_b, 'bb-distances')
    # cidx, cond, oidx = cc.cidx, cc.cond, cc.oidx
    # np.savetxt(
    #     f"{args.out_results}_{args.algorithm}_combined-cluster-indices_bb-distances.csv",
    #     np.array([cidx, cond, oidx], dtype=int).T,
    #     delimiter=', ', fmt='%i',
    #     header='Cluster, Condition, Index within condition'
    # )
    # if args.write:
    #     write_cluster_traj(
    #         cidx[cond == 0], args.ref_file_a, args.trj_file_a,
    #         f"{args.out_frames_a}_{args.algorithm}_bb-distances", args.start_frame
    #     )
    #     write_cluster_traj(
    #         cidx[cond == 1], args.ref_file_b, args.trj_file_b,
    #         f"{args.out_frames_b}_{args.algorithm}_bb-distances", args.start_frame
    #     )
