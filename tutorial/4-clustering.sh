#!/bin/bash

PENSA_PYTHON="/home/carine/PENSA-CLoNe/fork/.venv/bin/python"
PENSA_PATH="/home/carine/refactor-pensa/pensa"
SCRIPT="$PENSA_PATH/scripts/calculate_combined_clusters.py"

mkdir -p plots
mkdir -p clusters
mkdir -p results

# --- k-means clustering (WSS elbow to guide k selection) ---
PYTHONPATH="$PENSA_PATH" $PENSA_PYTHON $SCRIPT --write --wss \
	--ref_file_a 'traj/condition-a_receptor.gro' \
	--trj_file_a 'traj/condition-a_receptor.xtc' \
	--ref_file_b 'traj/condition-b_receptor.gro' \
	--trj_file_b 'traj/condition-b_receptor.xtc' \
	--label_a 'A' \
	--label_b 'B' \
	--out_plots 'plots/receptor' \
	--out_results 'results/receptor' \
	--out_frames_a 'clusters/receptor' \
	--out_frames_b 'clusters/receptor' \
	--start_frame 2000 \
	--feature_type 'bb-torsions' \
	--algorithm 'kmeans' \
	--max_num_clusters 12 \
	--write_num_clusters 2

# --- CLoNe clustering ---
# Step 1: pdc sensitivity sweep — shows how the number of detected clusters
#         varies with pdc. Use this plot to choose pdc before the main run.
PYTHONPATH="$PENSA_PATH" $PENSA_PYTHON $SCRIPT --no-write --pdc-sensitivity \
	--ref_file_a 'traj/condition-a_receptor.gro' \
	--trj_file_a 'traj/condition-a_receptor.xtc' \
	--ref_file_b 'traj/condition-b_receptor.gro' \
	--trj_file_b 'traj/condition-b_receptor.xtc' \
	--label_a 'A' \
	--label_b 'B' \
	--out_plots 'plots/receptor' \
	--start_frame 2000 \
	--feature_type 'bb-torsions' \
	--algorithm 'clone'

# Step 2: main CLoNe run with chosen pdc.
# Adjust --pdc based on the sensitivity plot above.
PYTHONPATH="$PENSA_PATH" $PENSA_PYTHON $SCRIPT --write --no-wss \
	--ref_file_a 'traj/condition-a_receptor.gro' \
	--trj_file_a 'traj/condition-a_receptor.xtc' \
	--ref_file_b 'traj/condition-b_receptor.gro' \
	--trj_file_b 'traj/condition-b_receptor.xtc' \
	--label_a 'A' \
	--label_b 'B' \
	--out_plots 'plots/receptor' \
	--out_results 'results/receptor' \
	--out_frames_a 'clusters/receptor' \
	--out_frames_b 'clusters/receptor' \
	--start_frame 2000 \
	--feature_type 'bb-torsions' \
	--algorithm 'clone' \
	--pdc 4 \
	--n_resize 4 \
	--filt 0.1

# --- TICC clustering ---
# Step 1: BIC sweep to choose number of clusters.
# PCA pre-reduction is applied automatically (90% variance rule) when
# n_features > 50. Override with --pca_variance_threshold or --pca_n_components.
PYTHONPATH="$PENSA_PATH" $PENSA_PYTHON $SCRIPT --no-write --ticc-bic \
	--ref_file_a 'traj/condition-a_receptor.gro' \
	--trj_file_a 'traj/condition-a_receptor.xtc' \
	--ref_file_b 'traj/condition-b_receptor.gro' \
	--trj_file_b 'traj/condition-b_receptor.xtc' \
	--label_a 'A' \
	--label_b 'B' \
	--out_plots 'plots/receptor' \
	--start_frame 2000 \
	--feature_type 'bb-torsions' \
	--algorithm 'ticc' \
	--ticc_min_clusters 2 \
	--ticc_max_clusters 8

# Step 2: main TICC run with chosen number of clusters.
# Adjust --ticc_clusters based on the BIC plot above.
PYTHONPATH="$PENSA_PATH" $PENSA_PYTHON $SCRIPT --write --no-wss \
	--ref_file_a 'traj/condition-a_receptor.gro' \
	--trj_file_a 'traj/condition-a_receptor.xtc' \
	--ref_file_b 'traj/condition-b_receptor.gro' \
	--trj_file_b 'traj/condition-b_receptor.xtc' \
	--label_a 'A' \
	--label_b 'B' \
	--out_plots 'plots/receptor' \
	--out_results 'results/receptor' \
	--out_frames_a 'clusters/receptor' \
	--out_frames_b 'clusters/receptor' \
	--start_frame 2000 \
	--feature_type 'bb-torsions' \
	--algorithm 'ticc' \
	--ticc_clusters 5 \
	--window_size 10 \
	--lambda_parameter 0.11 \
	--beta 400
