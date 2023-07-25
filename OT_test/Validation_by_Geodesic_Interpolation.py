import numpy as np
import pandas as pd

import wot

VAR_GENE_DS_PATH = '/home/yxy/working_space/data_repo/data/ExprMatrix.var.genes.h5ad'
LEARNED_GROWTH_SCORES_PATH = '/home/yxy/working_space/data_repo/tmaps/serum_g.txt'
BATCH_PATH = '/home/yxy/working_space/data_repo/data/batches.txt'
CELL_DAYS_PATH = '/home/yxy/working_space/data_repo/data/cell_days.txt'
SERUM_CELL_IDS_PATH = '/home/yxy/working_space/data_repo/data/serum_cell_ids.txt'

adata = wot.io.read_dataset(VAR_GENE_DS_PATH, obs=[CELL_DAYS_PATH, BATCH_PATH, LEARNED_GROWTH_SCORES_PATH], obs_filter=SERUM_CELL_IDS_PATH)

ot_model = wot.ot.OTModel(adata, growth_rate_field='g2',growth_iters = 1) 
interp_summary = wot.ot.compute_validation_summary(ot_model, day_triplets=[(17, 17.5, 18)])
wot.graphics.plot_ot_validation_summary_stats(interp_summary.groupby(['interval_mid', 'name'])['distance'].agg([np.mean, np.std]))

# read in and plot results
all_triplets_stats = pd.read_csv('/home/yxy/working_space/data_repo/data/serum_validation_summary_stats.txt')
wot.graphics.plot_ot_validation_summary_stats(all_triplets_stats)