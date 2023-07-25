import wot

FULL_DS_PATH = '/home/yxy/working_space/data_repo/data/ExprMatrix.h5ad'
CELL_DAYS_PATH = '/home/yxy/working_space/data_repo/data/cell_days.txt'
CELL_SETS_PATH = '/home/yxy/working_space/data_repo/data/major_cell_sets.gmt'
TFS_PATH = '/home/yxy/working_space/data_repo/data/TFs.txt'

# Load expression dataset and subset to transcription factors
adata = wot.io.read_dataset(FULL_DS_PATH, obs=[CELL_DAYS_PATH], var_filter=TFS_PATH)
# Load transport map model and cell sets
tmap_model = wot.tmap.TransportMapModel.from_directory('/home/yxy/working_space/data_repo/tmaps/serum')
major_cell_sets = wot.io.read_sets('/home/yxy/working_space/data_repo/data/major_cell_sets.gmt', as_dict=True)

# create indicator vector for IPS cell set at day 17
target_cell_set = tmap_model.population_from_cell_sets({'IPS':major_cell_sets['IPS']}, at_time=17)
# Compute fate matrix for IPS 
fate_ds = tmap_model.fates(target_cell_set)

#  Find differentially expressed genes at day 15
results = wot.tmap.diff_exp(adata[adata.obs['day'].isin([15])], fate_ds, compare='all')
print(results[(results['t_fdr']<0.01)&(results['name1']=='IPS')].sort_values('fraction_expressed_ratio', ascending=False).head(10))