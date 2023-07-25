import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import wot

###Step 1: Construct initial estimate of cell growth rates (optional)
# load proliferation and apoptosis scores
gene_set_scores = pd.read_csv('/home/yxy/working_space/data_repo/data/gene_set_scores.csv', index_col=0)
proliferation=gene_set_scores['Cell.cycle']
apoptosis = gene_set_scores['Apoptosis']

# apply logistic function to transform to birth rate and death rate
def logistic(x, L, k, x0=0):
    f = L / (1 + np.exp(-k * (x - x0)))
    return f
def gen_logistic(p, beta_max, beta_min, pmax, pmin, center, width):
    return beta_min + logistic(p, L=beta_max - beta_min, k=4 / width, x0=center)

def beta(p, beta_max=1.7, beta_min=0.3, pmax=1.0, pmin=-0.5, center=0.25):
    return gen_logistic(p, beta_max, beta_min, pmax, pmin, center, width=0.5)

def delta(a, delta_max=1.7, delta_min=0.3, amax=0.5, amin=-0.4, center=0.1):
    return gen_logistic(a, delta_max, delta_min, amax, amin, center,
                          width=0.2)

birth = beta(proliferation)
death = delta(apoptosis)


# growth rate is given by 
gr = np.exp(birth-death)
growth_rates_df = pd.DataFrame(index=gene_set_scores.index, data={'cell_growth_rate':gr})
growth_rates_df.to_csv('/home/yxy/working_space/data_repo/data/My_growth_gs_init.txt')


### Step 2: Compute transport maps
VAR_GENE_DS_PATH = '/home/yxy/working_space/data_repo/data/ExprMatrix.var.genes.h5ad'
CELL_DAYS_PATH = '/home/yxy/working_space/data_repo/data/cell_days.txt'
SERUM_CELL_IDS_PATH = '/home/yxy/working_space/data_repo/data/serum_cell_ids.txt'
CELL_GROWTH_PATH = '/home/yxy/working_space/data_repo/data/growth_gs_init.txt'

# load data
adata = wot.io.read_dataset(VAR_GENE_DS_PATH, obs=[CELL_DAYS_PATH, CELL_GROWTH_PATH], obs_filter=SERUM_CELL_IDS_PATH)
adata.shape

# create OTModel
ot_model = wot.ot.OTModel(adata,epsilon = 0.05, lambda1 = 1,lambda2 = 50) 
# Compute a single transport map from day 7 to 7.5
tmap_annotated = ot_model.compute_transport_map(7,7.5)
# row annotations include cell growth rates
# print(tmap_annotated.obs)

# columns annotated by cell barcodes
# tmap_annotated.var

# Visualize how growth rates change with growth iterations
plt.scatter(tmap_annotated.obs['g0'],tmap_annotated.obs['g1'])
plt.xlabel("g0")
plt.ylabel("g1")
plt.title("Input vs Output Growth Rates")
plt.show()

ot_model_gr2 = wot.ot.OTModel(adata,epsilon = 0.05, lambda1 = 1,lambda2 = 50,growth_iters=2) 
tmap_anno_gr2 = ot_model_gr2.compute_transport_map(7,7.5)
colsums = tmap_anno_gr2.X.sum(axis=0);
plt.hist(colsums)
plt.show()