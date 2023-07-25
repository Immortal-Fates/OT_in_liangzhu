import ipywidgets as widgets
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

import wot


# Path to input files
FLE_COORDS_PATH ="/home/yxy/working_space/data_repo/data/fle_coords.txt"
FULL_DS_PATH = '/home/yxy/working_space/data_repo/data/ExprMatrix.h5ad'
VAR_DS_PATH = '/home/yxy/working_space/data_repo/data/ExprMatrix.var.genes.h5ad'
CELL_DAYS_PATH = '/home/yxy/working_space/data_repo/data/cell_days.txt'
GENE_SETS_PATH = '/home/yxy/working_space/data_repo/data/gene_sets.gmx'
GENE_SET_SCORES_PATH = '/home/yxy/working_space/data_repo/data/gene_set_scores.csv'
CELL_SETS_PATH = '/home/yxy/working_space/data_repo/data/cell_sets.gmt'

coord_df = pd.read_csv(FLE_COORDS_PATH, index_col='id', sep='\t')
days_df = pd.read_csv(CELL_DAYS_PATH, index_col='id', sep='\t')

# Read expression matrix, cell days, and 2-d coordinates
adata = wot.io.read_dataset(FULL_DS_PATH, obs=[days_df,coord_df])
unique_days = adata.obs['day'].unique()
unique_days = unique_days[np.isnan(unique_days) == False]

# # plot visualization coordinates
# figure = plt.figure(figsize=(10, 10))
# plt.axis('off')
# plt.tight_layout()
# plt.scatter(adata.obs['x'], adata.obs['y'],c=adata.obs['day'],
#                s=4, marker=',', edgecolors='none', alpha=0.8)
# cb = plt.colorbar()
# cb.ax.set_title('Day')
# plt.show()

# #generate a movie that shows cells bu day
# from matplotlib import animation

# from IPython.display import Video
# movie_file_name = '/home/yxy/working_space/OT_test/cells_by_day.gif'
# coord_days_df = coord_df.join(days_df)
# figure = plt.figure(figsize=(10, 10))
# plt.axis('off')
# plt.tight_layout()

# def animate(i):
#     day_df = coord_days_df[coord_days_df['day']==unique_days[i]]
#     plt.suptitle('Day {}, {:,} cells'.format(unique_days[i], day_df.shape[0]), y=0.95)
#     plt.plot(coord_df['x'], coord_df['y'], ',', color='#f0f0f0', ms=4)
#     plt.plot(day_df['x'], day_df['y'], ',', color='black', ms=4)   
    
    
    
# anim = animation.FuncAnimation(figure, func=animate, frames=range(0, len(unique_days)), init_func=lambda **args:None, repeat=False, interval=400)
# anim.save(movie_file_name,writer='ffmpeg',fps=1000/50)
# plt.close(figure)
# Video(movie_file_name)

###compute gene signature scores
# gs = wot.io.read_sets(GENE_SETS_PATH, adata.var.index.values)
# gene_set_scores_df = pd.DataFrame(index=adata.obs.index)
# for j in range(gs.shape[1]):
#     gene_set_name = str(gs.var.index.values[j])
#     result = wot.score_gene_sets(ds=adata, gs=gs[:, [j]], permutations=0, method='mean_z_score')
#     gene_set_scores_df[gene_set_name] = result['score']
# gene_set_scores_df.to_csv(GENE_SET_SCORES_PATH, index_label='id')

# gene_set_scores_df = pd.read_csv(GENE_SET_SCORES_PATH,index_col='id')
# gene_set_dropdown = widgets.Dropdown(
#     options=gene_set_scores_df.columns,
#     description='Gene Set:'
# )

# gene_set_scores_df = gene_set_scores_df.join(coord_df).join(days_df)
# day_selector = widgets.SelectionRangeSlider(
#     options=unique_days,
#     continuous_update=False,
#     index=(0,len(unique_days)-1),
#     description='Days'
# )

# def update_gene_set_vis(name, days):
#     gene_set_score_coords = gene_set_scores_df[(gene_set_scores_df['day']>=days[0]) & (gene_set_scores_df['day']<=days[1])]
#     figure = plt.figure(figsize=(10, 10))
#     plt.axis('off')
#     plt.tight_layout()
#     plt.title(name + ', days {}-{}'.format(days[0], days[1]))
#     plt.scatter(coord_df['x'], coord_df['y'], c='#f0f0f0',
#                    s=4, marker=',', edgecolors='none', alpha=0.8)
#     plt.scatter(gene_set_score_coords['x'], gene_set_score_coords['y'], c=gene_set_score_coords[name],
#                    s=4, marker=',', edgecolors='none')
#     cb = plt.colorbar()
#     cb.ax.set_title('Signature')
#     figure2 = plt.figure(figsize=(10, 5))
#     plt.title(name + ', days {}-{}'.format(days[0], days[1]))
#     plt.hist(gene_set_score_coords[name])
#     # plt.show()
#     return figure, figure2

# widgets.interact(update_gene_set_vis, name=gene_set_dropdown, days=day_selector)

#没看懂这几个interactive
### Cell Sets
# # Load cell sets
# cell_sets = wot.io.read_sets(CELL_SETS_PATH)

# # Visualize Cell Sets 
# cell_set_dropdown = widgets.Dropdown(
#     options=cell_sets.var.index,
#     description='Cell Set:'
# )

# day_selector = widgets.SelectionRangeSlider(
#     options=unique_days,
#     continuous_update=False,
#     index=(0,len(unique_days)-1),
#     description='Days'
# )

# def update_cell_set_vis(name, days):
#     cell_set = cell_sets[:, name]
#     cell_set_coords = cell_set[cell_set.X>0].obs.join(coord_df).join(days_df)
#     cell_set_coords = cell_set_coords[(cell_set_coords['day']>=days[0]) & (cell_set_coords['day']<=days[1])]
#     figure = plt.figure(figsize=(10, 10))
#     plt.axis('off')
#     plt.tight_layout()
#     plt.title(name + ', days {}-{}, {:,} cells'.format(days[0], days[1], cell_set_coords.shape[0]))
#     plt.scatter(coord_df['x'], coord_df['y'], c='#f0f0f0',
#                    s=4, marker=',', edgecolors='none', alpha=0.8)
#     plt.scatter(cell_set_coords['x'], cell_set_coords['y'], c=cell_set_coords['day'],
#                    s=4, marker=',', edgecolors='none', vmin=unique_days[0],  vmax=unique_days[len(unique_days)-1])
#     cb = plt.colorbar()
#     cb.ax.set_title('Day')

# widgets.interact(update_cell_set_vis, name=cell_set_dropdown, days=day_selector)


# adata_var = wot.io.read_dataset(VAR_DS_PATH, obs=[days_df])
# import pegasus as pg
# pg.pca(adata_var, features=None)
# pg.neighbors(adata_var)
# pg.diffmap(adata_var)
# pg.fle(adata_var)

# figure = plt.figure(figsize=(10, 10))
# plt.axis('off')
# plt.tight_layout()
# coords = adata_var.obsm['X_fle']
# adata_var.obs = adata_var.obs.join(adata.obs)
# plt.scatter(coords[:, 0], coords[:, 1], c=adata_var.obs['day'], s=4, marker=',', edgecolors='none', alpha=0.8)
# cb = plt.colorbar()
# cb.ax.set_title('Day')