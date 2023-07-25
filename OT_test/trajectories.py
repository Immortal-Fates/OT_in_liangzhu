import ipywidgets as widgets
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

import wot

# input paths
FULL_DS_PATH = '/home/yxy/working_space/data_repo/data/ExprMatrix.h5ad'
CELL_DAYS_PATH = '/home/yxy/working_space/data_repo/data/cell_days.txt'
VAR_DS_PATH = '/home/yxy/working_space/data_repo/data/ExprMatrix.var.genes.h5ad'
TMAP_PATH = '/home/yxy/working_space/data_repo/tmaps/serum'
CELL_SETS_PATH = '/home/yxy/working_space/data_repo/data/major_cell_sets.gmt'
COORDS_PATH = '/home/yxy/working_space/data_repo/data/fle_coords.txt'

tmap_model = wot.tmap.TransportMapModel.from_directory(TMAP_PATH)
cell_sets = wot.io.read_sets(CELL_SETS_PATH, as_dict=True)
#The populations object is a list of these probability vectors for different cell sets.
populations = tmap_model.population_from_cell_sets(cell_sets, at_time=12)
# The method trajectories takes a list of probability vectors and r
# eturns the trajectories containing descendant distributions 
# at later time points and ancestor distributions at earlier time points.
trajectory_ds = tmap_model.trajectories(populations)

# ### Visualize
# # Load embedding coordinates
# coord_df = pd.read_csv(COORDS_PATH, sep='\t', index_col=0)
# nbins = 500
# xrange = coord_df['x'].min(), coord_df['x'].max()
# yrange = coord_df['y'].min(), coord_df['y'].max()
# coord_df['x'] = np.floor(
#     np.interp(coord_df['x'], [xrange[0], xrange[1]], [0, nbins - 1])).astype(int)
# coord_df['y'] = np.floor(
#     np.interp(coord_df['y'], [yrange[0], yrange[1]], [0, nbins - 1])).astype(int)
# trajectory_ds.obs = trajectory_ds.obs.join(coord_df)
# # Visualize trajectories
# trajectory_dropdown = widgets.Dropdown(
#     options=trajectory_ds.var.index,
#     description='Trajectory:'
# )

# def update_trajectory_vis(name):
#     figure = plt.figure(figsize=(10, 10))
#     plt.axis('off')
#     plt.tight_layout()
#     plt.title(name)
#     plt.scatter(coord_df['x'], coord_df['y'], c='#f0f0f0',
#                    s=4, marker=',', edgecolors='none', alpha=0.8)
#     binned_df = trajectory_ds.obs.copy()
#     binned_df['values'] = trajectory_ds[:, name].X
#     binned_df = binned_df.groupby(['x', 'y'], as_index=False).sum()
#     plt.scatter(binned_df['x'], binned_df['y'], c=binned_df['values'],
#                    s=6, marker=',', edgecolors='none', vmax=binned_df['values'].quantile(0.975))
#     plt.colorbar().ax.set_title('Trajectory')

# widgets.interact(update_trajectory_vis, name=trajectory_dropdown)

# ### We can also generate a trajectory movie across time
# from matplotlib import animation
# from IPython.display import Video
# name = 'IPS'
# movie_file_name = '{}_trajectory.gif'.format(name)
# figure = plt.figure(figsize=(10, 10))
# plt.axis('off')

# plt.tight_layout()
# unique_days = trajectory_ds.obs['day'].unique()

# binned_df = trajectory_ds.obs.copy()
# binned_df['values'] = trajectory_ds[:, name].X
# binned_df = binned_df.groupby(['x', 'y'], as_index=False).sum()
# vmax=binned_df['values'].quantile(0.975)

# def animate(i):
#     _trajectory_ds = trajectory_ds[trajectory_ds.obs['day']==unique_days[i]]
#     binned_df = _trajectory_ds.obs.copy()
#     plt.suptitle('Day {}, {:,} cells'.format(unique_days[i], _trajectory_ds.shape[0]), y=0.95)
#     plt.scatter(coord_df['x'], coord_df['y'], c='#f0f0f0', s=4, marker=',', edgecolors='none')
#     binned_df['values'] = _trajectory_ds[:, name].X
#     binned_df = binned_df.groupby(['x', 'y'], as_index=False).sum()
#     plt.scatter(binned_df['x'], binned_df['y'], c=binned_df['values'],
#                 s=6, marker=',', edgecolors='none', vmin=0, vmax=vmax)
 
    
# anim = animation.FuncAnimation(figure, func=animate, frames=range(0, len(unique_days)), init_func=lambda **args:None, repeat=False, interval=400)
# anim.save(movie_file_name,writer='pillow')
# plt.close(figure)
# Video(movie_file_name)

### Expression trends along trajectories
# We now show how to compute trends in expression along trajectories
#Load expression data
adata = wot.io.read_dataset(FULL_DS_PATH) 

#Compute trends for all genes
trajectory_trends = wot.tmap.trajectory_trends_from_trajectory(trajectory_ds, adata)

# Save each trajectory in a separate file
for i in range(len(trajectory_trends)):
    wot.io.write_dataset(trajectory_trends[i], trajectory_ds.var.index[i] + '_trends.txt')

# Read in trajectory trends

trajectory_trend_datasets = []
trajectory_names = []

for i in range(trajectory_ds.shape[1]):
    trajectory_names.append(trajectory_ds.var.index[i]) 
    trajectory_trend_datasets.append(wot.io.read_dataset(trajectory_ds.var.index[i] + '_trends.txt'))

# Shared Ancestry
adata_var = wot.io.read_dataset(VAR_DS_PATH) 
divergence_df = wot.tmap.trajectory_divergence(adata_var, trajectory_ds, distance_metric='total_variation')
# Plot divergence

divergence_df['name'] = divergence_df['name1'].str.split('/').str.get(0) + ' vs. ' + divergence_df['name2'].str.split('/').str.get(
        0)
plt.figure(figsize=(10, 10))
plt.xlabel("Day")
plt.ylabel("Distance")
for p, d in divergence_df.groupby('name'):
    plt.plot(d['day2'], d['distance'], '-o', label=p)
plt.legend(loc='best')
plt.show()