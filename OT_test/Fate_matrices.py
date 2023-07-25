import math

import ipywidgets as widgets
import numpy as np
from matplotlib import pyplot as plt

import wot


# Load transport map model and cell sets
tmap_model = wot.tmap.TransportMapModel.from_directory('/home/yxy/working_space/data_repo/tmaps/serum')
cell_sets = wot.io.read_sets('/home/yxy/working_space/data_repo/data/major_cell_sets.gmt', as_dict=True)
# create indicator vectors for each cell set
target_destinations = tmap_model.population_from_cell_sets(cell_sets, at_time=18)
# Finally, we compute the fate matrices for each earlier time point ti<tj
#  all at once with the fates function.

fate_ds = tmap_model.fates(target_destinations)

