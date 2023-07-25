import wot

tmap_model = wot.tmap.TransportMapModel.from_directory('/home/yxy/working_space/data_repo/tmaps/serum')

# We can now easily compute the coupling between any pair of time-points as follows:

gamma_8_10 = tmap_model.get_coupling(8, 10)
print(gamma_8_10)
print(gamma_8_10.obs)