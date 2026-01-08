import numpy as np
import matplotlib.pyplot as plt
import h5py

with h5py.File("example_scattering_gfp2.h5","r") as f:
    tally1 = f["tallies"]["user_tallies"]["test_tally_1"]
    energy_loss = np.array(tally1["counts"][:])
    energy_bins = np.array(tally1["energy_edges"][:])
    x_bins = np.array(tally1["x_edges"][:])

plt.figure(1)
print(energy_loss.shape)
print(np.sum(energy_loss))
plt.show()