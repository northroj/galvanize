import numpy as np
import matplotlib.pyplot as plt
import h5py


#"""
with h5py.File("example_scattering_gfp2.h5","r") as f:

    tally2 = f["tallies"]["user_tallies"]["escape_energy"]
    energy_bins_gfp2 = np.array(tally2["energy_edges"][:])
    escape_energy_gfp2 = np.array(tally2["counts"][:])
#"""


#"""
with h5py.File("example_scattering_csd.h5","r") as f:

    tally2 = f["tallies"]["user_tallies"]["escape_energy"]
    energy_bins_csd = np.array(tally2["energy_edges"][:])
    escape_energy_csd = np.array(tally2["counts"][:])
#"""


plt.figure(2)
plt.plot(energy_bins_csd[1:], escape_energy_csd[0,:], color="red", label="csd")
plt.plot(energy_bins_gfp2[1:], escape_energy_gfp2[0,:], color="green", label="gfp2")
plt.xlim([3090,3105])
#plt.ylim([0,0.2])
plt.xlabel("Energy [keV]")
plt.legend()
plt.yscale("log")
plt.title("3.5 MeV protons slowing down in 0.001 cm of tungsten")
plt.show()