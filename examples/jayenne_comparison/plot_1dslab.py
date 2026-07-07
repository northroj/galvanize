import h5py
import numpy as np
import matplotlib.pyplot as plt

with h5py.File("example_1dslab.h5","r") as f:
    tally1 = f["tallies"]["user_tallies"]["test_tally_2"]
    flux_time_bins = np.array(tally1["time_edges"][:])
    flux_x_bins = np.array(tally1["x_edges"][:])
    flux_counts = np.array(tally1["counts"][:])

flux_counts = flux_counts * 2.1e7

kendra_data = np.loadtxt("1d_flux_cartesian.txt")
kendra_data = kendra_data[:,1:]

plt.figure(1)
plt.plot(flux_x_bins[1:], np.sum(flux_counts[0,:,:], axis=0), label="galvanize")
plt.plot(flux_x_bins[1:], np.sum(kendra_data, axis=0), label="jayenne")
plt.title("Flux over x")
plt.xlabel("x [cm]")
plt.ylabel("flux [?]")
plt.yscale("log")
plt.legend()
plt.show()

plt.figure(2)
plt.plot(flux_time_bins[1:], np.sum(flux_counts[0,:,:], axis=1), label="galvanize")
plt.plot(flux_time_bins[1:], np.sum(kendra_data, axis=1), label="jayenne")
plt.title("Flux over t")
plt.xlabel("t [shk]")
plt.ylabel("flux [?]")
plt.yscale("log")
plt.legend()
plt.show()

plt.figure(3)
plt.pcolormesh(flux_x_bins[1:], flux_time_bins[1:], np.log(flux_counts[0,:,:]+1e-10), cmap="RdBu")
plt.colorbar()
plt.xlabel("x [cm]")
plt.ylabel("t [shk]")
plt.title("1d slab flux - galvanize")
plt.show()

plt.figure(4)
plt.pcolormesh(flux_x_bins[1:], flux_time_bins[1:], np.log(kendra_data+1e-10), cmap="RdBu")
plt.colorbar()
plt.xlabel("x [cm]")
plt.ylabel("t [shk]")
plt.title("1d slab flux - Jayenne")
plt.show()