import h5py
import numpy as np
import matplotlib.pyplot as plt

with h5py.File("example_1dslab.h5","r") as f:
    tally1 = f["tallies"]["user_tallies"]["test_tally_2"]
    flux_time_bins = np.array(tally1["time_edges"][:])
    flux_x_bins = np.array(tally1["x_edges"][:])
    flux_counts = np.array(tally1["counts"][:])

kendra_data = np.loadtxt("prot_flux.txt")
kendra_data = kendra_data[:,1:]

plt.figure(1)
plt.plot(flux_x_bins[1:], np.sum(flux_counts[0,:,:], axis=0), label="galvanize")
plt.plot(flux_x_bins[1:], np.sum(kendra_data, axis=0), label="jeyenne")
plt.title("Flux over x")
plt.xlabel("x [cm]")
plt.ylabel("flux [?]")
plt.yscale("log")
plt.legend()
plt.show()

plt.figure(2)
plt.plot(flux_time_bins[1:], np.sum(flux_counts[0,:,:], axis=1), label="galvanize")
plt.plot(flux_time_bins[1:], np.sum(kendra_data, axis=1), label="jeyenne")
plt.title("Flux over t")
plt.xlabel("t [shk]")
plt.ylabel("flux [?]")
plt.yscale("log")
plt.legend()
plt.show()

plt.figure(3)
plt.pcolormesh(flux_x_bins[1:], flux_time_bins[1:], flux_counts[0,:,:], cmap="RdBu")
plt.colorbar()
plt.xlabel("x [cm]")
plt.ylabel("t [shk]")
plt.title("1d slab flux")
plt.show()

plt.figure(4)
plt.pcolormesh(flux_x_bins[1:], flux_time_bins[1:], kendra_data, cmap="RdBu")
plt.colorbar()
plt.xlabel("x [cm]")
plt.ylabel("t [shk]")
plt.title("1d slab flux")
plt.show()