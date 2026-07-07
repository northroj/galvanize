import numpy as np
import matplotlib.pyplot as plt
import h5py


test_case = 4 # 1 = 20 cm, 2 = 5, 3 = 0.5, 4 = 0.1

if test_case == 1:
    csd_file = "thin_slab_20_csd.h5"
    gfp2_file = "thin_slab_20_gfp2.h5"
    gfp4_file = "thin_slab_20_gfp4.h5"
    analog_file = "thin_slab_20_analog.h5"
    x1, x2 = 1430000, 1510000
    slab = 20
if test_case == 2:
    csd_file = "thin_slab_5_csd.h5"
    gfp2_file = "thin_slab_5_gfp2.h5"
    gfp4_file = "thin_slab_5_gfp4.h5"
    analog_file = "thin_slab_5_analog.h5"
    x1, x2 = 1620000, 1665000
    slab = 5
if test_case == 3:
    csd_file = "thin_slab_05_csd.h5"
    gfp2_file = "thin_slab_05_gfp2.h5"
    gfp4_file = "thin_slab_05_gfp4.h5"
    analog_file = "thin_slab_05_analog.h5"
    x1, x2 = 1680000, 1700000
    slab = 0.5
if test_case == 4:
    csd_file = "thin_slab_01_csd.h5"
    gfp2_file = "thin_slab_01_gfp2.h5"
    gfp4_file = "thin_slab_01_gfp4.h5"
    analog_file = "thin_slab_01_analog.h5"
    x1, x2 = 1694000, 1700000
    slab = 0.1

def extract_tally(file_name):
    with h5py.File(file_name,"r") as f:
        tally = f["tallies"]["user_tallies"]["escape_energy"]
        energy_bins = np.array(tally["energy_edges"][:])
        escape_energy = np.array(tally["counts"][:])
    return energy_bins, escape_energy

csd_bins, csd_counts = extract_tally(csd_file)
gfp2_bins, gfp2_counts = extract_tally(gfp2_file)
gfp4_bins, gfp4_counts = extract_tally(gfp4_file)
analog_bins, analog_counts = extract_tally(analog_file)

test_bins, test_counts = extract_tally("test_newgfp.h5")

plt.figure(1)
plt.plot(csd_bins[1:], csd_counts[0,:], color="red", label="csd")
plt.plot(gfp2_bins[1:], gfp2_counts[0,:], color="blue", label="gfp2")
plt.plot(gfp4_bins[1:], gfp4_counts[0,:], color="green", label="gfp4")
plt.plot(analog_bins[1:], analog_counts[0,:], color="black", label="analog")
#plt.plot(test_bins[1:], test_counts[0,:], color="black", label="test")
plt.xlim([x1, x2])
#plt.ylim(bottom=1e-4, top=None)
plt.xlabel("Energy [keV]")
plt.title(f"{slab} cm tungsten w/ 1.7e7 keV protons")
plt.yscale("log")
plt.grid()
plt.legend()
plt.show()
