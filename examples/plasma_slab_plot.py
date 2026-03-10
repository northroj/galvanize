import numpy as np
import matplotlib.pyplot as plt
import h5py


test_case = 1 # 1 =  protons on deuterium

if test_case == 1:
    csd_file = "plasma_slab_d_csd.h5"
    gfp2_file = "plasma_slab_d_gfp2.h5"
    x1, x2 = 1430000, 1510000
    slab = 20

def extract_tally(file_name):
    with h5py.File(file_name,"r") as f:
        tally = f["tallies"]["user_tallies"]["flux"]
        bins = np.array(tally["x_edges"][:])
        counts = np.array(tally["counts"][:])
    return bins, counts

csd_bins, csd_counts = extract_tally(csd_file)
gfp2_bins, gfp2_counts = extract_tally(gfp2_file)

plt.figure(1)
plt.plot(csd_bins[1:], csd_counts[0,:]/1000, color="red", label="csd")
plt.plot(gfp2_bins[1:], gfp2_counts[0,:]/1000, color="blue", label="gfp2")
#plt.xlim([x1, x2])
#plt.ylim(bottom=1e-4, top=None)
plt.xlabel("x [cm]")
plt.ylabel("Dose [MeV/cc]")
#plt.title(f"{slab} cm tungsten w/ 1.7e7 keV protons")
#plt.yscale("log")
plt.grid()
plt.legend()
plt.show()
