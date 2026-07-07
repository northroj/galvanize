import numpy as np
import matplotlib.pyplot as plt
import h5py


s_s_file = "spitzer_spitzer.h5"
c_s_file = "cutoff_spitzer.h5"
c2_s_file = "cutoff2_spitzer.h5"
s_r_file = "spitzer_rutherford.h5"


def extract_tally(file_name):
    with h5py.File(file_name,"r") as f:
        tally = f["tallies"]["user_tallies"]["flux"]
        bins = np.array(tally["x_edges"][:])
        counts = np.array(tally["counts"][:])
    return bins, counts

s_s_bins, s_s_counts = extract_tally(s_s_file)
c_s_bins, c_s_counts = extract_tally(c_s_file)
c2_s_bins, c2_s_counts = extract_tally(c2_s_file)
s_r_bins, s_r_counts = extract_tally(s_r_file)
test_bins, test_counts = extract_tally("test.h5")


plt.figure(1)
plt.plot(s_s_bins[1:], s_s_counts[0,:], color="red", label="spitzer(CSD)-spitzer(CSD)")
#plt.plot(c_s_bins[1:], c_s_counts[0,:], color="blue", label="cutoff-spitzer")
#plt.plot(c2_s_bins[1:], c2_s_counts[0,:], color="black", label="cutoff(GFP2)-spitzer(CSD)")
#plt.plot(s_r_bins[1:], s_r_counts[0,:], color="green", label="spitzer-rutherford")
plt.plot(test_bins[1:], test_counts[0,:], color="blue", label = "test")

#plt.xlim([x1, x2])
#plt.ylim(bottom=1e-4, top=None)
plt.xlabel("x [cm]")
plt.ylabel("Energy deposition")
plt.title(f"10 MeV protons on 5 g/cc deuterium w/ \nvarious <ion>-<electron> soft scattering models")
#plt.yscale("log")
plt.grid()
plt.legend()
plt.show()
