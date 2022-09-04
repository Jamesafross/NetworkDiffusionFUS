import matplotlib.pyplot as plt
import matplotlib
import numpy as np
plt.rcParams["font.family"] = "Times New Roman"
import scipy.stats as scs
matplotlib.rcParams.update({'font.size': 32})
plt.rcParams["figure.figsize"] = (10,10)

sp_dist = np.load("sp_dist.npz")
sp_delay = sp_dist / 13
start_times = np.load("start_times.npz")
T = np.load("T.npz")
start_times_falff = np.load("start_times_falff.npz")
T_falff = np.load("T_falff.npz")
sp_dist_falff = np.load("sp_dist_falff.npz")


a = 120


plt.scatter(sp_dist,start_times,s=a)
plt.ylabel("ReHo - start time [mins]")
plt.xlabel("Shortest path length [mm]")
plt.title(r'$\rho$ = 0.88')
plt.savefig('start_time_vs_distance_ReHo.pdf')
plt.show()
plt.clf()

plt.scatter(sp_dist,T,s=a)
plt.ylabel("ReHo - peak effect (T)")
plt.xlabel("Shortest path length [mm]")
plt.title(r'$\rho$ = 0.84')
plt.savefig('T_vs_distance_ReHo.pdf')
plt.show()
plt.clf()


plt.scatter(sp_delay,start_times,s=a)
plt.ylabel("ReHo - start time [mins] ")
plt.xlabel("Minimum delay [msec]")
plt.title(r'$\rho$ = 0.88')
plt.savefig('start_time_vs_delay_ReHo.pdf')
plt.show()
plt.clf()

plt.scatter(sp_delay,T,s=a)
plt.ylabel("ReHo - peak effect (T)")
plt.xlabel("Minimum delay [msec]")
plt.title(r'$\rho$ = 0.84')
plt.savefig('T_vs_delay_ReHo.pdf')
plt.show()
plt.clf()


