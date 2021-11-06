import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['text.usetex'] = True

SMALL_SIZE = 10
MEDIUM_SIZE = 15
BIGGER_SIZE = 20

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

burn = 10000
monitor1 = np.fromfile('monitor1f.bin', np.float32)
monitor2 = np.fromfile('monitor2f.bin', np.float32)
monitor3 = np.fromfile('monitor3f.bin', np.float32)
monitor4 = np.fromfile('monitor3f_2.bin', np.float32)
monitor5 = np.fromfile('monitor4f_2.bin', np.float32)
monitor1_reshape = monitor1.reshape(int(monitor1.size/2), 2)
monitor2_reshape = monitor2.reshape(int(monitor2.size/2), 2)
monitor3_reshape = monitor3.reshape(int(monitor3.size/2), 2)
monitor4_reshape = monitor4.reshape(int(monitor4.size/2), 2)
monitor5_reshape = monitor5.reshape(int(monitor5.size/2), 2)
converged_dist = monitor1_reshape[burn:, 0]
converged_dist2 = monitor2_reshape[burn:, 0]
converged_dist3 = monitor3_reshape[burn:, 0]
converged_dist4 = monitor4_reshape[burn:, 0]
converged_dist5 = monitor5_reshape[burn:, 0]
plt.figure(1)
plt.plot(np.arange(1, monitor1_reshape.shape[0]+1, 1), monitor1_reshape[:, 0], marker='o', 
mfc='none', linestyle='none', color='red')
plt.xlabel('Number of Accepted moves')
plt.xscale('log')
plt.ylabel('Action')
plt.grid(True)
plt.savefig('plot_monitor1f.png', dpi=100)
plt.figure(2)
plt.plot(np.arange(1, monitor2_reshape.shape[0]+1, 1), monitor2_reshape[:, 0], marker='o', 
mfc='none', linestyle='none', color='red')
plt.xlabel('Number of Accepted moves')
plt.xscale('log')
plt.ylabel('Action')
plt.grid(True)
plt.savefig('plot_monitor2f.png', dpi=100)
plt.figure(3)
plt.plot(np.arange(1, monitor3_reshape.shape[0]+1, 1), monitor3_reshape[:, 0], marker='o', 
mfc='none', linestyle='none', color='red')
plt.xlabel('Number of Accepted moves')
plt.xscale('log')
plt.ylabel('Action')
plt.grid(True)
plt.savefig('plot_monitor3f.png', dpi=100)
plt.figure(4)
plt.plot(np.arange(1, monitor4_reshape.shape[0]+1, 1), monitor4_reshape[:, 0], marker='o', 
mfc='none', linestyle='none', color='red')
plt.xlabel('Number of Accepted moves')
plt.xscale('log')
plt.ylabel('Action')
plt.grid(True)
plt.savefig('plot_monitor3f_2.png', dpi=100)
plt.figure(5)
plt.plot(np.arange(1, monitor5_reshape.shape[0]+1, 1), monitor5_reshape[:, 0], marker='o', 
mfc='none', linestyle='none', color='red')
plt.xlabel('Number of Accepted moves')
plt.xscale('log')
plt.ylabel('Action')
plt.grid(True)
plt.savefig('plot_monitor4f_2.png', dpi=100)
plt.figure(6)
plt.hist(converged_dist, bins = 1000, label=r'$n_f=1$, setup: 1')
plt.hist(converged_dist2, bins = 1000, label=r'$n_f=2$, setup: 1')
plt.hist(converged_dist3, bins = 1000, label=r'$n_f=3$, setup: 1')
plt.hist(converged_dist4, bins = 1000, label=r'$n_f=3$, setup: 2')
plt.hist(converged_dist5, bins = 1000, label=r'$n_f=4$, setup: 2')
plt.xlabel('Action')
plt.ylabel('Occurance')
plt.legend()
plt.savefig('plot_distf.png', dpi=100)
plt.figure(7)
plt.hist(converged_dist4, bins = 1000, label=r'$n_f=3$, setup: 2')
plt.hist(converged_dist5, bins = 1000, label=r'$n_f=4$, setup: 2')
plt.xlabel('Action')
plt.ylabel('Occurance')
plt.legend()
plt.savefig('plot_dist2f.png', dpi=100)