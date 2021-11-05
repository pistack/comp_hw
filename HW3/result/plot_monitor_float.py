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
monitor1 = np.genfromtxt('monitor1f.txt')
monitor2 = np.genfromtxt('monitor2f.txt')
monitor3 = np.genfromtxt('monitor3f.txt')
monitor4 = np.genfromtxt('monitor3f_2.txt')
monitor5 = np.genfromtxt('monitor4f_2.txt')
converged_dist = monitor1[burn:, 1]
converged_dist2 = monitor2[burn:, 1]
converged_dist3 = monitor3[burn:, 1]
converged_dist4 = monitor4[burn:, 1]
converged_dist5 = monitor5[burn:, 1]
plt.figure(1)
plt.plot(monitor1[:, 0], monitor1[:, 1], marker='o', 
mfc='none', linestyle='none', color='red')
plt.xlabel('Number of Accepted moves')
plt.xscale('log')
plt.ylabel('Action')
plt.grid(True)
plt.savefig('plot_monitor1f.png', dpi=100)
plt.figure(2)
plt.plot(monitor2[:, 0], monitor2[:, 1], marker='o', 
mfc='none', linestyle='none', color='red')
plt.xlabel('Number of Accepted moves')
plt.xscale('log')
plt.ylabel('Action')
plt.grid(True)
plt.savefig('plot_monitor2f.png', dpi=100)
plt.figure(3)
plt.plot(monitor3[:, 0], monitor3[:, 1], marker='o', 
mfc='none', linestyle='none', color='red')
plt.xlabel('Number of Accepted moves')
plt.xscale('log')
plt.ylabel('Action')
plt.grid(True)
plt.savefig('plot_monitor3f.png', dpi=100)
plt.figure(4)
plt.plot(monitor4[:, 0], monitor4[:, 1], marker='o', 
mfc='none', linestyle='none', color='red')
plt.xlabel('Number of Accepted moves')
plt.xscale('log')
plt.ylabel('Action')
plt.grid(True)
plt.savefig('plot_monitor3f_2.png', dpi=100)
plt.figure(5)
plt.plot(monitor5[:, 0], monitor5[:, 1], marker='o', 
mfc='none', linestyle='none', color='red')
plt.xlabel('Number of Accepted moves')
plt.xscale('log')
plt.ylabel('Action')
plt.grid(True)
plt.savefig('plot_monitor4f.png', dpi=100)
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