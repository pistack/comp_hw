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
monitor = np.genfromtxt('monitor4_2.txt')
converged_dist = monitor[burn:, 1]
plt.figure(1)
plt.plot(monitor[:, 0], monitor[:, 1], marker='o', 
mfc='none', linestyle='none', color='red')
plt.xlabel('Number of Accepted moves')
plt.xscale('log')
plt.ylabel('Action')
plt.grid(True)
plt.figure(2)
plt.hist(converged_dist, bins = 100)
plt.xlabel('Action')
plt.ylabel('Occurance')
plt.show()