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

task_lst = [10, 20, 100, 1000, 10000]

r1 = np.genfromtxt('zeta10.txt')
r2 = np.genfromtxt('zeta20.txt')
r3 = np.genfromtxt('zeta100.txt')
r4 = np.genfromtxt('zeta1000.txt')
r5 = np.genfromtxt('zeta10000.txt')
r6 = np.genfromtxt('zeta100000.txt')

plt.figure(1)
plt.plot(r1[:, 0], r1[:, 1], marker='o', linestyle='none', label='n=10' )
plt.plot(r2[:, 0], r2[:, 1], marker='o', linestyle='none', label='n=20' )
plt.plot(r3[:, 0], r3[:, 1], label='n=100' )
plt.plot(r4[:, 0], r4[:, 1], label='n=1000' )
plt.plot(r5[:, 0], r5[:, 1], label='n=10000' )
plt.plot(r6[:, 0], r6[:, 1], label='n=100000' )
plt.grid(True)
plt.xlabel(r'$t$')
plt.ylabel(r'$\zeta$')
plt.xlim(0, 10)
plt.ylim(0.85, 1.15)
plt.legend()
plt.savefig('plot_converge_zeta_double.png', dpi=100)
plt.savefig('plot_converge_zeta_double.eps', dpi=300)

plt.figure(2)
plt.plot(r1[:, 0], r1[:, 2], marker='o', linestyle='none', label='n=10' )
plt.plot(r2[:, 0], r2[:, 2], marker='o', linestyle='none', label='n=20' )
plt.plot(r3[:, 0], r3[:, 2], label='n=100' )
plt.plot(r4[:, 0], r4[:, 2], label='n=1000' )
plt.plot(r5[:, 0], r5[:, 2], label='n=10000' )
plt.plot(r6[:, 0], r6[:, 2], label='n=100000' )
plt.grid(True)
plt.xlabel(r'$t$')
plt.ylabel(r'$\theta$')
plt.xlim(0, 10)
plt.ylim(0, 10)
plt.legend()
plt.savefig('plot_converge_theta_double.png', dpi=100)
plt.savefig('plot_converge_theta_double.eps', dpi=300)

plt.figure(3)
plt.plot(r1[:, 1]*np.cos(r1[:, 2]), r1[:, 1]*np.sin(r1[:, 2]), 
marker='o', linestyle='none', label='n=10')
plt.plot(r2[:, 1]*np.cos(r2[:, 2]), r2[:, 1]*np.sin(r2[:, 2]), 
marker='o', linestyle='none', label='n=20')
plt.plot(r3[:, 1]*np.cos(r3[:, 2]), r3[:, 1]*np.sin(r3[:, 2]), 
label='n=100')
plt.plot(r4[:, 1]*np.cos(r4[:, 2]), r4[:, 1]*np.sin(r4[:, 2]), 
label='n=1000')
plt.plot(r5[:, 1]*np.cos(r5[:, 2]), r5[:, 1]*np.sin(r5[:, 2]), 
label='n=10000')
plt.plot(r6[:, 1]*np.cos(r6[:, 2]), r6[:, 1]*np.sin(r6[:, 2]), 
label='n=100000')
plt.grid(True)
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.xlim(-1.15, 0.95)
plt.ylim(-1.1, 1.1)
plt.legend(loc='upper right')
plt.savefig('plot_converge_traj_double.png', dpi=100)
plt.savefig('plot_converge_traj_double.eps', dpi=300)