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

r = np.genfromtxt('zeta1.txt')
r2 = np.genfromtxt('zeta2.txt')
r3 = np.genfromtxt('zeta3.txt')
r4 = np.genfromtxt('zeta3_2.txt')
r5 = np.genfromtxt('zeta4_2.txt')
ref = np.genfromtxt('zeta_ref.txt')


plt.figure(1)
plt.plot(r[:, 1]*np.cos(r[:, 2]), r[:, 1]*np.sin(r[:, 2]), marker='o',
markersize=3, linestyle='none',
label=r'$n_f=1$, setup: 1')
plt.plot(r2[:, 1]*np.cos(r2[:, 2]), r2[:, 1]*np.sin(r2[:, 2]), marker='o',
markersize=3, linestyle='none',
label=r'$n_f=2$, setup: 1')
plt.plot(r3[:, 1]*np.cos(r3[:, 2]), r3[:, 1]*np.sin(r3[:, 2]), marker='o',
markersize=3, linestyle='none',
label=r'$n_f=3$, setup: 1')
plt.plot(r4[:, 1]*np.cos(r4[:, 2]), r4[:, 1]*np.sin(r4[:, 2]), marker='o',
markersize=3, linestyle='none',
label=r'$n_f=3$, setup: 2')
plt.plot(r4[:, 1]*np.cos(r4[:, 2]), r4[:, 1]*np.sin(r4[:, 2]), marker='o',
markersize=3, linestyle='none',
label=r'$n_f=4$, setup: 2')
plt.plot(ref[:, 1]*np.cos(ref[:, 2]), ref[:, 1]*np.sin(ref[:, 2]), 
label=r'reference', color='black')
plt.grid(True)
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.xlim(-1.175, 0.95)
plt.legend(loc='upper right')
plt.savefig('plot_traj.png', dpi=100)
plt.savefig('plot_traj.eps', dpi=300)

plt.figure(2)
plt.plot(r[:, 0], r[:, 1], marker='o',
markersize=3, linestyle='none',
label=r'$n_f=1$, setup: 1')
plt.plot(r2[:, 0], r2[:, 1], marker='o',
markersize=3, linestyle='none',
label=r'$n_f=2$, setup: 1')
plt.plot(r3[:, 0], r3[:, 1], marker='o',
markersize=3, linestyle='none',
label=r'$n_f=3$, setup: 1')
plt.plot(r4[:, 0], r4[:, 1], marker='o',
markersize=3, linestyle='none',
label=r'$n_f=3$, setup: 2')
plt.plot(r5[:, 0], r5[:, 1], marker='o',
markersize=3, linestyle='none',
label=r'$n_f=4$, setup: 2')
plt.plot(ref[:, 0], ref[:, 1], label=r'referece',
color='black')
plt.grid(True)
plt.xlabel(r'$t$')
plt.ylabel(r'$\zeta$')
plt.xlim(0.0, 3.2)
plt.legend(loc='upper right')
plt.savefig('plot_zeta.png', dpi=100)
plt.savefig('plot_zeta.eps', dpi=300)

plt.figure(3)
plt.plot(r[:, 0], r[:, 2], marker='o',
markersize=3, linestyle='none',
label=r'$n_f=1$, setup: 1')
plt.plot(r2[:, 0], r2[:, 2], marker='o',
markersize=3, linestyle='none',
label=r'$n_f=2$, setup: 1')
plt.plot(r3[:, 0], r3[:, 2], marker='o',
markersize=3, linestyle='none',
label=r'$n_f=3$, setup: 1')
plt.plot(r4[:, 0], r4[:, 2], marker='o',
markersize=3, linestyle='none',
label=r'$n_f=3$, setup: 2')
plt.plot(r5[:, 0], r5[:, 2], marker='o',
markersize=3, linestyle='none',
label=r'$n_f=4$, setup: 2')
plt.plot(ref[:, 0], ref[:, 2], label=r'reference',
color='black')
plt.grid(True)
plt.xlabel(r'$t$')
plt.ylabel(r'$\theta$')
plt.xlim(0.0, 3.2)
plt.legend(loc='upper right')
plt.savefig('plot_theta.png', dpi=100)
plt.savefig('plot_theta.eps', dpi=300)