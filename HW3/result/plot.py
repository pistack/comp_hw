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

r = np.genfromtxt('zeta_nf1.txt')
r2 = np.genfromtxt('zeta_nf2.txt')
r3 = np.genfromtxt('zeta_nf4.txt')
ref = np.genfromtxt('zeta_ref.txt')


plt.figure(1)
plt.plot(r[:, 1]*np.cos(r[:, 2]), r[:, 1]*np.sin(r[:, 2]), marker='o',
markersize=3, linestyle='none',
label='nf=1')
plt.plot(r2[:, 1]*np.cos(r2[:, 2]), r2[:, 1]*np.sin(r2[:, 2]), marker='o',
markersize=3, linestyle='none',
label='nf=2')
plt.plot(r3[:, 1]*np.cos(r3[:, 2]), r3[:, 1]*np.sin(r3[:, 2]), marker='o',
markersize=3, linestyle='none',
label='nf=4')
plt.plot(ref[:, 1]*np.cos(ref[:, 2]), ref[:, 1]*np.sin(ref[:, 2]), color='black')
plt.grid(True)
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
# plt.xlim(-1.15, 0.95)
# plt.ylim(-1.1, 1.1)
plt.legend(loc='upper right')

plt.figure(2)
plt.plot(r[:, 0], r[:, 1], marker='o',
markersize=3, linestyle='none',
label='nf=1')
plt.plot(r2[:, 0], r2[:, 1], marker='o',
markersize=3, linestyle='none',
label='nf=2')
plt.plot(r3[:, 0], r3[:, 1], marker='o',
markersize=3, linestyle='none',
label='nf=4')
plt.plot(ref[:, 0], ref[:, 1], label='ref',
color='black')
plt.grid(True)
plt.xlabel(r'$t$')
plt.ylabel(r'$\zeta$')
# plt.xlim(-1.15, 0.95)
# plt.ylim(-1.1, 1.1)
plt.legend(loc='upper right')

plt.figure(3)
plt.plot(r[:, 0], r[:, 2], marker='o',
markersize=3, linestyle='none',
label='nf=1')
plt.plot(r2[:, 0], r2[:, 2], marker='o',
markersize=3, linestyle='none',
label='nf=2')
plt.plot(r3[:, 0], r3[:, 2], marker='o',
markersize=3, linestyle='none',
label='nf=4')
plt.plot(ref[:, 0], ref[:, 2], label='ref',
color='black')
plt.grid(True)
plt.xlabel(r'$t$')
plt.ylabel(r'$\theta$')
# plt.xlim(-1.15, 0.95)
# plt.ylim(-1.1, 1.1)
plt.legend(loc='upper right')
plt.show()