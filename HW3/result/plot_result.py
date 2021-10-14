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
r2 = np.genfromtxt('zeta100.txt')
r3 = np.genfromtxt('zeta1000.txt')
r4 = np.genfromtxt('zeta10000.txt')
r5 = np.genfromtxt('zeta100000.txt')
r6 = np.genfromtxt('zeta1000000.txt')
ref = np.genfromtxt('zeta_ref.txt')

plt.figure(1)
plt.plot(r1[:, 0], r1[:, 1], marker='o', linestyle='none', label='iter=10' )
plt.plot(r2[:, 0], r2[:, 1], marker='o', linestyle='none', label='iter=100' )
plt.plot(r3[:, 0], r3[:, 1], marker='o', linestyle='none', label='iter=1000' )
plt.plot(r4[:, 0], r4[:, 1], marker='o', linestyle='none', label='iter=10000' )
plt.plot(r5[:, 0], r5[:, 1], marker='o', linestyle='none', label='iter=100000' )
plt.plot(r6[:, 0], r6[:, 1], marker='o', linestyle='none', label='iter=1000000' )
plt.plot(ref[:, 0], ref[:, 1], label='reference (finite difference)')
plt.grid(True)
plt.xlabel(r'$t$')
plt.ylabel(r'$\zeta$')
plt.xlim(0, 3.2)
plt.ylim(0.85, 1.15)
plt.legend()
plt.savefig('plot_converge_zeta.png', dpi=150)
plt.savefig('plot_converge_zeta.eps', dpi=300)

plt.figure(2)
plt.plot(r1[:, 0], r1[:, 2], marker='o', linestyle='none', label='iter=10' )
plt.plot(r2[:, 0], r2[:, 2], marker='o', linestyle='none', label='iter=100' )
plt.plot(r3[:, 0], r3[:, 2], marker='o', linestyle='none', label='iter=1000' )
plt.plot(r4[:, 0], r4[:, 2], marker='o', linestyle='none', label='iter=10000' )
plt.plot(r5[:, 0], r5[:, 2], marker='o', linestyle='none', label='iter=100000' )
plt.plot(r6[:, 0], r6[:, 2], marker='o', linestyle='none', label='iter=1000000' )
plt.grid(True)
plt.xlabel(r'$t$')
plt.ylabel(r'$\theta$')
plt.xlim(0, 3.2)
plt.ylim(-0.05, 3.19)
plt.legend()
plt.savefig('plot_converge_theta.png', dpi=150)
plt.savefig('plot_converge_theta.eps', dpi=300)