import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['text.usetex'] = True

def exp_decay(s, l, s_min):
    return l*np.exp(-l*(s-s_min))

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
monitor1 = np.fromfile('monitor1.bin', np.float64)
monitor2 = np.fromfile('monitor2.bin', np.float64)
monitor3 = np.fromfile('monitor3.bin', np.float64)
monitor4 = np.fromfile('monitor3_2.bin', np.float64)
monitor5 = np.fromfile('monitor4_2.bin', np.float64)
monitor6 = np.fromfile('monitor6_b.bin', np.float64)
monitor1_reshape = monitor1.reshape(int(monitor1.size/2), 2)
monitor2_reshape = monitor2.reshape(int(monitor2.size/2), 2)
monitor3_reshape = monitor3.reshape(int(monitor3.size/2), 2)
monitor4_reshape = monitor4.reshape(int(monitor4.size/2), 2)
monitor5_reshape = monitor5.reshape(int(monitor5.size/2), 2)
monitor6_reshape = monitor6.reshape(int(monitor6.size/2), 2)
converged_dist = monitor1_reshape[burn:, 0]
converged_dist2 = monitor2_reshape[burn:, 0]
converged_dist3 = monitor3_reshape[burn:, 0]
converged_dist4 = monitor4_reshape[burn:, 0]
converged_dist5 = monitor5_reshape[burn:, 0]
converged_dist6 = monitor6_reshape[burn:, 0]
prob1, action1 = np.histogram(converged_dist, bins=1000)
prob2, action2 = np.histogram(converged_dist2, bins=1000)
prob3, action3 = np.histogram(converged_dist3, bins=1000)
prob4, action4 = np.histogram(converged_dist4, bins=1000)
prob5, action5 = np.histogram(converged_dist5, bins=1000)
prob6, action6 = np.histogram(converged_dist6, bins=1000)
s1 = (action1[:1000]+action1[1:])/2
s2 = (action2[:1000]+action2[1:])/2
s3 = (action3[:1000]+action3[1:])/2
s4 = (action4[:1000]+action4[1:])/2
s5 = (action5[:1000]+action5[1:])/2
s6 = (action6[:1000]+action6[1:])/2

print(s1[np.argmax(prob1)])
print(s1[1]-s1[0])
print(s2[np.argmax(prob2)])
print(s2[1]-s2[0])
print(s3[np.argmax(prob3)])
print(s3[1]-s3[0])
print(s4[np.argmax(prob4)])
print(s4[1]-s4[0])
print(s5[np.argmax(prob5)])
print(s5[1]-s5[0])
print(s6[np.argmax(prob6)])
print(s6[1]-s6[0])

print(prob1[0]/np.max(prob1))
print(prob2[0]/np.max(prob2))
print(prob3[0]/np.max(prob3))
print(prob4[0]/np.max(prob4))
print(prob5[0]/np.max(prob5))
print(prob6[0]/np.max(prob6))


plt.figure(1)
plt.plot(np.arange(1, monitor1_reshape.shape[0]+1, 1), monitor1_reshape[:, 0], marker='o', 
mfc='none', linestyle='none', color='red')
plt.xlabel('Number of Accepted moves')
plt.xscale('log')
plt.ylabel('Action')
plt.grid(True)
plt.savefig('plot_monitor1.png', dpi=100)
plt.savefig('plot_monitor1_latex.png', dpi=300)
plt.figure(2)
plt.plot(np.arange(1, monitor2_reshape.shape[0]+1, 1), monitor2_reshape[:, 0], marker='o', 
mfc='none', linestyle='none', color='red')
plt.xlabel('Number of Accepted moves')
plt.xscale('log')
plt.ylabel('Action')
plt.grid(True)
plt.savefig('plot_monitor2.png', dpi=100)
plt.savefig('plot_monitor2_latex.png', dpi=300)
plt.figure(3)
plt.plot(np.arange(1, monitor3_reshape.shape[0]+1, 1), monitor3_reshape[:, 0], marker='o', 
mfc='none', linestyle='none', color='red')
plt.xlabel('Number of Accepted moves')
plt.xscale('log')
plt.ylabel('Action')
plt.grid(True)
plt.savefig('plot_monitor3.png', dpi=100)
plt.savefig('plot_monitor3_latex.png', dpi=300)
plt.figure(4)
plt.plot(np.arange(1, monitor4_reshape.shape[0]+1, 1), monitor4_reshape[:, 0], marker='o', 
mfc='none', linestyle='none', color='red')
plt.xlabel('Number of Accepted moves')
plt.xscale('log')
plt.ylabel('Action')
plt.grid(True)
plt.savefig('plot_monitor3_2.png', dpi=100)
plt.savefig('plot_monitor3_2_latex.png', dpi=300)
plt.figure(5)
plt.plot(np.arange(1, monitor5_reshape.shape[0]+1, 1), monitor5_reshape[:, 0], marker='o', 
mfc='none', linestyle='none', color='red')
plt.xlabel('Number of Accepted moves')
plt.xscale('log')
plt.ylabel('Action')
plt.grid(True)
plt.savefig('plot_monitor4_2.png', dpi=100)
plt.savefig('plot_monitor4_2_latex.png', dpi=300)
plt.figure(6)
plt.plot(np.arange(1, monitor6_reshape.shape[0]+1, 1), monitor6_reshape[:, 0], marker='o', 
mfc='none', linestyle='none', color='red')
plt.xlabel('Number of Accepted moves')
plt.xscale('log')
plt.ylabel('Action')
plt.grid(True)
plt.savefig('plot_monitor6_b.png', dpi=100)
plt.savefig('plot_monitor6_b_latex.png', dpi=300)
plt.figure(7)
plt.hist(converged_dist, bins = 1000, label=r'$n_f=1$, setup: 1')
plt.hist(converged_dist2, bins = 1000, label=r'$n_f=2$, setup: 1')
plt.hist(converged_dist3, bins = 1000, label=r'$n_f=3$, setup: 1')
plt.hist(converged_dist4, bins = 1000, label=r'$n_f=3$, setup: 2')
plt.hist(converged_dist5, bins = 1000, label=r'$n_f=4$, setup: 2')
plt.xlabel('Action')
plt.ylabel('Occurance')
plt.legend()
plt.savefig('plot_dist.png', dpi=100)
plt.savefig('plot_dist.eps', dpi=300)
plt.figure(8)
plt.hist(converged_dist4, bins = 1000, label=r'$n_f=3$, setup: 2')
plt.hist(converged_dist5, bins = 1000, label=r'$n_f=4$, setup: 2')
plt.xlabel('Action')
plt.ylabel('Occurance')
plt.legend()
plt.savefig('plot_dist2.png', dpi=100)
plt.savefig('plot_dist2.eps', dpi=300)
plt.figure(9)
plt.hist(converged_dist6, bins = 1000, label=r'$n=6$, setup: 2')
plt.xlabel('Action')
plt.ylabel('Occurance')
plt.legend()
plt.savefig('plot_dist_bezier.png', dpi=100)
plt.savefig('plot_dist_bezier.eps', dpi=300)