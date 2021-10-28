import numpy as np
from scipy.optimize import curve_fit
import matplotlib
import matplotlib.pyplot as plt
# matplotlib.rcParams['text.usetex'] = True

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

task_lst = [5, 10, 100, 1000, 10000]

exact = np.pi*((0.9+1.125)/2.0)**(1.5) # t_f
error = np.zeros(5)

r1 = np.genfromtxt('zeta5.txt')
r2 = np.genfromtxt('zeta10.txt')
r3 = np.genfromtxt('zeta100.txt')
r4 = np.genfromtxt('zeta1000.txt')
r5 = np.genfromtxt('zeta10000.txt')


error[0] = np.abs(exact-r1[-1, 0])
error[1] = np.abs(exact-r2[-1, 0])
error[2] = np.abs(exact-r3[-1, 0])
error[3] = np.abs(exact-r4[-1, 0])
error[4] = np.abs(exact-r5[-1, 0])


plt.plot(task_lst, error, marker='o')
plt.grid(True)
plt.xlabel(r'$n$')
plt.ylabel('error')
plt.xscale('log')
plt.yscale('log')
plt.savefig('plot_error.png', dpi=100)
plt.savefig('plot_error.eps', dpi=300)

log_n = np.log(np.array(task_lst))
log_error = np.log(error)

def f(x, a, b):
    return a*x+b
    
popt, pconv = curve_fit(f, log_n, log_error)
print('a \t b')
print(f'{popt[0]} \t {popt[1]}')