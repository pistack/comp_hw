#Error Analysis pick 10 pts
#t=1, 2, 3, 4, 5, 6, 7, 8, 9, 10
# Assume n=10^5 is exact!
import numpy as np
from scipy.optimize import curve_fit
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

task_list = [10, 20, 100, 1000, 10000]
t=list(range(1, 11, 1))
EA = np.zeros((10, 5))
for i in range(5):
    r = np.genfromtxt(f'zeta{task_list[i]}f.txt')
    ref = np.genfromtxt(f'zeta{task_list[i]}.txt')
    offset = int(task_list[i]/10)
    offset_ref = offset
    for j in range(1, 11, 1):
        EA[j-1, i] = \
            (r[j*offset, 1]*np.cos(r[j*offset, 2]) -
            ref[j*offset_ref, 1]*np.cos(ref[j*offset_ref, 2]))**2 + \
                (r[j*offset, 1]*np.sin(r[j*offset, 2]) -
                ref[j*offset_ref, 1]*np.sin(ref[j*offset_ref, 2]))**2

EA = np.sqrt(EA)

for j in range(5):
    print(f'n={task_list[j]}')
    str = ''
    for i in range(1, 11, 1):
        print(f't={i}: {EA[i-1, j]}')



