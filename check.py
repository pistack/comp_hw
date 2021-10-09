import subprocess
import numpy as np
import matplotlib.pyplot as plt
from HW1 import HW1

if __name__ == "__main__":

    t_ref, zeta_ref = HW1(0, 3.2, 1000, 0.9, 0)
    f = open('tmp.txt', 'w')
    subprocess.run(['./a.out'], stdout=f)
    f.close()
    A = np.genfromtxt('tmp.txt')
    plt.plot(t_ref, zeta_ref, label='HW python')
    plt.plot(A[:, 0], A[:, 1], marker='o', markersize=3,
             linestyle='none', label='HW cpp')
    plt.legend()
    plt.show()
