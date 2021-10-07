import numpy as np
import matplotlib.pyplot as plt
from HW1 import HW1


# Change of variable
# Final equation
# t = 
# 8/sqrt(alpha) int 0 to sqrt(zeta-9/10) sqrt(9/8-9/10-x^2) dx +
# 10/sqrt(alpha) int sqrt(9/8-zeta) to sqrt(9/8-9/10) sqrt(9/8-9/10-x^2) dx

def f(x):
    return np.sqrt(9/8-9/10-x**2)


def HW2(zeta_min, zeta_max, n):
    zeta = np.linspace(zeta_min, zeta_max, n)
    t = np.zeros(n)
    x1 = np.sqrt(zeta-9/10)
    x2 = np.sqrt(9/8-zeta)
    fx1 = f(x1)
    fx2 = f(x2)
    t[0] = 0
    for i in range(n-1):
        t[i+1] = t[i] + 8/np.sqrt(0.987654)*(x1[i+1]-x1[i])*fx1[i] - \
            10/np.sqrt(0.987654)*(x2[i+1]-x2[i])*fx2[i]

    return t, zeta


if __name__ == "__main__":

    t_ref, zeta_ref = HW1(0, 3.2, 1000, 0.9, 0)
    t_3, zeta_3 = HW2(0.9, 9/8, 200)
    plt.plot(t_ref, zeta_ref, label='HW1')
    plt.plot(t_3, zeta_3, marker='o', markersize=3,
             linestyle='none', label='HW2')
    plt.legend()
    plt.show()
        
