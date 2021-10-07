import numpy as np
import matplotlib.pyplot as plt


def HW1(t0, t1, n, y0, y0p):
    t = np.linspace(t0, t1, n)
    
    zeta = np.zeros(n)
    zeta[0] = y0
    zeta[1] = 1/n*(t1-t0)*y0p + zeta[0]

    for i in range(n-2):
        zeta[i+2] = ((t1-t0)/n)**2*(1/zeta[i]**3 - 1/zeta[i]**2) + \
            2*zeta[i+1] - zeta[i]

    return t, zeta


if __name__ == "__main__":
    t, zeta = HW1(0, 10, 1000, 0.9, 0)
    plt.plot(t, zeta)
    plt.show()
    
