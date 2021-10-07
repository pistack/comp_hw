'''
HW1
python code for homework1 in computer1 class in Yonsei University
Using explicit Euler Method to solve Kepler problem
'''

import numpy as np
import matplotlib.pyplot as plt


def HW1(t0: float, t1: float, n: int,
        y0: float, y0p: float) -> np.ndarray:

    '''
    Homework1:
    Solve :math:`\\ddots{\\zeta} = \\zeta^{-3} - \\zeta^{-2}`
    with initial condtion
     * y0 = y0
     * y'0 = y0p

    Using explict Euler Method

    Args:
     t0: intial time
     t1: finial time
      n: number of grid points between t0 and t1
     y0: initial condition for y
     y0p: initial condition for y'

    Returns:
     Solution of 2nd order ODE computed in n grid points between
     t0 and t1
    '''

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
