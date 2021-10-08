'''
HW2
python code for homework2 in computer1 class in Yonsei University
Using numerical integration to solve Kepler problem
'''
from typing import Tuple
import numpy as np
import matplotlib.pyplot as plt
from HW1 import HW1


def HW2(zeta_min: float, zeta_max: float,
        n: int) -> Tuple[np.ndarray, np.ndarray]:

    '''
    Integrate :math:`\\zeta\'(0.987654(\\zeta\'-0.9)(9/8-\\zeta\'))^{-1/2}`
    from 0.9 to :math:`\\zeta`
    To remove singularity in the integrand, I separates the integrand
    and changes
    :math:`\\zeta\'` to :math:`x=(\\zeta\'-0.9)^{1/2}` for the first part and
    to :math:`x=(9/8-\\zeta\')^{1/2}` for the second part.

    Args:
     zeta_min: minimum value of zeta
     zeta_max: maximum value of zeta
     n: number of points between zeta_min and zeta_max

    Returns:
     Integration evaluated in n points between zeta_min and zeta_max
    '''

    zeta = np.linspace(zeta_min, zeta_max, n+1)
    t = np.zeros(n+1)
    x1 = np.sqrt(zeta-zeta_min)
    x2 = np.sqrt(zeta_max-zeta)
    t[0] = 0
    for i in range(n):
        t[i+1] = t[i] + 8/np.sqrt(0.987654)*(x1[i+1]-x1[i])*x2[i] - \
            10/np.sqrt(0.987654)*(x2[i+1]-x2[i])*x1[i]

    return t, zeta


if __name__ == "__main__":

    t_ref, zeta_ref = HW1(0, 3.2, 1000, 0.9, 0)
    t, zeta = HW2(0.9, 9/8, 10)
    plt.plot(t_ref, zeta_ref, label='HW1')
    plt.plot(t, zeta, marker='o', markersize=3,
             linestyle='none', label='HW2')
    plt.legend()
    plt.show()
