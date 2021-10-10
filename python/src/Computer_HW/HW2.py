'''
HW2
python code for homework2 in computer1 class in Yonsei University
Using numerical integration to solve Kepler problem
'''
from typing import Tuple
import numpy as np


def HW2(zeta_min: float, t0: float,
        n: int) -> Tuple[np.ndarray, np.ndarray]:

    '''
    Integrate
    :math:`\\zeta\'(-a(\\zeta\'-\\zeta_{min})(\\zeta_{max}-\\zeta\'))^{-1/2}`
    from :math:`\\zeta_{min}` to :math:`\\zeta`
    To remove singularity in the integrand, I separates the integrand
    and changes
    :math:`\\zeta\'` to :math:`x=(\\zeta\'-\\zeta_{min})^{1/2}`
    for the first part and
    to :math:`x=(\\zeta_{max}-\\zeta\')^{1/2}` for the second part.

    Args:
     zeta_min: minimum value of zeta
     t0: initial time
     n: number of points between zeta_min and zeta_max

    Returns:
     Integration evaluated in n points between zeta_min and zeta_max

    Notes:
     * :math:`a = \\zeta_{min}^{-2}-2\\zeta_{min}^{-1}`
     * :math:`a = \\zeta_{max}^{-2}-2\\zeta_{max}^{-1}`
    '''

    tmp = 1/zeta_min
    a = tmp*(tmp-2)
    tmp = 1.0+np.sqrt(1.0+a)
    zeta_max = -tmp/a
    tmp = zeta_max - zeta_min
    tmp = 1/tmp
    c1 = 2*zeta_min*tmp
    c2 = 2*zeta_max*tmp
    tmp = np.sqrt(-a)
    c1 = c1/tmp
    c2 = c2/tmp

    zeta = np.linspace(zeta_min, zeta_max, n+1)
    t = np.zeros(n+1)
    x1 = np.sqrt(zeta-zeta_min)
    x2 = np.sqrt(zeta_max-zeta)
    t[0] = 0
    for i in range(n):
        t[i+1] = t[i] + c1*(x1[i+1]-x1[i])*x2[i] - \
            c2*(x2[i+1]-x2[i])*x1[i]

    return t+t0, zeta
