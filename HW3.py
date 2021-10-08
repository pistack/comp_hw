'''
Homework03
Solve Kepler problem
Using Least action principle
'''

from typing import Tuple
import numpy as np
import matplotlib.pyplot as plt
from HW1 import HW1


def sum_of_sine(t: np.ndarray, c: np.ndarray,
                num_sine: int) -> Tuple[np.ndarray, np.ndarray]:

    '''
    Evaluate the sum of sine function weighted by c and its derivative

    Args:
     t: time range
     c: coefficient
     num_sine: number of sine function to add

   Returns:
    value and the derivative of the sum of sine function weighted by c

   Note:
    The periode of ith sine function is equal to 4*(tf-t0)/i
    for 0 < i < num_sine+1.
   '''

    # time
    periode = 4*(t[-1]-t[0])
    n = t.shape[0]
    y = np.zeros(n)
    deriv = np.zeros(n)
    omega = 2*np.pi/periode
    for i in range(num_sine):
        tmp = omega*(i+1)
        y = y + c[i]*np.sin(tmp*t)
        deriv = deriv + c[i]*tmp*np.cos(tmp*t)
    return y, deriv


def eval_action(t: np.ndarray, zeta: np.ndarray,
                theta: np.ndarray,
                deriv_zeta: np.ndarray, deriv_theta: np.ndarray) -> float:

    '''
    Evaluate action of given path

    Args:
     t: time
     zeta: zeta of path
     theta: theta of path
     deriv_zeta: derivative of zeta
     deriv_theta: derivative of theta

    Returns:
     The action of given path
    '''

    zeta_inv = 1/zeta
    dt = np.zeros_like(t)
    dt[0] = 0
    dt[1:] = t[1:]-t[:-1]
    action = (0.5*(deriv_zeta**2+(zeta*deriv_theta)**2)+np.abs(zeta_inv))*dt

    return np.sum(action)


def move_step(init_guess: np.ndarray,
              lb: np.ndarray, ub: np.ndarray,
              step: float) -> np.ndarray:

    '''
    Move init_guess by step

    Args:
     init_guess: initial guess
     weight: weight
     lb: lower bound
     ub: upper bound
     step: step size

    Returns:
     initial guess moved by step
    '''

    c = init_guess + step*(np.random.rand(init_guess.shape[0])-0.5)
    c[c < lb] = lb[c < lb]
    c[c > ub] = ub[c > ub]

    return c


def HW3(zeta_min: float, t0: float, n: int, num_sine: int,
        num_iter: int) -> Tuple[float, np.ndarray, np.ndarray, np.ndarray]:

    '''
    Find the minimum action path using
    Markov Chain Monte Carlo Method

    Args:
     zeta_min: minimum value of zeta
     t0: initial time
     n: number of point to evaluate
     num_sine: number of sine function used for approximation of path
     num_iter: number of iteration

    Returns:
     Value of action of given path, time, zeta and theta.
    '''

    # define condition
    zeta_max = zeta_min/(2*zeta_min-1)
    a = 0.5*(zeta_min+zeta_max)
    t = np.linspace(0, np.pi*a**(1.5), n)
    lb = -0.5*np.ones(num_sine)
    ub = 0.5*np.ones(num_sine)
    weight = 1/np.linspace(1, num_sine, 1)

    # define step
    step = (1/num_iter)**(1/num_sine)

    # first guess
    c_zeta = (np.random.rand(num_sine)-0.5)*weight
    c_theta = (np.random.rand(num_sine)-0.5)*weight

    min_zeta, min_deriv_zeta = sum_of_sine(t, c_zeta, num_sine)
    min_theta, min_deriv_theta = sum_of_sine(t, c_theta, num_sine)

    while (min_zeta[-1] == 0) or (min_theta[-1] == 0):
        c_zeta = move_step(c_zeta, lb, ub, step)
        c_theta = move_step(c_theta, lb, ub, step)
        min_zeta, min_deriv_zeta = sum_of_sine(t, c_zeta, num_sine)
        min_theta, min_deriv_theta = sum_of_sine(t, c_theta, num_sine)

    scale_zeta = (zeta_max-zeta_min)/min_zeta[-1]
    scale_theta = np.pi/min_theta[-1]
    min_zeta = scale_zeta*min_zeta + zeta_min
    min_deriv_zeta = scale_zeta*min_deriv_zeta
    min_theta = scale_theta*min_theta
    min_deriv_theta = scale_theta*min_deriv_theta

    min_action = eval_action(t, min_zeta, min_theta,
                             min_deriv_zeta, min_deriv_theta)

    for i in range(num_iter):
        c_zeta = move_step(c_zeta, lb, ub, step)
        c_theta = move_step(c_theta, lb, ub, step)
        tmp_zeta, tmp_deriv_zeta = sum_of_sine(t, c_zeta, num_sine)
        tmp_theta, tmp_deriv_theta = sum_of_sine(t, c_theta, num_sine)

        while (min_zeta[-1] == 0) or (min_theta[-1] == 0):
            c_zeta = move_step(c_zeta, lb, ub, step)
            c_theta = move_step(c_theta, lb, ub, step)
            tmp_zeta, tmp_deriv_zeta = sum_of_sine(t, c_zeta, num_sine)
            tmp_theta, tmp_deriv_theta = sum_of_sine(t, c_theta, num_sine)

        scale_zeta = (zeta_max-zeta_min)/tmp_zeta[-1]
        scale_theta = np.pi/tmp_theta[-1]
        tmp_zeta = scale_zeta*tmp_zeta + zeta_min
        tmp_deriv_zeta = scale_zeta*tmp_deriv_zeta
        tmp_theta = scale_theta*tmp_theta
        tmp_deriv_theta = scale_theta*tmp_deriv_theta

        tmp_action = eval_action(t, tmp_zeta, tmp_theta,
                                 tmp_deriv_zeta, tmp_deriv_theta)

        if tmp_action < min_action:
            min_action = tmp_action
            min_zeta = tmp_zeta
            min_theta = tmp_theta
        else:
            pass

    return min_action, t+t0, np.abs(min_zeta), np.abs(min_theta)


if __name__ == "__main__":

    t_ref, zeta_ref = HW1(0, 3.2, 100, 0.9, 0)
    action, t, zeta, theta = HW3(0.9, 0, 15, 5, 2**20)
    print(action)
    plt.plot(t_ref, zeta_ref, label='HW1')
    plt.plot(t, zeta, marker='o', markersize=3,
             linestyle='none', label='HW3')
    plt.legend()
    plt.show()
