'''
Homework03
Solve Kepler problem
Using Least action principle
'''

from typing import Tuple
import numpy as np


def sum_of_fourier(t: np.ndarray, c: np.ndarray,
                   num_fourier: int) -> Tuple[np.ndarray, np.ndarray]:

    '''
    Evaluate the sum of sine and cosine function
    weighted by c and its derivative

    Args:
     t: time range
     c: coefficient
     num_sine: number of sine function to add

   Returns:
    value and the derivative of the sum of sine and cosine
    function weighted by c

   Note:
    The periode of ith sine and cosine function is equal to 2*(tf-t0)/i
    for 0 < i < num_fourier+1.
   '''

    # time
    periode = 2*(t[-1]-t[0])
    n = t.shape[0]
    y = np.zeros(n)
    deriv = np.zeros(n)
    omega = 2*np.pi/periode
    for i in range(0, 2*num_fourier, 2):
        tmp = omega*(i/2+1)
        y = y + c[i]*np.sin(tmp*t)
        deriv = deriv + c[i]*tmp*np.cos(tmp*t)
        y = y + c[i]*np.cos(tmp*t)
        deriv = deriv - c[i]*tmp*np.sin(tmp*t)
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
              step: float) -> np.ndarray:

    '''
    Move init_guess by step

    Args:
     init_guess: initial guess
     step: step size

    Returns:
     initial guess moved by step
    '''
    move = init_guess + step*np.random.uniform(-1, 1,
                                               size=init_guess.shape[0])
    move[move > 1] = 1
    move[move < -1] = -1

    return move


def HW3(zeta_min: float, t0: float, n: int, num_fourier: int,
        num_iter: int,
        step: float,
        lamda: float) -> Tuple[float, np.ndarray, np.ndarray, np.ndarray]:

    '''
    Find the minimum action path using
    Markov Chain Monte Carlo Method

    Args:
     zeta_min: minimum value of zeta
     t0: initial time
     n: number of point to evaluate
     num_fourier: number of sine and cosine function
                  used for approximation of path
     num_iter: number of iteration
     step: size of step to move coefficients
     lamda: parameter to adapt step size

    Returns:
     Value of action of minimum path, time, zeta and theta.
    '''

    # define condition
    zeta_max = zeta_min/(2*zeta_min-1)
    a = 0.5*(zeta_min+zeta_max)
    t = np.linspace(0, np.pi*a**(1.5), n)

    # first guess
    c_zeta = np.zeros(2*num_fourier)
    c_theta = np.zeros(2*num_fourier)
    c_zeta[0] = 0.5
    c_theta[0] = 0.5
    c_zeta[1] = -0.5
    c_theta[1] = -0.5

    # adapt step
    adapt_step = step

    min_zeta, min_deriv_zeta = sum_of_fourier(t, c_zeta, num_fourier)
    min_theta, min_deriv_theta = sum_of_fourier(t, c_theta, num_fourier)

    scale_zeta = (zeta_max-zeta_min)/(min_zeta[-1]-min_zeta[0])
    add_zeta = zeta_min - scale_zeta*min_zeta[0]
    scale_theta = np.pi/(min_theta[-1]-min_theta[0])
    add_theta = -scale_theta*min_theta[0]

    min_zeta = scale_zeta*min_zeta + add_zeta
    min_deriv_zeta = scale_zeta*min_deriv_zeta
    min_theta = scale_theta*min_theta + add_theta
    min_deriv_theta = scale_theta*min_deriv_theta

    min_action = eval_action(t, min_zeta, min_theta,
                             min_deriv_zeta, min_deriv_theta)

    for i in range(num_iter):
        while True:
            tmp_c_zeta = move_step(c_zeta, adapt_step)
            tmp_c_theta = move_step(c_theta, adapt_step)
            tmp_zeta, tmp_deriv_zeta = \
                sum_of_fourier(t, tmp_c_zeta, num_fourier)
            tmp_theta, tmp_deriv_theta = \
                sum_of_fourier(t, tmp_c_theta, num_fourier)
            if not np.abs(tmp_zeta[0]-tmp_zeta[-1]) < 1e-8 or \
               np.abs(tmp_theta[0] - tmp_theta[-1]) < 1e-8:
                break

        scale_zeta = (zeta_max-zeta_min)/(tmp_zeta[-1]-tmp_zeta[0])
        add_zeta = zeta_min - scale_zeta*tmp_zeta[0]
        scale_theta = np.pi/(tmp_theta[-1]-tmp_theta[0])
        add_theta = -scale_theta*tmp_theta[0]

        tmp_zeta = scale_zeta*tmp_zeta + add_zeta
        tmp_deriv_zeta = scale_zeta*tmp_deriv_zeta
        tmp_theta = scale_theta*tmp_theta + add_theta
        tmp_deriv_theta = scale_theta*tmp_deriv_theta
        tmp_action = eval_action(t, tmp_zeta, tmp_theta,
                                 tmp_deriv_zeta, tmp_deriv_theta)

        if tmp_action < min_action:
            min_action = tmp_action
            min_zeta = tmp_zeta
            min_theta = tmp_theta

        if adapt_step > np.exp(-lamda*adapt_step*(tmp_action-min_action)):
            adapt_step = np.exp(-lamda*adapt_step*(tmp_action-min_action))
            print(adapt_step)

    return min_action, t+t0, min_zeta, min_theta
