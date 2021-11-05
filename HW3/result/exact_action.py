import numpy as np
from scipy.integrate import simpson;

def deriv_mid(t, p):
    dp = np.zeros_like(p)
    for i in range(1, t.size-1):
        h = t[i+1] - t[i-1]
        dp[i, :] = (p[i+1, :] - p[i-1, :])/h
    return dp

def kepler_lag(p, dp):
    return 0.5*(dp[:, 0]**2.0+(1/p[:, 0])**2.0)+1/p[:, 0]


ref = np.genfromtxt('zeta_ref.txt')

t_ref = ref[:, 0]
p_ref = ref[:, 1:3]
dp_ref = deriv_mid(t_ref, p_ref)
lag_ref = kepler_lag(p_ref, dp_ref)
action = simpson(lag_ref, t_ref)
print(f'{action:.14f}')
