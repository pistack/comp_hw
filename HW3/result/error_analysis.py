import numpy as np
from scipy.integrate import simpson;

zeta1 = np.genfromtxt('zeta1.txt')
zeta2 = np.genfromtxt('zeta2.txt')
zeta3 = np.genfromtxt('zeta3.txt')
zeta3_2 = np.genfromtxt('zeta3_2.txt')
zeta4_2 = np.genfromtxt('zeta4_2.txt')
zeta_ref = np.genfromtxt('zeta_ref.txt')
zeta_ref_sub  = np.zeros((101, 2))
offset = 1000
for i in range(101):
    zeta_ref_sub[i, :] = zeta_ref[offset*i, 1:]

def error(p, p_ref):
    x = p[:, 0]*np.cos(p[:, 1])
    y = p[:, 0]*np.sin(p[:, 1])
    x_ref = p_ref[:, 0]*np.cos(p_ref[:, 1])
    y_ref = p_ref[:, 0]*np.sin(p_ref[:, 1])
    return (x-x_ref)**2 + (y-y_ref)**2

t = zeta1[:, 0]
error1 = error(zeta1[:, 1:], zeta_ref_sub)

error2 = error(zeta2[:, 1:], zeta_ref_sub)
error3 = error(zeta3[:, 1:], zeta_ref_sub)
error3_2 = error(zeta3_2[:, 1:], zeta_ref_sub)
error4_2 = error(zeta4_2[:, 1:], zeta_ref_sub)
l_2_error_1 = np.sqrt(simpson(error1, t))
l_2_error_2 = np.sqrt(simpson(error2, t))
l_2_error_3 = np.sqrt(simpson(error3, t))
l_2_error_3_2 = np.sqrt(simpson(error3_2, t))
l_2_error_4_2 = np.sqrt(simpson(error4_2, t))
print('setup: 1')
print(f'n_f: 1 L_2 error: {l_2_error_1:.14f}')
print(f'n_f: 2 L_2 error: {l_2_error_2:.14f}')
print(f'n_f: 3 L_2 error: {l_2_error_3:.14f}')
print('setup: 2')
print(f'n_f: 3 L_2 error: {l_2_error_3_2:.14f}')
print(f'n_f: 4 L_2 error: {l_2_error_4_2:.14f}')
