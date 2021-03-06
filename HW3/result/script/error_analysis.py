import numpy as np
from scipy.integrate import simps as simpson;

zeta1f = np.genfromtxt('zeta1f.txt')
zeta2f = np.genfromtxt('zeta2f.txt')
zeta3f = np.genfromtxt('zeta3f.txt')
zeta3f_2 = np.genfromtxt('zeta3f_2.txt')
zeta4f_2 = np.genfromtxt('zeta4f_2.txt')
zeta6f_b = np.genfromtxt('zeta6f_b.txt')

zeta1 = np.genfromtxt('zeta1.txt')
zeta2 = np.genfromtxt('zeta2.txt')
zeta3 = np.genfromtxt('zeta3.txt')
zeta3_2 = np.genfromtxt('zeta3_2.txt')
zeta4_2 = np.genfromtxt('zeta4_2.txt')
zeta6_b = np.genfromtxt('zeta6_b.txt')

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
t_f = zeta1f[:, 0]

l_2_norm = np.sqrt(simpson(zeta_ref_sub[:, 1]**2, t))
error1f = error(zeta1f[:, 1:], zeta_ref_sub)
error2f = error(zeta2f[:, 1:], zeta_ref_sub)
error3f = error(zeta3f[:, 1:], zeta_ref_sub)
error3f_2 = error(zeta3f_2[:, 1:], zeta_ref_sub)
error4f_2 = error(zeta4f_2[:, 1:], zeta_ref_sub)
error6f_b = error(zeta6f_b[:, 1:], zeta_ref_sub)
l_2_error_1f = np.sqrt(simpson(error1f, t))/l_2_norm
l_2_error_2f = np.sqrt(simpson(error2f, t))/l_2_norm
l_2_error_3f = np.sqrt(simpson(error3f, t))/l_2_norm
l_2_error_3f_2 = np.sqrt(simpson(error3f_2, t))/l_2_norm
l_2_error_4f_2 = np.sqrt(simpson(error4f_2, t))/l_2_norm
l_2_error_6f_b = np.sqrt(simpson(error6f_b, t))/l_2_norm

error1 = error(zeta1[:, 1:], zeta_ref_sub)
error2 = error(zeta2[:, 1:], zeta_ref_sub)
error3 = error(zeta3[:, 1:], zeta_ref_sub)
error3_2 = error(zeta3_2[:, 1:], zeta_ref_sub)
error4_2 = error(zeta4_2[:, 1:], zeta_ref_sub)
error6_b = error(zeta6_b[:, 1:], zeta_ref_sub)
l_2_error_1 = np.sqrt(simpson(error1, t))/l_2_norm
l_2_error_2 = np.sqrt(simpson(error2, t))/l_2_norm
l_2_error_3 = np.sqrt(simpson(error3, t))/l_2_norm
l_2_error_3_2 = np.sqrt(simpson(error3_2, t))/l_2_norm
l_2_error_4_2 = np.sqrt(simpson(error4_2, t))/l_2_norm
l_2_error_6_b = np.sqrt(simpson(error6_b, t))/l_2_norm

print('single precision')
print('setup: 1')
print(f'n_f: 1 L_2 error: {l_2_error_1f:.5f}')
print(f'n_f: 2 L_2 error: {l_2_error_2f:.5f}')
print(f'n_f: 3 L_2 error: {l_2_error_3f:.5f}')
print('setup: 2')
print(f'n_f: 3 L_2 error: {l_2_error_3f_2:.5f}')
print(f'n_f: 4 L_2 error: {l_2_error_4f_2:.5f}')
print('beizer curve')
print(f'n: 6 L_2 error: {l_2_error_6f_b:.5f}')

print('double precision')
print('setup: 1')
print(f'n_f: 1 L_2 error: {l_2_error_1:.5f}')
print(f'n_f: 2 L_2 error: {l_2_error_2:.5f}')
print(f'n_f: 3 L_2 error: {l_2_error_3:.5f}')
print('setup: 2')
print(f'n_f: 3 L_2 error: {l_2_error_3_2:.5f}')
print(f'n_f: 4 L_2 error: {l_2_error_4_2:.5f}')
print('beizer curve')
print(f'n: 6 L_2 error: {l_2_error_6_b:.5f}')
