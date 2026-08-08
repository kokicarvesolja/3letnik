#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import scipy.special as ss
import scipy.integrate as si

# Constants

g = 9.81
L = 1.0
v0 = 1.0
sigma = 0.05 * L
m = 1000
N = 100

def stevec(x, nicle):
    '''
    Funkcija kot parameter vzame x, array dimenzije (N, m), in nicle, array dimenzije (N, m).
    '''
    return v0 * np.exp(- (x - L) ** 2 / (2 * sigma ** 2)) * x * ss.j0(nicle * np.sqrt(1 - x / L))

def alt_stevec(x, nicle):
    '''
    f(x) = 1
    '''
    return x * ss.j0(nicle * x)

def imenovalec(x, nicle):
    '''
    Funkcija kot parametra sprejme x, array dimenzije m, in nicle, array dimenzije N, ki je N ničel Besslove funkcije J_0.
    Vrne array dimenzije (N, m).
    '''
    return (ss.j0(nicle * np.sqrt(1 - x / L))) ** 2 * x

def gauss(v0, sigma, x):
    return v0 * np.exp(- (x - L) ** 2 / (2 * sigma ** 2))

def stem_plot(n, koef):
    '''
    Makes a stem plot of |B_n| in relation to n
    '''
    x = np.arange(n, step=1)
    plt.stem(x, np.abs(koef), markerfmt='.')
    plt.xlabel(r'$n$')
    plt.ylabel(r'$\left| B_N \right|$')
    plt.show()
    pass

def cumulative_tail(n, koef):
    '''
    Makes a cumulative tail plot
    '''
    x = np.arange(n, step=1)
    plt.plot(x, np.cumsum(np.abs(koef)))
    plt.xlabel(r'$n$')
    plt.ylabel(r'$\sum_{n = 1}^{N}\left| B_N \right|^2$')
    plt.show()
    pass


# ničle Besselovih funkcij, ki jih shranimo v array dimenzij (N, m)
# (N,) preko broadcast_to v (m, N) in s .T v (N, m)

xi = np.broadcast_to(ss.jn_zeros(0, N), (m, N)).T

# frekvenca omega_n v arrayju dimenzije (N,)

omega = (1 / 2) * np.sqrt(g / L) * xi[:, 0]

# --------------- intermezzo ---------------

# za funkcijo f(x) vem, da imajo koeficienti vrednost
# B_n = 2/(xi_n * J_1 (xi_n))
# source: https://math.libretexts.org/Bookshelves/Differential_Equations/Introduction_to_Partial_Differential_Equations_(Herman)/05%3A_Non-sinusoidal_Harmonics_and_Special_Functions/5.05%3A_Fourier-Bessel_Series

# interval za integriranje

uniform = np.linspace(0, 1, m)

# izračun koeficientov preko integracije

interval = np.broadcast_to(uniform, (N, m))

# integracija

koef_int = np.trapezoid(alt_stevec(interval, xi), uniform, axis=1)

koef_int = (2 / ((ss.j1(xi[:, 0]) ** 2)) * koef_int) 
print('koef_int ', koef_int)

# izračun koeficientov preko identitete

koef_id = (2 / (xi[:, 0] * ss.j1(xi[:, 0])))
print('koef_id ', koef_id)

# prvih 50 členov
temp_int = 0
temp_id = 0
for i in range(N):
    temp_int += koef_int[i] * ss.j0(xi[i, 0] * uniform)
    temp_id += koef_id[i] * ss.j0(xi[i, 0] * uniform)

plt.plot(uniform, temp_int, label='int', alpha=0.5)
plt.plot(uniform, temp_id, label='id', alpha=0.5)
plt.legend()
plt.show()

plt.plot(uniform, np.abs(temp_int - temp_id))
plt.show()


# --------------- konec intermezza ---------------

# gradnja intervala velikosti (N, m)

non_uniform = np.r_[np.linspace(0, 0.5, 100, endpoint=False), np.linspace(0.5, 1, m)]
uniform = np.linspace(0, 1, m)

interval = np.broadcast_to(uniform, (N, m))

# gradnja matrike z rezultati, dimenzij (N, 3)

x = np.arange(3, step=1.0)
y = np.arange(N, step=1.0)

# matrika oblike [[0, 0, 0], [1, 1, 1], ..., [N, N, N]]
_, rez = np.meshgrid(x, y) 

# stevec
rez[:, 0] = np.trapezoid(alt_stevec(interval, xi), uniform, axis=1)
# imenovalec
rez[:, 1] = ss.j1(xi[:, 0]) ** 2
# koeficienti
rez[:, 2] = (2 / omega) * (rez[:, 0] / rez[:, 1]) 

# initial condition plot

x_kord = uniform

j0i = 0 # rez[N-1, 2] * ss.j0(xi[N-1, 0] * np.sqrt(1 - x_kord / L))

for i in range(N):
    j0i += rez[i, 2] * ss.j0(xi[i, 0] * uniform)

#plt.plot(x_kord, j0i)
#plt.plot(x_kord, gauss(v0, sigma, x_kord))

