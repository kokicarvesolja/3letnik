#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import scipy.special as ss
import scipy.integrate as si

# Constants

g = 9.81
L = 1.0
v0 = 1.0
m = 1000
N = 100

def stevec(x, nicle):
    '''
    Funkcija kot parameter vzame x, array dimenzije (N, m), in nicle, array dimenzije (N, m).
    '''
    return ss.j0(nicle * np.sqrt(1 - x / L))

def alt_stevec(x, nicle):
    '''
    f(x) = 1
    '''
    return x * ss.j0(nicle * 20 * x)

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

def u_n(n, koef):
    '''
    Returns the sum of n-terms in the Fourier Bessel series
    '''
    temp = 0
    for i in range(n):
        temp += koef[i] * ss.j0(xi[i, 0] * np.sqrt(1 - uniform/L))
    return temp
        

# ničle Besselovih funkcij, ki jih shranimo v array dimenzij (N, m)
# (N,) preko broadcast_to v (m, N) in s .T v (N, m)

xi = np.broadcast_to(ss.jn_zeros(0, N), (m, N)).T

# frekvenca omega_n v arrayju dimenzije (N,)

omega = (1 / 2) * np.sqrt(g / L) * xi[:, 0]

# gradnja intervala velikosti (N, m)

uniform = np.linspace(0.95, 1, m)

interval = np.broadcast_to(uniform, (N, m))

# gradnja matrike z rezultati, dimenzij (N, 3)

x = np.arange(3, step=1.0)
y = np.arange(N, step=1.0)

# matrika oblike [[0, 0, 0], [1, 1, 1], ..., [N, N, N]]
_, rez = np.meshgrid(x, y) 

# stevec
rez[:, 0] = np.trapezoid(stevec(interval, xi), uniform, axis=1)
# imenovalec
rez[:, 1] = ss.j1(xi[:, 0]) ** 2
# koeficienti
rez[:, 2] = (rez[:, 0] / rez[:, 1]) 

stem_plot(N, rez[:, 2])

# initial condition plot

j0i = 0 # rez[N-1, 2] * ss.j0(xi[N-1, 0] * np.sqrt(1 - x_kord / L))

for i in range(N):
    j0i += rez[i, 2] * ss.j0(xi[i, 0] * np.sqrt(1 - uniform/L))

x_kord = np.r_[np.linspace(0, 0.95, m), uniform]
j0i = np.r_[np.zeros(m), j0i]

plt.plot(x_kord, j0i)
plt.ylim(0, 2)
plt.show()
#plt.plot(x_kord, gauss(v0, sigma, x_kord))

