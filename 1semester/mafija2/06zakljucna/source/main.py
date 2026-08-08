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
m = 10000
N = 100

def stevec(x, nicle):
    '''
    Funkcija kot parameter vzame x, array dimenzije (N, m), in nicle, array dimenzije (N, m).
    '''
    return v0 * np.exp(- (x - L) ** 2 / (2 * sigma ** 2)) * ss.j0(nicle * np.sqrt(1 - x / L))

def imenovalec(x, nicle):
    '''
    Funkcija kot parametra sprejme x, array dimenzije m, in nicle, array dimenzije N, ki je N ničel Besslove funkcije J_0.
    Vrne array dimenzije (N, m).
    '''
    return (ss.j0(nicle * np.sqrt(1 - x / L))) ** 2

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

xi = np.broadcast_to(ss.jn_zeros(0, N), (m + 100, N)).T

# frekvenca omega_n v arrayju dimenzije (N,)

omega = (1 / 2) * np.sqrt(g / L) * xi[:, 0]

# gradnja intervala velikosti (N, m)

non_uniform = np.r_[np.linspace(0, 0.5, 100, endpoint=False), np.linspace(0.5, 1, m)]

interval = np.broadcast_to(non_uniform, (N, m+100))

# gradnja matrike z rezultati, dimenzij (N, 3)

x = np.arange(3, step=1.0)
y = np.arange(N, step=1.0)

# matrika oblike [[0, 0, 0], [1, 1, 1], ..., [N, N, N]]
_, rez = np.meshgrid(x, y) 

# stevec
rez[:, 0] = np.trapezoid(stevec(interval, xi), non_uniform, axis=1)
# imenovalec
rez[:, 1] = np.trapezoid(imenovalec(interval, xi), non_uniform, axis=1)
# koeficienti
rez[:, 2] = (1 / omega) * (rez[:, 0] / rez[:, 1]) 

# initial condition plot

x_kord = np.r_[np.linspace(0, 0.5, 100, endpoint=False), np.linspace(0.5, 1, m)]

j0i = 0 # rez[N-1, 2] * ss.j0(xi[N-1, 0] * np.sqrt(1 - x_kord / L))

for i in range(N):
    j0i += rez[i, 2] * ss.j0(xi[i, 0] * np.sqrt(1 - x_kord / L))

plt.plot(x_kord, j0i)
#plt.plot(x_kord, gauss(v0, sigma, x_kord))

plt.show()



'''

for i in range(N):
    num = np.trapezoid(stevec(interval, i + 1), interval)
    den = np.trapezoid(imenovalec(interval, i + 1), interval)
    b_n = (1 / omega(i + 1)) * (num / den)
    rez[i] = b_n 

#plt.plot(np.arange(N, step=1), np.cumsum(rez))

plt.stem(np.arange(N, step=1), np.abs(rez), markerfmt='.')

plt.show()

plt.plot(np.arange(N, step=1), np.cumsum(np.abs(rez)))

plt.show()
'''
    

    
    




