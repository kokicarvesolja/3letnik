#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import scipy.special as ss
import scipy.integrate as si

# Constants

g = 9.81
L = 1.0
v0 = 1.0
sigma = 0.5 * L
m = 1000
N = 1000

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

# ničle Besselovih funkcij, ki jih shranimo v array dimenzij (N, m)
# (N,) preko broadcast_to v (m, N) in s .T v (N, m)

xi = np.broadcast_to(ss.jn_zeros(0, N), (m, N)).T

# frekvenca omega_n v arrayju dimenzije (N,)

omega = (1 / 2) * np.sqrt(g / L) * xi[:, 0]

# preverimo ortogonalnost

n = [2, 5, 6, 7]
x_kord = np.linspace(0, 1, m)

# druga ničla

xi_2 = ss.jn_zeros(0, 2)[-1]

# tretja ničla

xi_3 = ss.jn_zeros(0, 3)[-1]

value2 = si.quad(lambda x: ss.j0(xi_2 * np.sqrt(1 - x / L)) * ss.j0(xi_2 * np.sqrt(1 - x / L)), 0, 1)
value3 = si.quad(lambda x: ss.j0(xi_2 * np.sqrt(1 - x / L)) * ss.j0(xi_3 * np.sqrt(1 - x / L)), 0, 1)

print("scipy.integrate ", value2, ' ', value3)

valueNP2 = np.trapezoid(ss.j0(xi_2 * np.sqrt(1 - x_kord/ L)) * ss.j0(xi_2 * np.sqrt(1 - x_kord / L)), x_kord)
valueNP3 = np.trapezoid(ss.j0(xi_2 * np.sqrt(1 - x_kord/ L)) * ss.j0(xi_3 * np.sqrt(1 - x_kord / L)), x_kord)

print("numpy ", valueNP2, ' ', valueNP3)

'''
# gradnja intervala velikosti (N, m)

interval = np.broadcast_to(np.linspace(0, 1, m), (N, m))

# gradnja matrike z rezultati, dimenzij (N, 3)

x = np.arange(3, step=1.0)
y = np.arange(N, step=1.0)

# matrika oblike [[0, 0, 0], [1, 1, 1], ..., [N, N, N]]
_, rez = np.meshgrid(x, y) 


rez[:, 0] = np.trapezoid(stevec(interval, xi), np.linspace(0, 1, m), axis=1)
rez[:, 1] = np.trapezoid(imenovalec(interval, xi), np.linspace(0, 1, m), axis=1)
rez[:, 2] = (1 / omega) * (rez[:, 0] / rez[:, 1]) 

# graf exp - sum, plot maksimalnega errorja razlike

# vmesni rezultat

temp = 0
for i in range(N):
    # vsota Besselovih funkcija
    temp += rez[i, 2] * ss.j0(xi[i] * np.sqrt(1 - x_kord / L)) 

    # razlika od Gaussovke
    diff = np.abs(gauss(v0, sigma, x_kord)[-1] - temp[-1])

    # maksimalna error, plot
    plt.scatter(i, np.max(diff), marker='.')

    
#plt.stem(np.arange(N, step=1), np.abs(omega * rez[:, 2]), markerfmt='.')
#plt.semilogy(np.arange(N, step=1), np.cumsum(np.abs(rez[:, 2])))

plt.xlabel('N')
plt.ylabel(r'$\left| \exp{- \frac{(x - L)^2}{2 \sigma^2}} - \sum_{n = 1}^{N} B_n J_n (\ldots) \right|$')
plt.title('Razlika vrednosti samo pri x = L = 1')
plt.show()

plt.plot(x_kord, temp)
plt.plot(x_kord, gauss(v0, sigma, x_kord), color='r')
plt.show()
'''



'''
# gradnja matrike Y dimenzij (N, M)

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
    

    
    




