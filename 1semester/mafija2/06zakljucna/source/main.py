#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import scipy.special as ss
import scipy.integrate as si
import cmasher as cmr


# Constants

g = 9.81
L = 1.0
v0 = 1.0
M = 1000
N = 100

def nicle(n, m):
    '''
    For easier integration, returns array of size (n, m) with the n-th zero of J0 along the axis =1 
    '''
    # ničle Besselovih funkcij, ki jih shranimo v array dimenzij (N, m)
    # (N,) preko broadcast_to v (m, N) in s .T v (N, m)
    xi = np.broadcast_to(ss.jn_zeros(0, n), (m, n)).T
    return xi
    

def stevec(x, n, m, zero=nicle):
    '''
    Funkcija kot parameter vzame x, array dimenzije (N, m), in nicle, array dimenzije (N, m).
    '''
    return ss.j0(zero(n,m) * np.sqrt(1 - x / L))

def alt_stevec(x, nicle):
    '''
    f(x) = 1
    '''
    return x * ss.j0(nicle * 20 * x)

def stem_plot(n, color, koef):
    '''
    Makes a stem plot of |B_n| in relation to n
    '''
    x = np.arange(n, step=1)
    plt.stem(x, np.abs(koef), markerfmt='.')
    plt.xlabel(r'$n$')
    plt.ylabel(r'$\left| B_N \right|$')
    plt.title(r'Vrednosti posameznih koeficientov do $N = 100$')
    plt.savefig('stemplot.png')
    plt.show()
    pass

def cumulative_tail(n, m, color, koef):
    '''
    Makes a cumulative tail plot
    '''
    x = np.arange(n, step=1)
    plt.plot(x, np.cumsum(np.abs(koef)))
    plt.xlabel(r'$n$')
    plt.ylabel(r'$\Sigma$')
    plt.title(r'Kumulativna vsota koeficientov do $N = 100$')
    plt.savefig('cumtail.png')
    plt.show()
    pass

def u_n(n, m, koef, zero=nicle):
    '''
    Returns the sum of n-terms in the Fourier Bessel series
    '''
    uniform = np.linspace(0.95 * L, 1 * L, m)
    # vsota
    sumJ = 0
    for i in range(n):
        sumJ += koef(n, m)[i] * ss.j0(zero(n, m)[i, 0] * np.sqrt(1 - uniform/L))
    return sumJ

def u_nt(n, m, cas, koef, zero=nicle):
    '''
    Returns the sum of n-terms in the Fourier Bessel series with time evolution
    '''
    uniform = np.linspace(0.95 * L, 1 * L, m)
    # frekvence
    omega = (1 / 2) * np.sqrt(g / L) * zero(n, m)[:, 0]
    # vsota
    sumJ = 0
    for i in range(n):
        sumJ += koef(n, m)[i] * ss.j0(zero(n, m)[i, 0] * np.sqrt(1 - uniform/L)) *\
            np.sin(omega * cas)
    return sumJ

def koef_n(n, m, zero=nicle, kick=0.95, L=1, g=9.81):
    '''
    Returns array of length n with coefficients of Fourier-Bessel series.
    '''
    # frekvenca omega_n v arrayju dimenzije (n,)
    omega = (1 / 2) * np.sqrt(g / L) * zero(n, m)[:, 0]

    # gradnja intervala velikosti (N, m)
    uniform = np.linspace(0.95 * L, 1 * L, m)
    interval = np.broadcast_to(uniform, (n, m))

    # gradnja matrike z rezultati, dimenzij (N, 3)
    x = np.arange(3, step=1.0)
    y = np.arange(n, step=1.0)

    # matrika oblike [[0, 0, 0], [1, 1, 1], ..., [N, N, N]]
    _, rez = np.meshgrid(x, y)
    # stevec
    rez[:, 0] = np.trapezoid(stevec(interval, n, m, zero), uniform, axis=1)
    # imenovalec
    rez[:, 1] = ss.j1(zero(n, m)[:, 0]) ** 2
    # koeficienti
    rez[:, 2] = (1 / omega) * (rez[:, 0] / rez[:, 1]) 
    return rez[:, 2]

def koef_n_omega(n, m, zero=nicle, kick=0.95, L=1, g=9.81):
    '''
    Returns array of length n with coefficients of Fourier-Bessel series.
    '''
    # frekvenca omega_n v arrayju dimenzije (n,)
    omega = (1 / 2) * np.sqrt(g / L) * zero(n, m)[:, 0]

    # gradnja intervala velikosti (N, m)
    uniform = np.linspace(0.95 * L, 1 * L, m)
    interval = np.broadcast_to(uniform, (n, m))

    # gradnja matrike z rezultati, dimenzij (N, 3)
    x = np.arange(3, step=1.0)
    y = np.arange(n, step=1.0)

    # matrika oblike [[0, 0, 0], [1, 1, 1], ..., [N, N, N]]
    _, rez = np.meshgrid(x, y)
    # stevec
    rez[:, 0] = np.trapezoid(stevec(interval, n, m, zero), uniform, axis=1)
    # imenovalec
    rez[:, 1] = ss.j1(zero(n, m)[:, 0]) ** 2
    # koeficienti
    rez[:, 2] =  (rez[:, 0] / rez[:, 1]) 
    return rez[:, 2]


c1, c2, c3 = cmr.take_cmap_colors("cmr.cosmic", 3, cmap_range=(0.2, 0.8),
                                  return_fmt="hex")


# Gibbs phenomenon

NG = 200

plt.plot(np.linspace(0.95 * L, 1 * L, M), u_n(NG, M, koef=koef_n_omega), color=c2)
plt.xlabel(r'$x$')
plt.ylabel(r'$\left. \partial_t u_N(x, t)\right|_{t = 0}$')
plt.title(r'Gibbsov fenomen za $N = 200$')
plt.savefig('GibbsPhenomenon.png')
plt.show()

# time evolution



'''
# stem plot for \omega and \omega_n

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 6))

x = np.arange(N, step=1)
ax1.stem(x, np.log(np.abs(koef_n(N, M))), linefmt=c2, basefmt=c3, markerfmt='.')
ax1.set_ylabel(r'$\left| B_n \right|$')

ax2.stem(x, np.log(np.abs(koef_n_omega(N, M))), linefmt=c2, basefmt=c3, markerfmt='.')
ax2.set_ylabel(r'$\left|\omega_n B_n \right|$')

fig.supxlabel(r'$n$')
fig.suptitle(r'Vrednosti koeficientov do $N = 100$')
plt.savefig('stemplot.png')

'''

'''
# cumulative tail plot

fig, (ax1, ax2) = plt.subplots(2, 1)

N1 = 100
x1 = np.arange(N1, step=1)
ax1.plot(x1, np.cumsum(np.abs(koef_n(N1, M))), color=c2)


N2 = 1000
x2 = np.arange(N2, step=1)
ax2.plot(x2, np.cumsum(np.abs(koef_n(N2, M))), color=c2)

fig.supxlabel(r'$n$')
fig.supylabel(r'$\sum_n {\left|B_n \right|}^2$')
fig.suptitle(r'Konvergenca kumulativne vsote koeficientov')
fig.tight_layout()
plt.savefig('cumtail.png')
'''

# convergence in L2 space

mja = [1000, 10000]
sums = [i for i in range(25, 275, 25)]

# colors

colors = cmr.take_cmap_colors("cmr.cosmic", len(sums), cmap_range=(0.2, 1),
                              return_fmt="hex")
'''
for u, ax in zip(mja, [ax1, ax2]):
    for s, c in zip(sums, colors):
        print("Racunam za ", s)
        uniform = np.linspace(0.95, 1, u)
        error = np.sqrt(np.trapezoid((1 - u_n(s, u, koef=koef_n_omega)) ** 2, uniform))
        ax.scatter(s, np.max((error)), marker='.', color=c)

ax1.set_title(r'M = 1000')
ax2.set_title(r'M = 10000')
fig.supxlabel(r'$n$')
fig.supylabel(r'$I_{L_2}$')
fig.suptitle(r'Konvergenca v $L_2$ normi')
fig.tight_layout
fig.savefig('konvergencaL2.png')
fig.show()

'''

'''
# pointwise convergence

for u, ax in zip(mja, [ax1, ax2]):
    for s, c in zip(sums, colors):
        print("Racunam za ", s)
        error = np.max(np.abs(1 - u_n(s, u, koef_n)[200:]))
        ax.scatter(s, error, marker='.', color=c)
    
ax1.set_title(r'M = 1000')
ax2.set_title(r'M = 10000')
fig.supxlabel(r'$n$')
fig.supylabel(r'$I_{\text{tocke}}$')
fig.suptitle(r'Konvergenca po točkah')
fig.tight_layout
fig.savefig('konvergencaTocke.png')
plt.show()

'''
