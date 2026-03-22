#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import cmasher as cmr

# colours

c1, c2, c3 = cmr.take_cmap_colors("cmr.bubblegum", 3, cmap_range=(0,1), return_fmt="hex")


def rhs(x, k = 1):
    return (x / (1 - k  * x ** 2))

axis_x = np.linspace(0, 4 * np.pi, num = 10000)
print("axis_x ", axis_x)

value1 = rhs(axis_x)
value2 = np.tan(axis_x)

plt.plot(axis_x, value1, c=c1, label="RHS")
plt.plot(axis_x, value2, c=c3, label="LHS")

# setting ticks

tick_values = [0, 1 * np.pi, 2 * np.pi, 3 * np.pi, 4 * np.pi, ]
tick_labels = ["0", r"$\pi$", r"$2\pi$", r"$3\pi$", r"$4\pi$", ]

plt.xticks(tick_values, tick_labels)

plt.ylim(-10, 10)
plt.grid()
plt.hlines(y = 0, xmin = 0, xmax=4 * np.pi, color='k', zorder=10, alpha = 0.5)

plt.title(r"Grafične rešitve $\tan x = \frac{x}{1 - kx ^2}$")
plt.legend()
plt.savefig("12mehk_eigenfreq.png")
