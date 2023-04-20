#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Realización de una variedad es un espacio tridimensional
"""

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

def sphere(r):
    u = np.linspace(0, np.pi/8, 25)
    v = np.linspace(0, np.pi/8, 25)
    x = r * np.outer(np.cos(u), np.sin(v))
    y = r * np.outer(np.sin(u), np.sin(v))
    z = r * np.outer(np.ones(np.size(u)), np.cos(v))
    return x, y, z

x, y, z = sphere(10)

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
ax.set_zlim(9.3, 10.2)
ax.set_ylim(-0.2, 0.8)

# Desplazarlo
ax.plot_wireframe(x, y, z, rstride=8, cstride=8, color='b')

ax.quiver(x[8][8], y[8][8], z[8][8], 1, 0.5, 0.5, length=0.25,
          normalize=True, color="r", arrow_length_ratio=0.2)
ax.quiver(x[8][8], y[8][8], z[8][8], 1, -0.1, 0.5, length=0.25,
          normalize=True, color="r", arrow_length_ratio=0.2)

ax.plot([x[8][8]+0.11, x[8][8]+0.12],
        [y[8][6]+0.03, y[8][12]-0.03],
        z[8][8]+0.05,
        color='k')

ax.quiver(x[16][16], y[16][16], z[16][16], 1, 0.2, 0.5, length=0.25,
          normalize=True, color="g", arrow_length_ratio=0.2)
