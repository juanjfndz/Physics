"""
======================
Triangular 3D surfaces
======================

Plot a 3D surface with a triangular mesh.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


n_radii = 8
n_angles = 36

# Make radii and angles spaces (radius r=0 omitted to eliminate duplication).
radii = np.linspace(0.125, 1.0, n_radii)
angles = np.linspace(0, 2 * np.pi, n_angles, endpoint=False)[..., np.newaxis]

# Convert polar (radii, angles) coords to cartesian (x, y) coords.
# (0, 0) is manually added at this stage, so there will be no duplicate
# points in the (x, y) plane.
x = np.append(0, (radii * np.cos(angles)).flatten())
y = np.append(0, (radii * np.sin(angles)).flatten())

# Compute z to make the pringle surface.
z = x ** 3 - 3 * x * y ** 2

fig = plt.figure()
ax = fig.gca(projection="3d")

# Hide grid lines.
ax.grid(False)
ax.axis("off")
ax.view_init(20, 165)

ax.set_zlim(-2, 1)

ax.plot_trisurf(
    x, y, z, linewidth=0.2, antialiased=True, cmap="viridis"
)

plt.show()
plt.savefig("Fondo_poster_2.png", dpi=300)

