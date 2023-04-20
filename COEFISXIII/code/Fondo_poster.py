import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Hide grid lines
ax.grid(False)
ax.axis('off')

def f(x, y):
    r_s = 0.1
    r = np.sqrt(x ** 2 + y ** 2)
    return 2 * r_s * np.sqrt(r / r_s - 1)

x = np.linspace(-160, 160, 30)
y = np.linspace(-160, 160, 30)
z = f(x, y)
X, Y = np.meshgrid(x, y)
Z = f(X, Y)

ax.contour3D(X, Y, Z, 50, cmap='cividis')

plt.show()

