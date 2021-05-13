from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from itertools import product, combinations


fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_aspect("equal")


for i in range(0, 1000):
    inclination = np.arccos(np.random.uniform(-1, 1))
    azimuth = 2.0 * np.pi * np.random.uniform(0, 1)

    ax.scatter(np.sin(inclination) * np.cos(azimuth), np.sin(inclination)
               * np.sin(azimuth), np.cos(inclination), color="b", s=5)
plt.show()
