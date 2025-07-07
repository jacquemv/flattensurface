import flattensurface as flat

vertices = [
    [0, 0, 0],
    [1, 0, 0],
    [1, 1, 0],
    [0, 1, 0],
    [0, 0, 1],
    [1, 0, 1],
    [1, 1, 1],
    [0, 1, 1],
]
triangles = [
    [0, 2, 1],
    [0, 3, 2],
    [0, 7, 3],
    [0, 4, 7],
    [4, 6, 7],
    [4, 5, 6],
    [1, 6, 5],
    [1, 2, 6],
    [3, 6, 2],
    [3, 7, 6],
]

uv = flat.global_parameterization(
    vertices, triangles, method='meanvalue'
)


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize=(8, 3.5))
ax = fig.add_subplot(121, projection='3d')
vertices = np.array(vertices)
ax.plot_trisurf(vertices[:, 0], vertices[:, 1], vertices[:, 2], 
                triangles=triangles, 
                color="0.85", edgecolor='C0', linewidth=1.5)
ax.view_init(-30, -60)
plt.axis('off')
plt.title('Box with boundary')

ax = fig.add_subplot(122)
plt.triplot(uv[:, 0], uv[:, 1], triangles=triangles)
plt.axis('equal')
plt.axis('off')
plt.title('Planar parametrization')
plt.tight_layout()
#plt.savefig('example.png')
plt.show()