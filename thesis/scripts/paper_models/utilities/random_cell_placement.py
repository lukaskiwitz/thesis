import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from scipy.spatial import KDTree, distance_matrix

def random_3D_unit_vector():
    """
    Generates a random 3D unit vector (direction) with a uniform spherical distribution
    Algo from http://stackoverflow.com/questions/5408276/python-uniform-spherical-distribution
    :return:
    """
    phi = np.random.uniform(0,np.pi*2)
    costheta = np.random.uniform(-1,1)

    theta = np.arccos( costheta )
    x = np.sin( theta) * np.cos( phi )
    y = np.sin( theta) * np.sin( phi )
    z = np.cos( theta )
    return np.array([x,y,z])

rho = 5
penalty = 1

x, y, z = np.random.uniform(0, 140, (3,1000))
# xv, yv, zv = np.meshgrid(x, y, z)

positions = np.vstack(list(zip(x.ravel(), y.ravel(), z.ravel())))

# randomizer = lambda r: r + 10 * random_3D_unit_vector()
# r_positions = np.array([randomizer(pos) for pos in positions])

dist_m = distance_matrix(positions, positions)
dist_m[np.tril_indices(len(dist_m[0]), k = 0)] = None

iterations = 1000

for i in range(iterations):
    too_close = np.argwhere(dist_m < rho * 2)
    for tc in too_close:
        cell_1 = np.round(positions[tc[0]], 2)
        cell_2 = np.round(positions[tc[1]], 2)
        c_c_v = ((cell_2 - cell_1) / 2) * 0.1
        positions[tc[0]] += - c_c_v
        positions[tc[1]] += + c_c_v

    dist_m = distance_matrix(positions, positions)
    dist_m[np.tril_indices(len(dist_m[0]), k = 0)] = None

    too_close = np.argwhere(dist_m < rho * 2)
    if len(too_close) == 0:
        print(i)
        break

dist_list = []
for coord in range(len(dist_m)):
    # dist_list.append(np.nanmin(dist_m[coord]))
    dist_list.append(np.sort(dist_m[coord])[6])

flat_dist_array = np.array(dist_m).flatten()
plt.hist(dist_list, bins = 46)
plt.axvline(np.nanmean(dist_list), 0, 130, color="orange")
plt.axvline(20, 0, 130, color="orange")
plt.xlim((0, None))
plt.show()
#
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(positions[:,0],positions[:,1],positions[:,2], marker=".", s=20)
# ax.scatter(x, y, z, marker="o")
plt.show()

