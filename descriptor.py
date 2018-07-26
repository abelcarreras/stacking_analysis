import numpy as np
from scipy.optimize import leastsq


def get_plane(coords, direction=None):
    # Generate initial guess
    p0 = np.cross(coords[0] - coords[2], coords[0] - coords[-1]) / np.linalg.norm(np.cross(coords[0] - coords[2],
                                                                                           coords[0] - coords[-1]))

    # Fitting function to a plane
    def fitfunc(p, coords):
        average = np.average(coords, axis=0)
        return np.array([np.dot(p, average - c) for c in coords])

    # Error function (including force norm(normal) = 1)
    errfunc = lambda p, x: fitfunc(p, x)**2 + (np.linalg.norm(p) - 1.0)**2

    p1, flag = leastsq(errfunc, p0, args=(coords,))


    # Check final result
    point = np.average(coords, axis=0)
    normal = np.array(p1)/np.linalg.norm(p1)

    if direction is not None:
        vector = coords[direction[1]] - coords[direction[0]]
        # proj = vector - np.dot(vector, normal)*normal
        projected = np.cross(normal,  np.cross(vector, normal))
        projected /= np.linalg.norm(projected)

        return point, normal, projected

    return point, normal


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D


    # generate test data
    point0 = np.array([0, 3, 4])
    normal0 = np.array([1, -3, 4.2])

    d = -point0.dot(normal0)
    ix = int(point0[0]) - 10
    iy = int(point0[1]) - 10
    xx, yy = np.meshgrid(range(ix, ix+20, 2), range(iy, iy +20, 2))
    zz = (-normal0[0] * xx - normal0[1] * yy - d) * 1. / normal0[2] + np.random.random(len(xx[0])) - 0.5

    coords = np.array([np.ndarray.flatten(xx), np.ndarray.flatten(yy), np.ndarray.flatten(zz)]).T

    plt.plot(coords)
    plt.show()

    vec = [23, 60]
    #point, normal = get_plane(coords)
    point, normal, vector = get_plane(coords, direction=vec)

    print('center', point)
    print('normal', normal)
    print('vector', vector)
    v2 = np.cross(np.cross(normal, np.array(coords[vec[1]] - coords[vec[0]])),normal)


    # Plot data
    d = -point.dot(normal)
    ix = int(point[0])-10
    iy = int(point[1])-10
    xx, yy = np.meshgrid(range(ix, ix+20), range(iy, iy +20))
    zz = (-normal[0] * xx - normal[1] * yy - d) * 1. /normal[2]


    plt3d = plt.figure().gca(projection='3d')
    plt3d.plot_surface(xx, yy, zz, alpha=0.2)
    plt3d.scatter(*coords.T)
    plt3d.scatter(*point)
    plt3d.scatter(*coords[vec[0]], color='red')
    plt3d.scatter(*coords[vec[1]], color='red')
    plt3d.quiver(*point.tolist() + (vector*10).tolist())
    #plt3d.quiver(*point.tolist() + v2.tolist(), color='orange')
    plt3d.quiver(*point.tolist() + (normal*10).tolist(), color='green')
    #plt3d.quiver(*point.tolist() + np.array(coords[vec[1]] - coords[vec[0]]).tolist(), color='red')

    plt.show()
