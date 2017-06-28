import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

# TODO: Need functions that could read config.json and create corresponding
# plots of the physical parameters.
def scatter_3d(x1, x2, y, xlabel='u', ylabel='v', zlabel='flux', xlim3d=None,
               ylim3d=None, zlim3d=None):
    """
    Do 3d scatter plot.
    :param x1:
        Array-like of first coordinate.
    :param x2:
        Array-like of second coordinate.
    :param y:
        Array-like of values.
    :param xlabel:
    :param ylabel:
    :param zlabel:
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x1, x2, y, c='r', marker='o')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    # Configure axes
    if xlim3d is not None:
        ax.set_xlim3d(xlim3d[0], xlim3d[1])
    if ylim3d is not None:
        ax.set_ylim3d(ylim3d[0], ylim3d[1])
    if zlim3d is not None:
        ax.set_zlim3d(zlim3d[0], zlim3d[1])
    plt.show()


def beta(x, y, z, gamma0=10.):
    """
    Velocity field of jet.
    :param x, y, z:
        Rectangular coordinates. (3, N,) N-number of points
    :param gamma0:
        Lorentz-factor of jet at ``z_0`` & ``theta_0``
    :return:
        (3, N,) array of vectors of velocity of jet at given N point ``xyz`` in
        rectangular coordinates.
    """
    # r-component of velocity in spherical coordinates
    value = np.sqrt(1. - 1. / gamma0 ** 2.)
    result =  np.array([value * x / np.sqrt(x * x + y * y + z * z),
                        value * y / np.sqrt(x * x + y * y + z * z),
                        value * z / np.sqrt(x * x + y * y + z * z)])
    return result


if __name__ == '__main__':
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    x, y, z = np.meshgrid(np.arange(-5., 5., 1.0),
                          np.arange(-5., 5., 1.0),
                          np.arange(0.1, 10., 1.0))
    from
    result = beta(x, y, z)

    ax.quiver(x, y, z, result[0], result[1], result[2], length=1.0)

    plt.show()