import numpy as np
import matplotlib.pyplot as plt


def plot_simulations(sim_fname, imsize, mas_in_pix, delta=100, jet_delta=100,
                     contr_delta=10, each=1, rstride=6, cstride=6):
    """
    Plot simulation results in 3D projection.

    :param sim_fname:
        Path to file with simulations result of Stokes I.
    :param delta: (optional)
        Cut image from sides. Image will be shorter by ``2*delta`` pixels from
        sides. (default: ``100``)
    :param contr_delta: (optional)
        Space to leave on contr-jet side. (default: ``10``)

    :return:
        Instance of ``Figure``.
    """
    # This import is needed for 3D projections
    from mpl_toolkits.mplot3d import Axes3D
    # Making transparent color map
    from matplotlib import cm
    # theCM = cm.Blues
    # theCM._init()
    # alphas = np.abs(np.linspace(-1.0, 1.0, theCM.N))
    # theCM._lut[:-3, -1] = alphas

    image = np.loadtxt(sim_fname)
    print(sum(image))
    y = np.arange(-imsize/2+delta, imsize/2-delta, each, dtype=float)
    x = np.arange(-contr_delta, imsize/2-jet_delta, each, dtype=float)
    x *= mas_in_pix*each
    y *= mas_in_pix*each
    xx, yy = np.meshgrid(x, y)

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.grid(False)
    plt.hold(True)

    # mask = image < 0.0001
    # image = np.ma.array(image, mask=mask)
    surf = ax.plot_surface(xx, yy,
                           image[delta:-delta:each, imsize/2-contr_delta:imsize-jet_delta:each],
                           cmap='OrRd', linewidth=2, rstride=rstride,
                           cstride=cstride,
                           antialiased=True)
    # fig.colorbar(surf, shrink=0.5, aspect=5)
    ax.set_xlabel("[mas]", fontsize=14)
    ax.set_ylabel("[mas]", fontsize=14)
    ax.set_zlabel("Flux [Jy/pix]", fontsize=14)

    # ax.set_zlim([0.0002, None])
    fig.tight_layout()
    plt.show()
    return fig


simulation_file = '/home/ilya/github/bck/jetshow/uvf_mf_adds/map_i_09_15.4.txt'
pixel_size_mas = 0.002533
imsize = 756
from matplotlib import rcParams
rcParams[u'patch.force_edgecolor'] = True
rcParams[u'figure.edgecolor'] = 'black'
fig = plot_simulations(simulation_file, imsize, pixel_size_mas, each=1,
                       delta=300, jet_delta=250, rstride=6, cstride=6)
# fig.savefig('/home/ilya/github/bck/jetshow/uvf/bk_3d.png',
#             bbox_inches='tight', dpi=300)