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



# FIXME: Don't use ``Components`` as arguments here
def plot(contours=None, colors=None, vectors=None, vectors_values=None, x=None,
         y=None, blc=None, trc=None, cmap='hsv', abs_levels=None,
         rel_levels=None, min_abs_level=None, min_rel_level=None, k=2, vinc=2,
         show_beam=False, beam_corner='ll', beam=None, contours_mask=None,
         colors_mask=None, vectors_mask=None, plot_title=None, color_clim=None,
         outfile=None, outdir=None, ext='png', close=False, slice_points=None,
         beam_place='ll', colorbar_label=None, show=True, contour_color='k',
         beam_edge_color='black', beam_face_color='green', beam_alpha=0.3,
         show_points=None, components=None, slice_color='black',
         plot_colorbar=True):
    """
    Plot image(s). Need much work to clean it!

    :param contours: (optional)
        Numpy 2D array (possibly masked) that should be plotted using contours.
    :param colors: (optional)
        Numpy 2D array (possibly masked) that should be plotted using colors.
    :param vectors: (optional)
        Numpy 2D array (possibly masked) that should be plotted using vectors.
    :param vectors_values: (optional)
        Numpy 2D array (possibly masked) that should be used as vector's lengths
        when plotting ``vectors`` array.
    :param x: (optional)
        Iterable of x-coordinates. It's length must be comparable to that part
        of image to display. If ``None`` then don't plot coordinates - just
        pixel numbers. (default=``None``)
    :param y: (optional)
        Iterable of y-coordinates. It's length must be comparable to that part
        of image to display. If ``None`` then don't plot coordinates - just
        pixel numbers. (default=``None``)
    :param blc: (optional)
        Iterable of two values for Bottom Left Corner (in pixels). Must be in
        range ``[1, image_size]``. If ``None`` then use ``(1, 1)``. (default:
        ``None``)
    :param trc: (optional)
        Iterable of two values for Top Right Corner (in pixels). Must be in
        range ``[1, image_size]``. If ``None`` then use ``(image_size,
        image_size)``. (default: ``None``)
    :param cmap: (optional)
        Colormap to use for plotting colors. Available color maps could be
        printed using ``sorted(m for m in plt.cm.datad if not
        m.endswith("_r"))`` where ``plt`` is imported ``matplotlib.pyplot``.
        For further details on plotting available colormaps see
        http://matplotlib.org/1.2.1/examples/pylab_examples/show_colormaps.html.
        (default: ``hsv``)
    :param abs_levels: (optional)
        Iterable of absolute levels. If ``None`` then construct levels in other
        way. (default: ``None``)
    :param min_abs_level: (optional)
        Values of minimal absolute level. Used with conjunction of ``factor``
        argument for building sequence of absolute levels. If ``None`` then
        construct levels in other way. (default: ``None``)
    :param rel_levels: (optional)
        Iterable of relative levels. If ``None`` then construct levels in other
        way. (default: ``None``)
    :param min_rel_level: (optional)
        Values of minimal relative level. Used with conjunction of ``factor``
        argument for building sequence of relative levels. If ``None`` then
        construct levels in other way. (default: ``None``)
    :param k: (optional)
        Factor of incrementation for levels. (default: ``2.0``)
    :param show_beam: (optional)
        Convertable to boolean. Should we plot beam in corner? (default:
        ``False``)
    :param beam_corner: (optional)
        Place (corner) where to plot beam on map. One of ('ll', 'lr', 'ul',
        'ur') where first letter means lower/upper and second - left/right.
        (default: ``ll'')
    :param beam: (optional)
        If ``show_beam`` is True then ``beam`` should be iterable of major axis,
        minor axis [mas] and beam positional angle [deg]. If no coordinates are
        supplied then beam parameters must be in pixels.
    :param colorbar_label: (optional)
        String to label colorbar. If ``None`` then don't label. (default:
        ``None``)
    :param slice_points: (optional)
        Iterable of 2 coordinates (``y``, ``x``) [mas] to plot slice. If
        ``None`` then don't plot slice. (default: ``None``)
    :param show_points: (optional)
        Iterable of 2 coordinates (``y``, ``x``) [mas] to plot points. If
        ``None`` then don't plot points. (default: ``None``)
    :param plot_colorbar: (optional)
        If colors is set then should we plot colorbar? (default: ``True``).

    :note:
        ``blc`` & ``trc`` are AIPS-like (from 1 to ``imsize``). Internally
        converted to python-like zero-indexing. If none are specified then use
        default values. All images plotted must have the same shape.
    """
    label_size = 14
    matplotlib.rcParams['xtick.labelsize'] = label_size
    matplotlib.rcParams['ytick.labelsize'] = label_size
    matplotlib.rcParams['axes.titlesize'] = label_size
    matplotlib.rcParams['axes.labelsize'] = label_size
    matplotlib.rcParams['font.size'] = label_size
    matplotlib.rcParams['legend.fontsize'] = label_size
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42

    image = None
    if contours is not None:
        image = contours
    elif colors is not None and image is None:
        image = colors
    elif vectors is not None and image is None:
        image = vectors

    if image is None:
        raise Exception("No images to plot!")
    if x is None:
        x = np.arange(image.shape[0])
        factor_x = 1
    else:
        factor_x = 1. / mas_to_rad
    if y is None:
        y = np.arange(image.shape[1])
        factor_y = 1
    else:
        factor_y = 1. / mas_to_rad

    # Set BLC & TRC
    blc = blc or (1, 1,)
    trc = trc or image.shape
    # Use ``-1`` because user expect AIPS-like behavior of ``blc`` & ``trc``
    x_slice = slice(blc[1] - 1, trc[1], None)
    y_slice = slice(blc[0] - 1, trc[0],  None)

    # Create coordinates
    imsize_x = x_slice.stop - x_slice.start
    imsize_y = y_slice.stop - y_slice.start
    # In mas (if ``x`` & ``y`` were supplied in rad) or in pixels (if no ``x`` &
    # ``y`` were supplied)
    x_ = x[x_slice] * factor_x
    y_ = y[y_slice] * factor_y
    # With this coordinates are plotted as in Zhenya's map
    # x_ *= -1.
    # y_ *= -1.
    # Coordinates for plotting
    x = np.linspace(x_[0], x_[-1], imsize_x)
    y = np.linspace(y_[0], y_[-1], imsize_y)

    # Optionally mask arrays
    if contours is not None and contours_mask is not None:
        contours = np.ma.array(contours, mask=contours_mask)
    if colors is not None and colors_mask is not None:
        colors = np.ma.array(colors, mask=colors_mask)
    if vectors is not None and vectors_mask is not None:
        vectors = np.ma.array(vectors, mask=vectors_mask)

    # Actually plotting
    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    ax.set_xlabel(u'Relative R.A. (mas)')
    ax.set_ylabel(u'Relative Decl. (mas)')

    # Plot contours
    if contours is not None:
        # If no absolute levels are supplied then construct them
        if abs_levels is None:
            print "constructing absolute levels for contours..."
            max_level = contours[x_slice, y_slice].max()
            # from given relative levels
            if rel_levels is not None:
                print "from relative levels..."
                # Build levels (``pyplot.contour`` takes only absolute values)
                abs_levels = [-max_level] + [max_level * i for i in rel_levels]
                # If given only min_abs_level & increment factor ``k``
            else:
                # from given minimal absolute level
                if min_abs_level is not None:
                    print "from minimal absolute level..."
                    n_max = int(math.ceil(math.log(max_level / min_abs_level, k)))
                # from given minimal relative level
                elif min_rel_level is not None:
                    print "from minimal relative level..."
                    min_abs_level = min_rel_level * max_level / 100.
                    n_max = int(math.ceil(math.log(max_level / min_abs_level, k)))
                abs_levels = [-min_abs_level] + [min_abs_level * k ** i for i in
                                                 range(n_max)]
            print "Constructed absolute levels are: ", abs_levels
        # return y, x, contours[x_slice, y_slice]
        co = ax.contour(y, x, contours[x_slice, y_slice], abs_levels,
                        colors=contour_color)
        print "OK"
    if colors is not None:
        im = ax.imshow(colors[x_slice, y_slice], interpolation='none',
                       origin='lower', extent=[y[0], y[-1], x[0], x[-1]],
                       cmap=plt.get_cmap(cmap), clim=color_clim)
    if vectors is not None:
        if vectors_values is not None:
            # TODO: Does "-" sign because of RA increases to the left actually?
            # VLBIers do count angles from North to negative RA.
            u = -vectors_values[x_slice, y_slice] * np.sin(vectors[x_slice,
                                                                   y_slice])
            v = vectors_values[x_slice, y_slice] * np.cos(vectors[x_slice,
                                                                  y_slice])
        else:
            u = -np.sin(vectors[x_slice, y_slice])
            v = np.cos(vectors[x_slice, y_slice])

        if vectors_mask is not None:
            u = np.ma.array(u, mask=vectors_mask[x_slice, y_slice])
            v = np.ma.array(v, mask=vectors_mask[x_slice, y_slice])
        vec = ax.quiver(y[::vinc], x[::vinc], u[::vinc, ::vinc],
                        v[::vinc, ::vinc], angles='uv',
                        units='xy', headwidth=0., headlength=0., scale=0.005,
                        width=0.05, headaxislength=0.)
    # Set equal aspect
    ax.set_aspect('equal')

    if slice_points is not None:
        for single_slice in slice_points:
            ax.plot([single_slice[0][0], single_slice[1][0]],
                    [single_slice[0][1], single_slice[1][1]], color=slice_color)

    if show_points is not None:
        for point in show_points:
            ax.plot(point[0], point[1], '.k')

    if plot_title:
        title = ax.set_title(plot_title, fontsize='large')
    # Add colorbar if plotting colors
    if colors is not None:
        if plot_colorbar:
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="10%", pad=0.00)
            cb = fig.colorbar(im, cax=cax)
            if colorbar_label is not None:
                cb.set_label(colorbar_label)

    if show_beam:
        from matplotlib.patches import Ellipse
        e_height = beam[0]
        e_width = beam[1]
        r_min = e_height / 2
        if beam_place == 'lr':
            y_c = y[0] + r_min
            x_c = x[-1] - r_min
        elif beam_place == 'll':
            if y[0] > 0:
                y_c = y[0] - r_min
            else:
                y_c = y[0] + r_min
            if x[0] > 0:
                x_c = x[0] - r_min
            else:
                x_c = x[0] + r_min
        elif beam_place == 'ul':
            y_c = y[-1] - r_min
            x_c = x[0] + r_min
        elif beam_place == 'ur':
            y_c = y[-1] - r_min
            x_c = x[-1] - r_min
        else:
            raise Exception

        # FIXME: check how ``bpa`` should be plotted
        e = Ellipse((y_c, x_c), e_width, e_height, angle=beam[2],
                    edgecolor=beam_edge_color, facecolor=beam_face_color,
                    alpha=beam_alpha)
        ax.add_patch(e)

    if components:
        for comp in components:
            y_c = -comp.p[1]
            x_c = comp.p[2]
            if len(comp) == 6:
                e_height = comp.p[3]
                e_width = comp.p[3] * comp.p[4]
                e = Ellipse((y_c, x_c), e_width, e_height,
                            angle=90+180*comp.p[5]/np.pi,
                            edgecolor=beam_edge_color, facecolor='red',
                            alpha=beam_alpha)
            elif len(comp) == 4:
                c_size = comp.p[3]
                e = Circle((y_c, x_c), c_size,
                            edgecolor=beam_edge_color, facecolor='red',
                            alpha=beam_alpha)
            else:
                raise Exception("Only Circle or Ellipse components are plotted")
            ax.add_patch(e)


    # Saving output
    if outfile:
        if outdir is None:
            outdir = '.'
        # If the directory does not exist, create it
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        path = os.path.join(outdir, outfile)
        print "Saving to {}.{}".format(path, ext)
        plt.savefig("{}.{}".format(path, ext), bbox_inches='tight', dpi=200)

    if show:
        fig.show()
    if close:
        plt.close()

    return fig



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