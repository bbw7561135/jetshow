from abc import ABC
import functools
from astropy import units as u, constants as const, cosmology
import numpy as np
from fourier import FINUFFT_NUNU
from pyjetshow import get_i_image as get_i_image_tangled
from pyjetshow_toroidal import get_i_image as get_i_image_toroidal


class JetModelZoom(ABC):
    cosmo = cosmology.WMAP9

    def __init__(self, nu, redshift, pixsize_array=None, stokes="I", tau_max=10**5, central_vfield=True,
                 tangled_bfield=True, ft_class=FINUFFT_NUNU):
        self.nu = nu
        self.z = redshift
        # 2D array of u.angle pixsizes
        self.pixsize_array = pixsize_array
        if pixsize_array is not None:
            assert pixsize_array.ndim == 2, "Array with pixel sizes must be 2D!"
            assert pixsize_array.unit.physical_type == "angle", "pixsize_array" \
                                                                " should be angle!"
        self.stokes = stokes
        self.ft_class = ft_class
        self.ft_instance = None
        self.calculate_grid()
        self._image = None
        self.tau_max = tau_max
        self.central_vfield = central_vfield
        self.tangled_bfield = tangled_bfield
        if tangled_bfield:
            self.get_i_image = get_i_image_tangled
        else:
            self.get_i_image = get_i_image_toroidal

    def calculate_grid(self):
        """
        Calculate grid of ``(r_ob, d)`` - each point is in the center of the
        corresponding pixel. Thus, ``JetModelZoom.ft`` method should not shift
        phases on a half of a pixel
        """
        if self.pixsize_array is None:
            assert self.npix[0] % 2 == 0 and self.npix[1] % 2 == 0
            pc_x = (np.arange(0, self.npix[0])+0.5)*self.pix_to_pc
            pc_y = (np.arange(-self.npix[1]//2, self.npix[1]//2)+0.5)*self.pix_to_pc
            self.r_ob, self.d = np.meshgrid(pc_x, pc_y)
            self.r_ob = self.r_ob.T
            # FIXME: In analogy with older convention we need ``-`` here
            self.d = -self.d.T[:, ::-1]
        else:
            pc_x = np.cumsum(self.pix_to_pc, axis=0)-self.pix_to_pc/2
            pc_y_up = self.pix_to_pc[:, :self.pix_to_pc.shape[1]//2][::-1]
            pc_y_low = self.pix_to_pc[:, self.pix_to_pc.shape[1]//2:]
            pc_y_up = (np.cumsum(pc_y_up, axis=1) - pc_y_up/2)[::-1]
            pc_y_low = np.cumsum(pc_y_low, axis=1) - pc_y_low/2
            pc_y = np.hstack((pc_y_up[:, ::-1], -pc_y_low))
            self.r_ob = pc_x
            # FIXME: In analogy with older convention we need ``-`` here
            self.d = -pc_y

    @property
    def imgsize(self):
        if self.pixsize_array is None:
            return self.pixsize * self.npix[0], self.pixsize * self.npix[1]
        else:
            return np.max(np.cumsum(self.pixsize_array, axis=0)), \
                   np.max(np.cumsum(self.pixsize_array, axis=1))

    @property
    def img_extent(self):
        s0, s1 = self.imgsize[0].value, self.imgsize[1].value
        return 0, s0, -s1/2, s1/2

    @property
    def nparams(self): return 10

    def set_params_vec(self, vec):
        self.dx, self.dy, self.rot, self.theta, self.phi, B1, K1, Gamma, m, n = vec
        self.B1 = np.exp(B1)
        self.K1 = np.exp(K1)
        self.gamma = np.exp(Gamma)
        self.m = np.exp(m)
        self.n = np.exp(n)

    def image(self):
        return np.atleast_2d(self.get_i_image(self.theta, self.z, self.npix[0], self.npix[1],
                                              self.pixsize.value, self.phi, self.B1,
                                              self.m, self.K1, self.n, self.gamma,
                                              (self.nu/u.GHz).value, self.tau_max,
                                              self.central_vfield))

    def ft(self, uv):
        mas_to_rad = u.mas.to(u.rad)
        rot = np.array([[np.cos(self.rot), -np.sin(self.rot)],
                        [np.sin(self.rot), np.cos(self.rot)]])

        # No need in half pixel shifts cause coordinates are already at pixel
        # centers
        shift = [self.dx*mas_to_rad, self.dy*mas_to_rad]
        result = np.exp(-2.0*np.pi*1j*(uv @ shift))
        uv = uv @ rot

        x = (self.d*u.pc/self.ang_to_dist).to(u.rad).value
        y = (self.r_ob*u.pc/self.ang_to_dist).to(u.rad).value
        ft_instance = self.ft_class(uv, x.ravel(), y.ravel())
        img = self.image()
        result *= ft_instance.forward(img.T.ravel())
        del ft_instance, x, y

        return result

    @property
    @functools.lru_cache()
    def ang_to_dist(self):
        return self.cosmo.kpc_proper_per_arcmin(self.z)

    @property
    @functools.lru_cache()
    def pix_to_pc(self):
        """
        2D array of pixel sizes in parsecs or single value.
        """
        return (self.pixsize_array * self.ang_to_dist).to(u.pc).value

    @property
    @functools.lru_cache()
    def T_to_flux(self):
        """
        2D array of coefficients for each pixel or single value if
        ``pixsize_array`` is None.
        """
        omega = self.pixsize_array**2
        mul = (u.K / ((const.c / self.nu)**2 / (2 * const.k_B * omega)) / u.rad**2).to(u.Jy).value
        return mul

    def plot_contours(self, nlevels=15, outfile=None):
        # Factor that accounts non-uniform pixel size in plotting
        # FIXME: self.image have fluxes already fixed by T_to_flux
        factor = (self.pixsize_array/np.min(self.pixsize_array))**2
        image = self.image()/factor.T
        import matplotlib.pyplot as plt
        from matplotlib.ticker import MaxNLocator
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        fig, axes = plt.subplots(1, 1)
        image = np.ma.array(image, mask=image == 0)
        levels = MaxNLocator(nbins=nlevels).tick_values(image.min(),
                                                        image.max())
        cmap = plt.get_cmap('plasma')

        # Contours are *point* based plots, so it is suitable for ``d`` and
        # ``r_ob`` that are centers of pixels.
        cf = axes.contourf(self.r_ob, self.d, image.T, levels=levels, cmap=cmap)
        axes.set_ylabel(r"$d$, pc")
        axes.set_xlabel(r"$r_{\rm ob}$, pc")

        # Make a colorbar with label and units
        divider = make_axes_locatable(axes)
        cax = divider.append_axes("right", size="10%", pad=0.00)
        cb = fig.colorbar(cf, cax=cax)
        # Intensity is in Jy per minimal pixel
        cb.set_label("Intensity, Jy/pix")
        if outfile:
            fig.savefig(outfile, dpi=600, bbox_inches="tight")
        return fig

    def plot(self, nlevels=15, outfile=None):
        # Factor that accounts non-uniform pixel size in plotting
        factor = (self.pixsize_array/np.min(self.pixsize_array))**2
        image = self.image()/factor.T
        import matplotlib.pyplot as plt
        from matplotlib.ticker import MaxNLocator
        from matplotlib.colors import BoundaryNorm
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        fig, axes = plt.subplots(1, 1)
        image = np.ma.array(image, mask=image == 0)
        levels = MaxNLocator(nbins=nlevels).tick_values(image.min(),
                                                        image.max())
        cmap = plt.get_cmap('plasma')
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

        # Here X and Y are 2D arrays of bounds, so ``image`` should be the value
        # *inside* those bounds. Therefore, we should remove the last value from
        # the ``image`` array. Currently we are not doing it.
        im = axes.pcolormesh(self.r_ob, self.d, image.T, norm=norm, cmap=cmap)
        axes.set_ylabel(r"$d$, pc")
        axes.set_xlabel(r"$r_{\rm ob}$, pc")

        # Make a colorbar with label and units
        divider = make_axes_locatable(axes)
        cax = divider.append_axes("right", size="10%", pad=0.00)
        cb = fig.colorbar(im, cax=cax)
        # Intensity is in Jy per minimal pixel
        cb.set_label("Intensity, Jy/pix")
        if outfile:
            fig.savefig(outfile, dpi=600, bbox_inches="tight")
        return fig


if __name__ == "__main__":
    resolutions = np.logspace(-3, -1, 10)
    pixsize_array = np.tile(resolutions, 6).reshape(6, 10).T*u.mas
    cosmo = cosmology.WMAP9
    z = 0.0165
    print((pixsize_array*cosmo.kpc_proper_per_arcmin(z)).to(u.pc))