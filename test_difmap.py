import os
import numpy as np
import sys
sys.path.insert(0, '/home/ilya/github/vlbi_errors/vlbi_errors')
from spydiff import modelfit_difmap
from uv_data import UVData
from components import CGComponent, ImageComponent
from model import Model
from example_simple import gaussian


imsize = 1024

mas_in_pix = 0.01
y = np.arange(-imsize / 2, imsize / 2, dtype=float)
x = np.arange(-imsize / 2, imsize / 2, dtype=float)
imsize = (imsize, imsize)
x *= mas_in_pix
y *= mas_in_pix
xx, yy = np.meshgrid(x, y)
flux = 1.0
bmaj = bmin = 0.2
sigma = bmin/(2.0*np.sqrt(2.0*np.log(2.0)))
amp = flux/(2.0*np.pi*(sigma/mas_in_pix)**2)
g = gaussian(amp, 0, 0, bmaj)
image_g = g(xx, yy)
image_g[image_g < 0.00001] = 0


uvdata = UVData('/home/ilya/github/bck/jetshow/uvf/0235+164.u.2006_06_15.uvf')
mas_to_rad = 4.8481368 * 1E-09
y, z = np.meshgrid(np.arange(imsize[0]), np.arange(imsize[1]))
y = y - imsize[0] / 2. + 0.5
z = z - imsize[0] / 2. + 0.5
y_mas = y * mas_in_pix
z_mas = z * mas_in_pix
y_rad = mas_to_rad * y_mas
z_rad = mas_to_rad * z_mas

icomp = ImageComponent(image_g, y_rad[0, :], z_rad[:, 0])

noise = uvdata.noise(use_V=True)
for key, value in noise.items():
    noise[key] = 0.1 * value
model = Model(stokes='I')
model.add_component(icomp)

# jet_comp = CGComponent(0.5, 1., 0., 0.3)
# model.add_component(jet_comp)

uvdata.substitute([model])
uvdata.noise_add(noise)
uvdata.save('/home/ilya/github/bck/jetshow/uvf/test.fits', rewrite=True)

modelfit_difmap('test.fits', 'initial_cg.mdl', 'out_test.mdl',
                niter=300, path='/home/ilya/github/bck/jetshow/uvf',
                mdl_path='/home/ilya/github/bck/jetshow',
                out_path='/home/ilya/github/bck/jetshow/uvf')

