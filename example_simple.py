import os
os.chdir('/home/ilya/github/vlbi_errors/vlbi_errors')
from components import ImageComponent
from uv_data import UVData
from model import Model
from utils import mas_to_rad
uvdata = UVData('0235+164.x.2006_06_15.uvf')
imsize = (500,500)
mas_in_pix = 0.002
y, z = np.meshgrid(np.arange(imsize[0]), np.arange(imsize[1]))
y = y - imsize[0] / 2. + 0.5
z = z - imsize[0] / 2. + 0.5
y_mas = y*mas_in_pix
z_mas = z*mas_in_pix
y_rad = mas_to_rad * y_mas
z_rad = mas_to_rad * z_mas
image_i = np.loadtxt('/home/ilya/github/bck/jetshow/cmake-build-debug/map_i.txt')
image_tau = np.loadtxt('/home/ilya/github/bck/jetshow/cmake-build-debug/map_tau.txt')
matshow(image_i);colorbar()
contour(log10(image_tau), [-1, -0.5, 0, 0.5, 1, 2], cmap='tab10');colorbar()
icomp = ImageComponent(image_i, y_rad[0,:], z_rad[:,0])
noise = uvdata.noise(use_V=True)
model = Model(stokes='I')
model.add_component(icomp)
uvdata.substitute([model])
uvdata.uvplot()

