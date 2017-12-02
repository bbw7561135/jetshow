import os
import sys
sys.path.insert(0, '/home/ilya/github/vlbi_errors/vlbi_errors')
import numpy as np
from skimage.transform import rotate
import matplotlib.pyplot as plt
from uv_data import UVData
from components import ImageComponent
from model import Model


mas_to_rad = 4.8481368 * 1E-09
# uv_file = '/home/ilya/github/bck/jetshow/uvf/0716+714_raks01xg_C_LL_0060s_uva.fits'
uv_file = '/home/ilya/github/bck/jetshow/uvf/2200+420_K_SVLBI.uvf'
uvdata = UVData(uv_file)

# fig = uvdata.uvplot(stokes=["LL"])
fig = uvdata.uvplot()

images = list()
angles = range(0, 180, 30)
# image = '/home/ilya/github/bck/jetshow/uvf/map_i_09_C.txt'
image = '/home/ilya/github/bck/jetshow/cmake-build-debug/map_i.txt'
image = np.loadtxt(image)
images.append(image)

# imsize = 1096
imsize = 1734
imsize = (imsize, imsize)
# mas_in_pix = 0.005
mas_in_pix = 0.00253
y, z = np.meshgrid(np.arange(imsize[0]), np.arange(imsize[1]))
y = y - imsize[0] / 2. + 0.5
z = z - imsize[0] / 2. + 0.5
y_mas = y * mas_in_pix
z_mas = z * mas_in_pix
y_rad = mas_to_rad * y_mas
z_rad = mas_to_rad * z_mas

for angle in angles:
    rotated_image = rotate(image, angle, order=1)
    images.append(rotated_image)
# plt.matshow(image); colorbar()
# plt.matshow(rotated_image); colorbar()

colors = [a[u'color'] for a in list(plt.rcParams[u'axes.prop_cycle'])]
for i, image in enumerate(images):
    image[image < 0] = 0
    image[image > 10.0] = 0
    icomp = ImageComponent(image, y_rad[0, :], z_rad[:, 0])
    model = Model(stokes='I')
    model.add_component(icomp)
    uvdata.substitute([model])
    fig = uvdata.uvplot(fig=fig, color=colors[i], alpha=0.05)

# fig.savefig('/home/ilya/github/bck/jetshow/uvf_mf_adds/ra.png',
#             bbox_inches='tight', dpi=300)