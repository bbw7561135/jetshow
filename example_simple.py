imsize = (200,200)
mas_in_pix = 0.02
y, z = np.meshgrid(np.arange(imsize[0]), np.arange(imsize[1]))
y = y - imsize[0] / 2. + 0.5
z = z - imsize[0] / 2. + 0.5
y_mas = y*mas_in_pix
z_mas = z*mas_in_pix
from utils import mas_to_rad
y_rad = mas_to_rad * y_mas
z_rad = mas_to_rad * z_mas
yy, xx = np.where(image != 0)
x_ = y_rad[0,:]; y_ = z_rad[:,0]
xx = x_[xx]
yy = y_[yy]
uvdata = UVData('0235+164.x.2006_06_15.uvf')
u = uvdata.uv[:, 0]
v = uvdata.uv[:, 1]

visibilities = list()
for u0, v0 in zip(u, v):
    visibility = (image.compressed() * np.exp(-2.0 * math.pi * 1j *
                                              (u0 * xx + v0 * yy))).sum()
    visibilities.append(visibility)

amps = [np.sqrt(v*v.conjugate()).real for v in visibilities]
phases = [angle(v) for v in visibilities]
uv_rad = list()
for u_, v_ in zip(u,v):
    uv_rad.append(np.sqrt(u_*u_+v_*v_))

plot(uv_rad, amps, '.k')
plot(uv_rad, phases, '.k')


icomp = ImageComponent(image_i, y_rad[0,:], z_rad[:,0])
noise = uvdata.noise(use_V=True)
from model import Model
model = Model(stokes='I')
model.add_component(icomp)
uvdata.substitute([model])
uvdata.uvplot()
