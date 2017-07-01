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
uvdata.noise_add(noise)
uvdata.save('/home/ilya/github/bck/jetshow/bk.fits')
uvdata.uvplot()


from spydiff import modelfit_difmap, clean_difmap, import_difmap_model
from image import plot as iplot
from image import find_bbox
from from_fits import create_clean_image_from_fits_file
from image_ops import rms_image

path_to_script = '/home/ilya/github/vlbi_errors/difmap/final_clean_nw'
data_dir = '/home/ilya/github/bck/jetshow/'
clean_difmap('bk.fits', 'bk_cc.fits', 'I', (1024, 0.1), path=data_dir,
             path_to_script=path_to_script, show_difmap_output=True,
             outpath=data_dir)
modelfit_difmap('bk.fits', 'initial.mdl', 'bk.mdl', niter=100, path=data_dir,
                mdl_path=data_dir, out_path=data_dir)
             
ccimage = create_clean_image_from_fits_file(os.path.join(data_dir, 'bk_cc.fits'))
beam = ccimage.beam
rms = rms_image(ccimage)
blc, trc = find_bbox(ccimage.image, rms, 10)
comps = import_difmap_model('bk.mdl', data_dir)
iplot(ccimage.image, x=ccimage.x, y=ccimage.y, min_abs_level=3*rms, beam=beam,
      show_beam=True, blc=blc, trc=trc, components=comps)


import collections


def deep_update(source, overrides):
    """Update a nested dictionary or similar mapping.
    https://stackoverflow.com/a/30655448

    Modify ``source`` in place.
    """
    for key, value in overrides.iteritems():
        if isinstance(value, collections.Mapping) and value:
            returned = deep_update(source.get(key, {}), value)
            source[key] = returned
        else:
            source[key] = overrides[key]
    return source
    
    

def update_config(cfg_in, update_dict, cfg_out=None):
    """Update json config file with some nested structure e.g. change frequency
    to 2.0 GHz:

    ``
    update_dict = {u'observation' : {u'frequency_ghz': 2.0}}
    ``
    """
    with open(cfg_in, "r") as jsonFile:
        data = json.load(jsonFile)
    deep_update(data, update_dict)
    if cfg_out is None:
        cfg_out = cfg_in
    with open(cfg_out, "w") as jsonFile:
        json.dump(data, jsonFile)
    

