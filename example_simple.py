import os
import sys
sys.path.insert(0, '/home/ilya/github/vlbi_errors/vlbi_errors')
import subprocess
import collections
import json
import numpy as np
from components import ImageComponent
from uv_data import UVData
from model import Model
from utils import mas_to_rad


# uvdata = UVData('0235+164.x.2006_06_15.uvf')
# imsize = (500,500)
# mas_in_pix = 0.002
# y, z = np.meshgrid(np.arange(imsize[0]), np.arange(imsize[1]))
# y = y - imsize[0] / 2. + 0.5
# z = z - imsize[0] / 2. + 0.5
# y_mas = y*mas_in_pix
# z_mas = z*mas_in_pix
# y_rad = mas_to_rad * y_mas
# z_rad = mas_to_rad * z_mas
# image_i = np.loadtxt('/home/ilya/github/bck/jetshow/cmake-build-debug/map_i.txt')
# image_tau = np.loadtxt('/home/ilya/github/bck/jetshow/cmake-build-debug/map_tau.txt')
# matshow(image_i);colorbar()
# contour(log10(image_tau), [-1, -0.5, 0, 0.5, 1, 2], cmap='tab10');colorbar()
# icomp = ImageComponent(image_i, y_rad[0,:], z_rad[:,0])
# noise = uvdata.noise(use_V=True)
# model = Model(stokes='I')
# model.add_component(icomp)
# uvdata.substitute([model])
# uvdata.noise_add(noise)
# uvdata.save('/home/ilya/github/bck/jetshow/bk.fits')
# uvdata.uvplot()


from spydiff import modelfit_difmap, clean_difmap, import_difmap_model
from image import plot as iplot
from image import find_bbox
from from_fits import create_clean_image_from_fits_file
from image_ops import rms_image

# path_to_script = '/home/ilya/github/vlbi_errors/difmap/final_clean_nw'
# data_dir = '/home/ilya/github/bck/jetshow/'
# clean_difmap('bk.fits', 'bk_cc.fits', 'I', (1024, 0.1), path=data_dir,
#              path_to_script=path_to_script, show_difmap_output=True,
#              outpath=data_dir)
# modelfit_difmap('bk.fits', 'initial.mdl', 'bk.mdl', niter=100, path=data_dir,
#                 mdl_path=data_dir, out_path=data_dir)
#
# ccimage = create_clean_image_from_fits_file(os.path.join(data_dir, 'bk_cc.fits'))
# beam = ccimage.beam
# rms = rms_image(ccimage)
# blc, trc = find_bbox(ccimage.image, rms, 10)
# comps = import_difmap_model('bk.mdl', data_dir)
# iplot(ccimage.image, x=ccimage.x, y=ccimage.y, min_abs_level=3*rms, beam=beam,
#       show_beam=True, blc=blc, trc=trc, components=comps)
#



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


def run_simulations(cfg_file, out_dir, initial_dfm_model,
                    path_to_executable, uv_fits_template,
                    uv_fits_save_fname='bk.fits', out_dfm_model_fn='bk.mdl',
                    update_params_dict=None):
    """
    Run simulation of BK jet, fit results.

    :param cfg_file:
        Path to configuration file in JSON format.
    :param out_dir:
        Directory to store results.
    :param initial_dfm_model:
        Path to difmap model used as initial guess.
    :param path_to_executable:
        Path to compiled executable of ``jetshow``.
    :param uv_fits_template:
        Path to UV-FITS file with data being substituted by results of the
        simulations.
    :param uv_fits_save_fname: (optional)
        Filename of UV-FITS file with artificial data that will be saved in
        ``out_dir``. (default: ``bk.fits``)
    :param out_dfm_model_fn: (optional)
        Filename of difmap model of artificial data hat will be saved in
        ``out_dir``. (default: ``bk.mdl``)
    :param update_params_dict: (optional)
        Dictionary (possible nested) with updates of the parameters that
        corresponds to config.json configuration file.
    """
    if update_params_dict is not None:
        update_config(cfg_file, update_params_dict)
    # Load parameters
    with open(cfg_file, "r") as jsonFile:
        params = json.load(jsonFile)

    exe_dir, exe = os.path.split(path_to_executable)
    cwd = os.getcwd()
    os.chdir(exe_dir)
    subprocess.call(["./{}".format(exe)])
    os.chdir(cwd)

    uv_fits_template_dir, uv_fits_template_fn = os.path.split(uv_fits_template)
    uvdata = UVData(uv_fits_template)
    imsize = params[u'image'][u'number_of_pixels']
    imsize = (imsize, imsize)
    mas_in_pix = params[u'image'][u'pixel_size_mas']
    y, z = np.meshgrid(np.arange(imsize[0]), np.arange(imsize[1]))
    y = y - imsize[0] / 2. + 0.5
    z = z - imsize[0] / 2. + 0.5
    y_mas = y*mas_in_pix
    z_mas = z*mas_in_pix
    y_rad = mas_to_rad * y_mas
    z_rad = mas_to_rad * z_mas
    image_i_file = os.path.join(exe_dir, params[u'output'][u'file_i'])
    image_tau_file = os.path.join(exe_dir, params[u'output'][u'file_tau'])
    image_i = np.loadtxt(image_i_file)
    image_tau = np.loadtxt(image_tau_file)
    icomp = ImageComponent(image_i, y_rad[0,:], z_rad[:,0])
    noise = uvdata.noise(use_V=True)
    model = Model(stokes='I')
    model.add_component(icomp)
    uvdata.substitute([model])
    uvdata.noise_add(noise)
    uvdata.save(os.path.join(out_dir, uv_fits_save_fname))

    initial_dfm_dir, initial_dfm_model_fn = os.path.split(initial_dfm_model)
    # Fit model in difmap
    modelfit_difmap(uv_fits_save_fname, initial_dfm_model_fn, out_dfm_model_fn,
                    niter=100, path=out_dir,
                    mdl_path=initial_dfm_dir,
                    out_path=out_dir)



if __name__ == '__main__':
    main_dir = '/home/ilya/github/bck/jetshow'
    out_dir = os.path.join(main_dir, 'test')
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    cfg_file = os.path.join(main_dir, 'config.json')
    initial_dfm_model = os.path.join(main_dir, 'initial_eg.mdl')
    path_to_executable = os.path.join(main_dir, 'cmake-build-debug', 'jetshow')
    uv_fits_template = '/home/ilya/github/vlbi_errors/vlbi_errors/'
    uv_fits_template = os.path.join(uv_fits_template,
                                    '0235+164.x.2006_06_15.uvf')
    run_simulations(cfg_file, out_dir, initial_dfm_model,
                    path_to_executable, uv_fits_template,
                    uv_fits_save_fname='bk.fits', out_dfm_model_fn='bk.mdl',
                    update_params_dict=None)