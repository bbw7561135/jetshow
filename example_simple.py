import os
import sys
sys.path.insert(0, '/home/ilya/github/vlbi_errors/vlbi_errors')
import subprocess
import collections
import json
import numpy as np
from scipy.optimize import curve_fit
from components import ImageComponent, CGComponent
from uv_data import UVData
from model import Model
from utils import mas_to_rad
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting


mas_to_rad = 4.8481368 * 1E-09
rad_to_mas = 1. / mas_to_rad


class FailedFindBestImageParamsException(Exception):
    pass

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


def analyze_tau_stripe(fname, log10_tau_small=-4.0, border_tau1=(0.001, 0.3),
                       border_tau_min=(0.5, 1), min_pixels_to_tau1=20):
    tau = np.loadtxt(fname)
    length = len(tau)
    idx_tau1 = np.argmin(np.abs(np.log10(tau) - 0.0))
    idx_tau_ = np.argmin(np.abs(np.log10(tau) - log10_tau_small))

    increase_pixel_size = False
    decrease_pixels_size = False
    increase_number_of_pixels = False
    decrease_number_of_pixels = False

    if float(idx_tau1)/length > border_tau1[0]:
        if float(idx_tau1)/length > border_tau1[1]:
            increase_pixel_size = True
    else:
        decrease_pixels_size = True

    if idx_tau_ < length*border_tau_min[1] - 1:
        if idx_tau_ < border_tau_min[0] * length:
            decrease_number_of_pixels = True
    else:
        increase_number_of_pixels = True

    if idx_tau1 < min_pixels_to_tau1:
        decrease_pixels_size = True

    return {"increase_pixel_size": increase_pixel_size,
            "decrease_pixels_size": decrease_pixels_size,
            "increase_number_of_pixels": increase_number_of_pixels,
            "decrease_number_of_pixels": decrease_number_of_pixels}


def run_jetshow(cfg_file, path_to_executable, update_params_dict=None):
    """
    Run simulation of BK jet.

    :param cfg_file:
        Path to configuration file in JSON format.
    :param path_to_executable:
        Path to compiled executable of ``jetshow``.
    :param update_params_dict: (optional)
        Dictionary (possible nested) with updates of the parameters that
        corresponds to config.json configuration file. If ``None`` then doesn't
        update configuration JSON-file. (default: ``None``)
    :return:
        Dictionary with parameters used.
    """
    if update_params_dict is not None:
        update_config(cfg_file, update_params_dict)
    # Load parameters
    with open(cfg_file, "r") as jsonFile:
        params = json.load(jsonFile)
    exe_dir, exe = os.path.split(path_to_executable)
    cwd = os.getcwd()
    os.chdir(exe_dir)
    print("Running jetshow with"
          " parameters N={}, size={} mas".format(params[u'image'][u'number_of_pixels'],
                                                 params[u'image'][u'pixel_size_mas']))
    print("Calculating mode = {}".format(params[u'calculate']))
    subprocess.call(["./{}".format(exe)])
    os.chdir(cwd)
    return params


def find_image_params(cfg_file, path_to_executable):
    """
    Iteratively find parameters of image (number of pixels & size of a pixel)
    that makes image of jet with physical parameters specified in ``cfg_file``
    looks good.

    :param cfg_file:
        Path to configuration file in JSON format.
    :param path_to_executable:
        Path to compiled executable of ``jetshow``.
    :return:
        Number of pixels & size of a single pixel [mas].
    """
    exe_dir, exe = os.path.split(path_to_executable)

    with open(cfg_file, "r") as jsonFile:
        params = json.load(jsonFile)

    # First, calculate stripe of optical depth
    try:
        os.unlink(os.path.join(exe_dir, 'stripe_tau.txt'))
    except OSError:
        pass
    params = run_jetshow(cfg_file, path_to_executable,
                         update_params_dict={"calculate": "tau"})
    # Check resulting stripe and find what we need to do with image parameters
    # to make image looking nice
    decision_dict = analyze_tau_stripe(os.path.join(exe_dir, 'stripe_tau.txt'))
    updated = False

    # Until image parameters are OK for us
    updates_steps_size = np.linspace(1.5, 2, 20)[::-1]
    updates_steps_n = np.linspace(1, 1.25, 20)[::-1]
    update_step_n = 0
    update_step_size = 0

    while True in decision_dict.values():
        if updated:
            decision_dict = analyze_tau_stripe(os.path.join(exe_dir,
                                                            'stripe_tau.txt'))
        if decision_dict["increase_pixel_size"]:
            params[u'image'][u'pixel_size_mas'] *= updates_steps_size[update_step_size]
            update_step_size += 1
        elif decision_dict["decrease_pixels_size"]:
            params[u'image'][u'pixel_size_mas'] /= updates_steps_size[update_step_size]
            update_step_size += 1
        elif decision_dict["increase_number_of_pixels"]:
            new_number = int(updates_steps_n[update_step_n]*params[u'image'][u'number_of_pixels'])
            update_step_n += 1
            if new_number % 2:
                new_number += 1
            params[u'image'][u'number_of_pixels'] = new_number
        elif decision_dict["decrease_number_of_pixels"]:
            new_number = int(params[u'image'][u'number_of_pixels']/updates_steps_n[update_step_n])
            update_step_n += 1
            if new_number % 2:
                new_number += 1
            params[u'image'][u'number_of_pixels'] = new_number
        else:
            updated = False
        updated =True

        os.unlink(os.path.join(exe_dir, 'stripe_tau.txt'))
        params = run_jetshow(cfg_file, path_to_executable,
                             update_params_dict=params)

    return (params[u'image'][u'number_of_pixels'],
            params[u'image'][u'pixel_size_mas'])


def run_simulations(cfg_file, path_to_executable, map_size=None):
    """
    Run simulation of BK jet, fit results.

    :param cfg_file:
        Path to configuration file in JSON format.
    :param path_to_executable:
        Path to compiled executable of ``jetshow``.
    :param map_size: (optional)
        Iterable of number of pixels and pixel size in map. If ``None`` then
        determine using optical depth calculation. (default: ``None``)
    """
    exe_dir, exe = os.path.split(path_to_executable)
    if map_size is not None:
        number_of_pixels, pixel_size_mas = map_size
    else:
        try:
            number_of_pixels, pixel_size_mas = find_image_params(cfg_file,
                                                                 path_to_executable)
        # Hope second pass will help
        except IndexError:
            try:
                number_of_pixels, pixel_size_mas = find_image_params(cfg_file,
                                                                     path_to_executable)
            except IndexError:
                raise FailedFindBestImageParamsException

    update_params_dict = {"calculate": "full",
                          "image": {"number_of_pixels": number_of_pixels,
                                    "pixel_size_mas": pixel_size_mas}}
    update_config(cfg_file, update_params_dict)
    if map_size is None:
        try:
            os.unlink(os.path.join(exe_dir, 'stripe_tau.txt'))
        except OSError:
            pass
    params = run_jetshow(cfg_file, path_to_executable)
    return params


def modelfit_simulation_result(exe_dir, initial_dfm_model, noise_factor,
                               out_dfm_model_fn, out_dir,
                               params, uv_fits_save_fname,
                               uv_fits_template, niter=300):
    """
    Modelfit result of simulations FT to uv-plane by substituting some real
    data.

    :param exe_dir:
        Directory with results of the simulations (it should contain exe-file).
    :param initial_dfm_model:
        Initial model used in difmap fitting.
    :param noise_factor:
        Float value that scales noise from real data set added to FT of
        simulation result.
    :param out_dfm_model_fn:
        File name of the modelfit result txt-file.
    :param out_dir:
        Directory to save modelfit result txt-file.
    :param params:
        Dictionary with parameters used for simulation.
    :param uv_fits_save_fname:
        File name of generated uv-data to save.
    :param uv_fits_template:
        Path to fits file used as template for substituting by FT of simulation
        result.
    :param niter: (optional)
        Number of iterations for difmap fitting. (default: ``300``)
    """
    uvdata = UVData(uv_fits_template)

    # Create coordinate grid
    imsize = params[u'image'][u'number_of_pixels']
    imsize = (imsize, imsize)
    mas_in_pix = params[u'image'][u'pixel_size_mas']
    y, z = np.meshgrid(np.arange(imsize[0]), np.arange(imsize[1]))
    y = y - imsize[0] / 2. + 0.5
    z = z - imsize[0] / 2. + 0.5
    y_mas = y * mas_in_pix
    z_mas = z * mas_in_pix
    y_rad = mas_to_rad * y_mas
    z_rad = mas_to_rad * z_mas

    image_i_file = os.path.join(exe_dir, params[u'output'][u'file_i'])
    image_i = np.loadtxt(image_i_file)
    image_i[image_i < 0] = 0
    image_i[image_i > 10.0] = 0
    icomp = ImageComponent(image_i, y_rad[0, :], z_rad[:, 0])

    noise = uvdata.noise(use_V=True)
    for key, value in noise.items():
        noise[key] = noise_factor * value
    model = Model(stokes='I')
    model.add_component(icomp)

    # jet_comp = CGComponent(0.5, 1., 0., 0.3)
    # model.add_component(jet_comp)

    uvdata.substitute([model])
    uvdata.noise_add(noise)
    uvdata.save(os.path.join(out_dir, uv_fits_save_fname), rewrite=True)
    initial_dfm_dir, initial_dfm_model_fn = os.path.split(initial_dfm_model)
    # Fit model in difmap
    modelfit_difmap(uv_fits_save_fname, initial_dfm_model_fn, out_dfm_model_fn,
                    niter=niter, path=out_dir, mdl_path=initial_dfm_dir,
                    out_path=out_dir)


def nested_dict():
    return collections.defaultdict(nested_dict)


def shift_model(nu, a, k):
    return a * nu**(-1./k)


def shift_model_dr(nu, a, k, b):
    return a * nu**(-1./k) + b


def n1_equipartition(b1, gama_min=1):
    """
    Equipartition particle density for given magnetic field.

    :param b1:
        Magnetic field value on 1 pc [G].
    :param gama_min: (optional)
        Minimum particles lorentz factor. (default: ``1``)
    :return:
        Value of particle density [cm**(-3)] at 1 pc.

    :note:
        Assumptions used:

            * ``n = 2*m``
            * ``\alpha=-0.5``, ``gamma_max = 10**4.34 * gamma_min`` - needs for
              ``K = 1/log(gamma_max/gamma_min) = 0.1``. See 2005ApJ...619...73H
    """
    return 0.1*b1**2/(gama_min*9.11*10**(-28)*9.0*10**20*8*np.pi)


def b1_equipartition(n1, gama_min=1):
    """
    Equipartition particle density for given magnetic field.

    :param n1:
        Particle density value [cm**(-3)] on 1 pc.
    :param gama_min: (optional)
        Minimum particles lorentz factor. (default: ``1``)
    :return:
        Magnetic field value on 1 pc [G].

    :note:
        Assumptions used:

            * ``n = 2*m``
            * ``\alpha=-0.5``, ``gamma_max = 10**4.34 * gamma_min`` - needs for
              ``K = 1/log(gamma_max/gamma_min) = 0.1``. See 2005ApJ...619...73H
    """
    return np.sqrt(n1*gama_min*9.11*10**(-28)*9.0*10**20*8*np.pi/0.1)


def log_base(x, base):
    return np.log(x)/np.log(base)


def b_to_n_energy_ratio(b1, n1, gamma_min=1):
    """
    Ratio of magnetic energy density to particle energy density.

    :param b1:
        Magnetic field value on 1 pc [G].
    :param n1:
        Particle density value [cm**(-3)] on 1 pc.
    :param gamma_min: (optional)
        Minimum particles lorentz factor. (default: ``1``)
    :return:
        Value of ratio.

    :note:
    Assumptions used:

        * ``n = 2*m``
        * ``\alpha=-0.5``, ``gamma_max = 10**4.34 * gamma_min`` - needs for
          ``K = 1/log(gamma_max/gamma_min) = 0.1``. See 2005ApJ...619...73H
    """
    return 0.1*b1**2/(8*np.pi*n1*gamma_min*9.11*10**(-28)*9.0*10**20)


def t_syn(b, nu, D=1.0):
    """
    Synchrotron lifetime for given frequency.

    :param b:
    :param D:
    :param nu:
    :return:
        Lifetime [s].
    """
    # Thompson cross-section [cm**2]
    sigma_t = 6.65 * 10**(-25)
    m_e = 9.109382*1E-28
    q_e = 4.8*1E-10
    c = 3.0*10**10
    return 3./sigma_t * np.sqrt(2*np.pi*c*m_e*q_e/(b**3*D))*nu**(-0.5)


def tb(flux, freq, size, z=0., D=1.):
    """
    Brightness temperature.

    :param flux:
        Flux in Jy.
    :param freq:
        Frequency in GHz.
    :param size:
        Size in mas.
    :return:
        Value of true brightness temperature (corrected for Doppler and
        redhsift).
    """
    k = 1.38 * 10**(-16)
    c = 3.0 * 10 ** 10
    mas_to_rad = 4.8481368 * 1E-09
    freq *= 10**9
    size *= mas_to_rad
    flux *= 10**(-23)
    Tb = c**2*flux/(2.*np.pi*k*size**2*freq**2)
    return (1.+z)*Tb/D


def tb_comp(flux, bmaj, freq, z=0, bmin=None, D=1):
    mas_to_rad = 4.8481368 * 1E-09
    c = 3.0 * 10 ** 10
    k = 1.38 * 10 ** (-16)
    bmaj *= mas_to_rad
    if bmin is None:
        bmin = bmaj
    else:
        bmin *= mas_to_rad
    freq *= 10**9
    flux *= 10**(-23)
    return 2.*np.log(2)*(1.+z)*flux*c**2/(freq**2*np.pi*k*bmaj*bmin*D)


def _find_shift_from_difmap_models(freq_difmap_models_dict):
    """
    Find shift using difmap model files of core.

    :param freq_difmap_models_dict:
        Dictionary with keys - frequency [GHz] and values - paths to difmap
        model files with core.
    :return:
        Dictionary with keys - names of the value used to measure shift and
        values - tuple of measured coefficient ``k`` and per-frequency values.
    """
    drs = list()
    bmajs = list()
    bmins = list()
    freqs = list()
    from collections import OrderedDict
    result_dict = OrderedDict()
    for freq, difmap_model in freq_difmap_models_dict.items():
        difmap_dir, difmap_fn = os.path.split(difmap_model)
        comps = import_difmap_model(difmap_fn, difmap_dir)
        core = comps[0]
        # Distance from SMBH (0,0) in mas
        dr = core.p[1]
        bmaj = core.p[3]
        bmin = bmaj * core.p[4]
        drs.append(dr)
        bmajs.append(bmaj)
        bmins.append(bmin)
        freqs.append(freq)
    for name, container in zip(('dr', 'bmaj', 'bmin'), (drs, bmajs, bmins)):
        try:
            res = curve_fit(shift_model, freqs, container, p0=[1.0, 1.0])
        except RuntimeError:
            res=([0.0, 0.0], None)
        result_dict[name] = (res[0], container)
    return result_dict


def plot_stripe(sim_fname, difmap_model, simulations_params, g=None,
                savefig=None, close=False):
    """
    Plot 1D stripe along jet axis.

    :param sim_fname:
        Path to file with simulations result of Stokes I.
    :param difmap_model:
        File with model in difmap format.
    :param simulations_params:
        Dictionary with simulations parameters.
    """
    difmap_dir, difmap_fn = os.path.split(difmap_model)
    comps = import_difmap_model(difmap_fn, difmap_dir)
    core = comps[0]

    image = np.loadtxt(sim_fname)
    image[image < 0] = 0
    image[image > 10.0] = 0
    imsize = simulations_params[u'image'][u'number_of_pixels']
    mas_in_pix = simulations_params[u'image'][u'pixel_size_mas']

    if g is not None:
        y = np.arange(-imsize/2, imsize/2, dtype=float)
        x = np.arange(-imsize/2, imsize/2, dtype=float)
        x *= mas_in_pix
        y *= mas_in_pix
        xx, yy = np.meshgrid(x, y)
        g_image = g(xx, yy)

    if len(core) == 6:
        e = core.p[4]
    else:
        e = 1.
    # put 2 instead of 4 before sqrt and remove 0.5 before p[3]
    amp = core.p[0] / (2.*np.pi*e*(core.p[3]/(2.*np.sqrt(2.*np.log(2))*mas_in_pix))**2)
    p = [amp, core.p[1], core.p[2], core.p[3]] + list(core.p[4:])
    gaus = gaussian(*p)
    gaus_image = gaus(xx, yy)

    x = np.arange(-20, imsize/2)*mas_in_pix
    # Plot simulation results
    plt.figure()
    plt.plot(x, image[imsize/2, imsize/2-20:], lw=2, label="Simulations")
    # Plot difmap model
    plt.plot(x, gaus_image[imsize/2, imsize/2-20:], label="Difmap fit")
    # Optionally plot model fitted in image plane
    if g is not None:
        plt.plot(x, g_image[imsize/2, imsize/2-20:],
                 label="Sim. image fit")
    # # Distance from SMBH (0,0) in pixels
    # dr = (np.unravel_index(image.argmax(), image.shape)[1]-imsize/2)
    # # In mas
    # dr *= mas_in_pix
    plt.xlabel("Distance along jet, [mas]")
    plt.ylabel("Flux, [Jy/pixel]")
    plt.legend(loc="best")
    plt.tight_layout()
    if savefig:
        plt.savefig(savefig)
    if close:
        plt.close()


def find_core_separation_from_jet_using_difmap_model(difmap_model):
    """
    Find separation between core and jet component using difmap modelfit result
    file.

    :param difmap_model:
        File with model in difmap format.
    :return:
        Value of distance between core and jet component. If only core is
        present then distance between phase center and core
    """

    difmap_dir, difmap_fn = os.path.split(difmap_model)
    comps = import_difmap_model(difmap_fn, difmap_dir)
    core = comps[0]
    try:
        jet = comps[1].p
    except IndexError:
        jet = [0, 0, 0]
    dr = (jet[1] - core.p[1])**2. + (jet[2] - core.p[2])**2.
    return np.sqrt(dr)


def find_core_separation_from_center_using_simulations(fname,
                                                       simulations_params):
    """
    Find separation of core from center (SMBH).

    :param fname:
        Path to file with simulations result of Stokes I.
    :param simulations_params:
        Dictionary with simulations parameters.
    :return:
        Value of distance between core and center (SMBH).
    """
    image = np.loadtxt(fname)
    imsize = simulations_params[u'image'][u'number_of_pixels']
    mas_in_pix = simulations_params[u'image'][u'pixel_size_mas']
    # Distance from SMBH (0,0) in pixels
    dr = (np.unravel_index(image.argmax(), image.shape)[1]-imsize/2)
    # In mas
    dr *= mas_in_pix
    return dr


def plot_simulations_3d(sim_fname, simulations_params, each=2, delta=100,
                        contr_delta=10, core=None):
    """
    Plot simulation results in 3D projection.

    :param sim_fname:
        Path to file with simulations result of Stokes I.
    :param simulations_params:
        Dictionary with simulations parameters.
    :param each: (optional)
        Plot each ``each`` pixel. (default: ``2``)
    :param delta: (optional)
        Cut image from sides. Image will be shorter by ``2*delta`` pixels from
        sides. (default: ``100``)
    :param contr_delta: (optional)
        Space to leave on contr-jet side. (default: ``10``)
    :param core: (optional)
        Instance of ``Component`` class to overplot. If ``None`` then don't plot
        component. (default: ``None``)

    :return:
        Instance of ``Figure``.
    """
    # This import is needed for 3D projections
    from mpl_toolkits.mplot3d import Axes3D
    # Making transparent color map
    from matplotlib import cm
    theCM = cm.Blues
    # theCM._init()
    # alphas = np.abs(np.linspace(-1.0, 1.0, theCM.N))
    # theCM._lut[:-3, -1] = alphas

    image = np.loadtxt(sim_fname)
    print(sum(image))
    imsize = simulations_params[u'image'][u'number_of_pixels']
    mas_in_pix = simulations_params[u'image'][u'pixel_size_mas']
    y = np.arange(-imsize/2+delta, imsize/2-delta, each, dtype=float)
    x = np.arange(-contr_delta, imsize/2, each, dtype=float)
    x *= mas_in_pix*each
    y *= mas_in_pix*each
    xx, yy = np.meshgrid(x, y)

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    plt.hold(True)

    # mask = image < 0.0001
    # image = np.ma.array(image, mask=mask)
    surf = ax.plot_surface(xx, yy,
                           image[delta:-delta:each, imsize/2-contr_delta::each],
                           rstride=1, cstride=1, cmap=theCM, linewidth=2,
                           antialiased=True)
    fig.colorbar(surf, shrink=0.5, aspect=5)

    if core is not None:
        theCM = cm.Reds
        # theCM._init()
        # alphas = np.abs(np.linspace(-1.0, 1.0, theCM.N))
        # theCM._lut[:-3, -1] = alphas
        amp = core.p[0] / (2.*np.pi*(core.p[3]/(4.*np.sqrt(2.*np.log(2))*each*mas_in_pix))**2)
        p = [amp, core.p[1], core.p[2], core.p[3]/2] + list(core.p[4:])
        gaus = gaussian(*p)
        gaus_image = gaus(xx, yy)
        print(sum(gaus_image))
        # mask = gaus_image < 0.0001
        # gaus_image = np.ma.array(gaus_image, mask=mask)
        ax.plot_surface(xx, yy, gaus_image, rstride=1, cstride=1, cmap=theCM,
                        linewidth=2, antialiased=True)
    # ax.set_zlim([0.0002, None])
    plt.show()
    return fig


def plot_simulations_2d(sim_fname, simulations_params, each=1, side_delta=None,
                        jet_delta=None, contr_delta=10, core=None, g=None,
                        savefig=None, close=False):
    """
    Plot simulation results in 2D projection.

    :param sim_fname:
        Path to file with simulations result of Stokes I.
    :param simulations_params:
        Dictionary with simulations parameters.
    :param each: (optional)
        Plot each ``each`` pixel. (default: ``2``)
    :param side_delta: (optional)
        Cut image from sides. Image will be shorter by ``2*delta`` pixels from
        sides. (default: ``100``)
    :param contr_delta: (optional)
        Space to leave on contr-jet side. (default: ``10``)
    :param core: (optional)
        Instance of ``Component`` class to overplot. If ``None`` then don't plot
        component. (default: ``None``)

    :return:
        Instance of ``Figure``.
    """
    # Making transparent color map
    from matplotlib import cm
    theCM = cm.Blues
    image = np.loadtxt(sim_fname)
    image[image < 0] = 0
    image[image > 10.0] = 0
    low = np.percentile(image[image > 0].flatten(), 50)
    high = np.percentile(image[image > 0].flatten(), 99.99)
    levels = np.logspace(np.log2(low), np.log2(high), 10, base=2)
    imsize = image.shape[0]

    blc, trc = find_bbox(image, low)
    if side_delta is None:
        side_delta = int(np.mean([blc[1], -trc[1]+imsize]))-50
    if jet_delta is None:
        jet_delta = imsize - trc[0]

    print("Flux of simulated image :")
    print(image.sum())
    imsize = simulations_params[u'image'][u'number_of_pixels']
    mas_in_pix = simulations_params[u'image'][u'pixel_size_mas']
    y = np.arange(-imsize/2+side_delta, imsize/2-side_delta, each, dtype=float)
    x = np.arange(-contr_delta, imsize/2-jet_delta, each, dtype=float)
    x *= mas_in_pix*each
    y *= mas_in_pix*each
    xx, yy = np.meshgrid(x, y)

    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')
    plt.hold(True)

    cont = ax.contour(xx, yy,
                      image[side_delta:-side_delta:each, imsize/2-contr_delta:imsize-jet_delta:each],
                      levels=levels, colors="blue", label="simulations")
    # fig.colorbar(cont, shrink=0.5, aspect=5)
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="10%", pad=0.00)
    cb = fig.colorbar(cont, cax=cax)
    cb.set_label("Jy/pixel")

    if core is not None:
        theCM = cm.viridis
        if len(core) == 6:
            e = core.p[4]
        else:
            e = 1.
        # put 2 instead of 4 before sqrt and remove 0.5 before p[3]
        amp = core.p[0] / (2.*np.pi*e*(core.p[3]/(2.*np.sqrt(2.*np.log(2))*each*mas_in_pix))**2)
        p = [amp, core.p[1], core.p[2], core.p[3]] + list(core.p[4:])
        gaus = gaussian(*p)
        gaus_image = gaus(xx, yy)
        print("Flux of difmap model image (should be = {}) : ".format(core.p[0]))
        print(gaus_image.sum())
        ax.contour(xx, yy, gaus_image, levels=levels, colors="green",
                   label="difmap")

    if g is not None:
        theCM = cm.viridis
        g_image = g(xx, yy)
        print("Flux of fit of simulated image :")
        print(g_image.sum())
        ax.contour(xx, yy, g_image, levels=levels, colors="red",
                   label="sim. image fit")

    ax.set_xlabel("Distance along jet, [mas]")
    ax.set_ylabel("Distance, [mas]")
    # ax.set_title("")
    fig.tight_layout()
    plt.show()

    if savefig:
        fig.savefig(savefig)
    if close:
        plt.close()

    return fig


def generate_sample_points(uniform_borders, n_samples=100):
    """
    Generate
    :param uniform_borders:
        Iterable of two-elements containers with borders for uniform
        distributions.
    :param n_samples: (optional)
        Number of samples to get. (default: ``100``)
    :return:
        2D numpy array with shape ``(n_samples, n_dim)`` where ``n_dim`` -
        number of dimensions (this corresponds to number of borders in
        ``uniform_borders``).
    """
    import mcerp
    import scipy.stats as ss
    distributions = list()
    for borders in uniform_borders:
        # Should be iterable with 2 elements (min & max)
        assert len(borders) == 2
        # Low border goes first
        assert borders[0] < borders[1]
        distributions.append(ss.uniform(loc=borders[0],
                                        scale=borders[1]-borders[0]))
    return mcerp.lhd(dist=distributions, size=n_samples)


def find_bbox(array, level, delta=0.):
    """
    Find bounding box for part of image containing source.

    :param array:
        Numpy 2D array with image.
    :param level:
        Level at which threshold image in image units.
    :param delta:
        Extra space to add symmetrically [pixels].
    :return:
        Tuples of BLC & TRC.
    """
    from scipy.ndimage.measurements import label
    from scipy.ndimage.morphology import generate_binary_structure
    from skimage.measure import regionprops

    signal = array > level
    s = generate_binary_structure(2, 2)
    labeled_array, num_features = label(signal, structure=s)
    props = regionprops(labeled_array, intensity_image=array)

    max_prop = props[0]
    for prop in props[1:]:
        if prop.max_intensity > max_prop.max_intensity:
            max_prop = prop

    bbox = max_prop.bbox
    blc = (int(bbox[1] - delta), int(bbox[0] - delta))
    trc = (int(bbox[3] + delta), int(bbox[2] + delta))
    return blc, trc


def gaussian(height, x0, y0, bmaj, e=1.0, bpa=0.0):
    """
    Returns a gaussian function with the given parameters.

    :example:
    create grid:
        x, y = np.meshgrid(x, y)
        imshow(gaussian(x, y))

    """
    import math
    # This brings PA to VLBI-convention (- = from North counter-clockwise)
    bpa = -bpa
    # Now bmaj is sigma
    bmaj = bmaj / (2. * np.sqrt(2. * np.log(2)))
    bmin = e * bmaj
    a = math.cos(bpa) ** 2. / (2. * bmaj ** 2.) + \
        math.sin(bpa) ** 2. / (2. * bmin ** 2.)
    b = math.sin(2. * bpa) / (2. * bmaj ** 2.) - \
        math.sin(2. * bpa) / (2. * bmin ** 2.)
    c = math.sin(bpa) ** 2. / (2. * bmaj ** 2.) + \
        math.cos(bpa) ** 2. / (2. * bmin ** 2.)
    return lambda x, y: height * np.exp(-(a * (x - x0) ** 2 +
                                          b * (x - x0) * (y - y0) +
                                          c * (y - y0) ** 2))


def fit_simulations_in_image_plane(sim_fname, simulations_params,
                                   core=None):
    from astropy.modeling import models, fitting
    imsize = simulations_params[u'image'][u'number_of_pixels']
    mas_in_pix = simulations_params[u'image'][u'pixel_size_mas']
    y = np.arange(-imsize/2, imsize/2, dtype=float)
    x = np.arange(-imsize/2, imsize/2, dtype=float)
    x *= mas_in_pix
    y *= mas_in_pix
    xx, yy = np.meshgrid(x, y)
    image = np.loadtxt(sim_fname)
    image[image < 0] = 0
    image[image > 10.0] = 0

    p = [1.0, 0.0, 0.0, 0.1, 0.1, 0.0]
    if core is not None:
        p = [core.p[0], core.p[1], core.p[2], core.p[3], core.p[3], 0.0]

    def tiedfunc(g_init):
        y_stddev = g_init.x_stddev
        return y_stddev

    if len(core) == 4:
        tied = {'y_stddev': tiedfunc}
        g_init = models.Gaussian2D(amplitude=p[0], x_mean=p[1], y_mean=p[2],
                                   x_stddev=p[3], y_stddev=p[3], theta=None,
                                   tied=tied)
    elif len(core) == 6:
        g_init = models.Gaussian2D(amplitude=p[0], x_mean=p[1], y_mean=p[2],
                                   x_stddev=p[3], y_stddev=p[3], theta=None)
    else:
        raise Exception

    fit_g = fitting.LevMarLSQFitter()
    return fit_g(g_init, xx, yy, image)


def plot_fitted_model(simulated_uv_fits, comps, savefig=None):
    uvdata = UVData(simulated_uv_fits)
    fig = uvdata.uvplot()
    model = Model(stokes="I")
    model.add_components(*comps)
    uvdata.substitute([model])
    fig = uvdata.uvplot(color='r', fig=fig)
    if savefig is not None:
        fig.savefig(savefig, dpi=300)
    plt.close()
    return fig


def parse_history(history):
    with open(history, "r") as fo:
        data = json.load(fo)

    keys = sorted(data.keys())
    X = list()
    y = list()
    for key in keys:
        X.append([data[key][u'parameters'][u'b'],
                  data[key][u'parameters'][u'n'],
                  data[key][u'parameters'][u'los']])
        y.append(data[key][u'results'][u'dr_obs'] -
                 data[key][u'results'][u'dr_true'])

    X = np.atleast_2d(X)
    y = np.array(y)


def parse_history_mf(history, freqs=(8, 12, 15)):
    """
    Parse history json-file of simulations.

    :param history:
        Location of json history file.
    :param freqs: (optional)
        Iterable of simulation frequencies (from low to high).
    :return:
        Two 2D numpy array. First one with ``b, n, los`` for each simulation and
        second one with ``bias_8_15, bias_12_15, bias_8_15_frac,
        bias_12_15_frac, k_obs, k_true, bias_k`` for each simulation.
    """
    with open(history, "r") as fo:
        data = json.load(fo)

    def key_func(simulation_number, freq):
        return simulation_number+'_{}'.format(freq)

    keys = sorted(data.keys())
    simulation_numbers = sorted(list(set([key.split('_')[0] for key in keys])))

    X = list()
    y = list()
    for simulation_number in simulation_numbers:
        key = key_func(simulation_number, freqs[0])
        X.append([data[key][u'parameters'][u'b'],
                  data[key][u'parameters'][u'n'],
                  data[key][u'parameters'][u'los']])

        # For each mf simulations find ``bias_8_15 = dr_8_15_obs -
        # dr_8_15_true`` and ``bias_8_12 = dr_8_12_obs - dr_8_12_true`` and
        # their relative values.
        dr_8_15_obs = data[key_func(simulation_number, freqs[0])][u'results'][u'dr_obs'] -\
                      data[key_func(simulation_number, freqs[-1])][u'results'][u'dr_obs']
        dr_8_15_true = data[key_func(simulation_number, freqs[0])][u'results'][u'dr_true'] - \
                       data[key_func(simulation_number, freqs[-1])][u'results'][u'dr_true']

        dr_12_15_obs = data[key_func(simulation_number, freqs[1])][u'results'][u'dr_obs'] - \
                      data[key_func(simulation_number, freqs[2])][u'results'][u'dr_obs']
        dr_12_15_true = data[key_func(simulation_number, freqs[1])][u'results'][u'dr_true'] - \
                       data[key_func(simulation_number, freqs[2])][u'results'][u'dr_true']

        bias_8_15 = dr_8_15_obs - dr_8_15_true
        bias_12_15 = dr_12_15_obs - dr_12_15_true
        bias_8_15_frac = bias_8_15 / dr_8_15_true
        bias_12_15_frac = bias_12_15 / dr_12_15_true

        # Calculate ``k``
        drs_obs = [dr_8_15_obs, dr_12_15_obs, 0.0]
        drs_true = [dr_8_15_true, dr_12_15_true, 0.0]
        res_obs = curve_fit(shift_model_dr, freqs, drs_obs, p0=[1.0, 1.0, -0.1],
                            bounds=((0.0, 0.0, -5.0), (10.0, 5.0, 0.0)),
                            method='trf')
        res_true = curve_fit(shift_model_dr, freqs, drs_true,
                             p0=[1.0, 1.0, -0.1], method='trf',
                             bounds=((0.0, 0.0, -5.0), (10.0, 5.0, 0.0)))

        y.append([bias_8_15, bias_12_15, bias_8_15_frac, bias_12_15_frac,
                  res_obs[0][1], res_true[0][1], res_obs[0][1]-res_true[0][1]])

    return np.atleast_2d(X), np.atleast_2d(y)


def parse_history_mf_(history, freqs=(1.665, 8.1, 12.1, 15.4)):
    """
    Parse history json-file of simulations.

    :param history:
        Location of json history file.
    :param freqs: (optional)
        Iterable of simulation frequencies (from low to high).
    :return:
        Two 2D numpy array. First one with ``b, n, los`` for each simulation and
        second one with ``k_obs, k_true, bias_k`` for each simulation.
    """
    with open(history, "r") as fo:
        data = json.load(fo)

    def key_func(simulation_number, freq):
        return simulation_number+'_{}'.format(freq)

    keys = sorted(data.keys())
    simulation_numbers = sorted(list(set([key.split('_')[0] for key in keys])))

    X = list()
    y = list()
    for simulation_number in simulation_numbers:
        key = key_func(simulation_number, freqs[0])
        X.append([data[key][u'parameters'][u'b'],
                  data[key][u'parameters'][u'n'],
                  data[key][u'parameters'][u'los']])

        dr_1_obs = data[key_func(simulation_number, freqs[0])][u'results'][u'dr_obs']
        dr_8_obs = data[key_func(simulation_number, freqs[1])][u'results'][u'dr_obs']
        dr_12_obs = data[key_func(simulation_number, freqs[2])][u'results'][u'dr_obs']
        dr_15_obs = data[key_func(simulation_number, freqs[3])][u'results'][u'dr_obs']

        dr_1_true = data[key_func(simulation_number, freqs[0])][u'results'][u'dr_true']
        dr_8_true = data[key_func(simulation_number, freqs[1])][u'results'][u'dr_true']
        dr_12_true = data[key_func(simulation_number, freqs[2])][u'results'][u'dr_true']
        dr_15_true = data[key_func(simulation_number, freqs[3])][u'results'][u'dr_true']


        # Calculate ``k``
        drs_obs = [dr_1_obs, dr_8_obs, dr_12_obs, dr_15_obs]
        drs_true = [dr_1_true, dr_8_true, dr_12_true, dr_15_true]
        res_obs = curve_fit(shift_model, freqs, drs_obs, p0=[1.0, 1.0],
                            bounds=((0.0, 0.0), (10.0, 5.0)),
                            method='trf')
        res_true = curve_fit(shift_model, freqs, drs_true,
                             p0=[1.0, 1.0], method='trf',
                             bounds=((0.0, 0.0), (10.0, 5.0)))

        y.append([res_obs[0][1], res_true[0][1], res_obs[0][1]-res_true[0][1]])

    return np.atleast_2d(X), np.atleast_2d(y)


def parse_history_mf__(history, freqs=(1.665, 8.1, 12.1, 15.4)):
    """
    Parse history json-file of simulations.

    :param history:
        Location of json history file.
    :param freqs: (optional)
        Iterable of simulation frequencies (from low to high).
    :return:
        Two 2D numpy array. First one with ``b, n, los`` for each simulation and
        second one with ``k_obs, k_true, bias_k`` for each simulation.
    """
    with open(history, "r") as fo:
        data = json.load(fo)

    def key_func(simulation_number, freq):
        return simulation_number+'_{}'.format(freq)

    keys = sorted(data.keys())
    simulation_numbers = sorted(list(set([key.split('_')[0] for key in keys])))

    X = list()
    y = list()
    fluxes = list()
    for simulation_number in simulation_numbers:
        key = key_func(simulation_number, freqs[0])
        X.append([data[key][u'parameters'][u'b'],
                  data[key][u'parameters'][u'n'],
                  data[key][u'parameters'][u'los']])

        dr_1_obs = data[key_func(simulation_number, freqs[0])][u'results'][u'dr_obs']
        dr_8_obs = data[key_func(simulation_number, freqs[1])][u'results'][u'dr_obs']
        dr_12_obs = data[key_func(simulation_number, freqs[2])][u'results'][u'dr_obs']
        dr_15_obs = data[key_func(simulation_number, freqs[3])][u'results'][u'dr_obs']

        dr_1_true = data[key_func(simulation_number, freqs[0])][u'results'][u'dr_true']
        dr_8_true = data[key_func(simulation_number, freqs[1])][u'results'][u'dr_true']
        dr_12_true = data[key_func(simulation_number, freqs[2])][u'results'][u'dr_true']
        dr_15_true = data[key_func(simulation_number, freqs[3])][u'results'][u'dr_true']


        # Calculate ``k``
        drs_obs = [1 - dr_1_obs, 1 - dr_8_obs, 1 - dr_12_obs, 1 - dr_15_obs]
        drs_true = [1 - dr_1_true, 1 - dr_8_true, 1 - dr_12_true, 1 - dr_15_true]
        res_obs = curve_fit(shift_model_dr, freqs, drs_obs, p0=[-1.0, 1.0, 1],
                            bounds=((-5.0, 0.0, 0.0), (0.0, 5.0, 5.0)),
                            method='trf')
        res_true = curve_fit(shift_model_dr, freqs, drs_true,
                             p0=[-1.0, 1.0, 1], method='trf',
                             bounds=((-5.0, 0.0, 0.0), (0.0, 5.0, 5.0)))

        y.append([res_obs[0][1], res_true[0][1], res_obs[0][1]-res_true[0][1]])

    return np.atleast_2d(X), np.atleast_2d(y)


def find_parameters_with_flux_between(json_history, flux_min, flux_max,
                                      freq=None):
    """
    Find simulation parameters that results in flux higher than specified.

    :param json_history:
        Location of json-format history file.
    :param flux_min:
        Minimum value of total flux.
    :param flux_max:
        Maximum value of total flux.
    :param freq: (optional)
        Frequency to consider. If ``None`` then consider all frequencies.
        (default: ``None``)
    :return:
        2D numpy array of simulation parameters [``b``, ``n``, ``los``].
    """
    assert flux_min < flux_max
    if freq is not None:
        assert str(freq) in ("15.4", "12.1", "8.1", "1.665")
    with open(json_history, "r") as fo:
        data = json.load(fo)
    params = list()
    keys = data.keys()
    if freq is not None:
        keys = select_keys_with_given_frequency(keys, freq)
    for key in keys:
        if flux_min < data[key][u'results'][u'flux'] < flux_max:
            params.append([data[key][u'parameters'][u'b'],
                           data[key][u'parameters'][u'n'],
                           data[key][u'parameters'][u'los']])
    return np.atleast_2d(params)


def find_params_for_given_keys(json_history, keys):
    with open(json_history, "r") as fo:
        data = json.load(fo)
    params = list()
    for key in keys:
        if key not in data.keys():
            raise Exception("No key {} in {} keys!".format(key, json_history))
        params.append([data[key][u'parameters'][u'b'],
                       data[key][u'parameters'][u'n'],
                       data[key][u'parameters'][u'los']])
    return np.atleast_2d(params)


def find_keys_with_flux_between(json_history, flux_min, flux_max):
    """
    Find keys of results in flux higher than specified.

    :param json_history:
        Location of json-format history file.
    :param flux_min:
        Minimum value of total flux.
    :param flux_max:
        Maximum value of total flux.
    :return:
        Iterable of keys.
    """
    assert flux_min < flux_max
    with open(json_history, "r") as fo:
        data = json.load(fo)
    keys = list()
    for key in data.keys():
        if flux_min < data[key][u'results'][u'flux'] < flux_max:
            keys.append(key)
    return keys


def select_keys_with_given_frequency(keys, freq):
    """
    Among iterable of keys select only those corresponding to given frequency.

    :param keys:
        Iterable of keys.
    :param freq:
        Frequency (``15.4``, ``12.1``, ``8.1`` or ``1.665``).
    :return:
        Iterable of keys.
    """
    new_keys = list()
    for key in keys:
        if str(freq) in key:
            new_keys.append(key)
    return new_keys


def find_tb_for_given_keys(json_history, keys):
    """
    Retrieve brightness temperatures for given keys.

    :param json_history:
        Location of json-format history file.
    :param keys:
        Iterable of keys.
    :return:
        2D numpy array with brightness temperatures - for simulation image
        pixels and for difmap models.
    """

    with open(json_history, "r") as fo:
        data = json.load(fo)
    tb = list()
    for key in keys:
        tb.append([data[key][u'results'][u'tb_pix'],
                   data[key][u'results'][u'tb_difmap']])
    return np.atleast_2d(tb)


def find_flux_for_given_keys(json_history, keys):
    with open(json_history, "r") as fo:
        data = json.load(fo)
    fluxes = list()
    for key in keys:
        if key not in data.keys():
            raise Exception("No key {} in {} keys!".format(key, json_history))
        fluxes.append(data[key][u'results'][u'flux'])
    return np.array(fluxes)


def strip_frequency_from_keys(keys):
    """
    Convert full keys (with frequency) to two-decimal keys, e.g.
    ``32_15.4`` -> ``32``.

    :param keys:
        Iterable of keys with frequencies.
    :return:
        Iterable of keys without frequencies.
    """
    return list(set([key.split('_')[0] for key in keys]))


def find_shifts_and_fluxes_for_given_keys_and_freqs(json_history, keys,
                                                    freq_low, freq_high):
    """
    Get shifts between given frequencies for given keys.

    :param json_history:
        Location of json-format history file.
    :param keys:
        Iterable of keys that consists of only 2 digit integer, e.g. ``23`` or
        ``keys = [str(i).zfill(2) for i in range(1, 34)]``
    :param freq_low:
        Lowest frequency to consider. (``15.4``, ``12.1``, ``8.1`` or ``1.665``)
    :param freq_high:
        Highest frequency to consider. (``15.4``, ``12.1``, ``8.1`` or
        ``1.665``)
    :return:
        2D numpy array with shifts between ``freq_low`` and ``freq_high`` for
        given keys - observed value, true, value and their difference.
    """
    assert float(freq_low) < float(freq_high)
    with open(json_history, "r") as fo:
        data = json.load(fo)
    shifts = list()
    fluxes = list()
    for key in keys:
        key_low = u'{}_{}'.format(key, freq_low)
        key_high = u'{}_{}'.format(key, freq_high)

        dr_obs_low = data[key_low][u'results'][u'dr_obs']
        dr_obs_high = data[key_high][u'results'][u'dr_obs']
        shift_obs = dr_obs_low - dr_obs_high

        dr_true_low = data[key_low][u'results'][u'dr_true']
        dr_true_high = data[key_high][u'results'][u'dr_true']
        shift_true = dr_true_low - dr_true_high

        shifts.append([shift_obs, shift_true, shift_obs-shift_true])
        fluxes.append(data[key_high][u'results'][u'flux'])
    return np.atleast_2d(shifts), np.array(fluxes)


def find_shifts_and_los_for_given_keys_and_freqs(json_history, keys,
                                                 freq_low, freq_high):
    """
    Get shifts between given frequencies for given keys.

    :param json_history:
        Location of json-format history file.
    :param keys:
        Iterable of keys that consists of only 2 digit integer, e.g. ``23`` or
        ``keys = [str(i).zfill(2) for i in range(1, 34)]``
    :param freq_low:
        Lowest frequency to consider. (``15.4``, ``12.1``, ``8.1`` or ``1.665``)
    :param freq_high:
        Highest frequency to consider. (``15.4``, ``12.1``, ``8.1`` or
        ``1.665``)
    :return:
        2D numpy array with shifts between ``freq_low`` and ``freq_high`` for
        given keys - observed value, true, value and their difference.
    """
    assert float(freq_low) < float(freq_high)
    with open(json_history, "r") as fo:
        data = json.load(fo)
    shifts = list()
    los = list()
    for key in keys:
        key_low = u'{}_{}'.format(key, freq_low)
        key_high = u'{}_{}'.format(key, freq_high)

        dr_obs_low = data[key_low][u'results'][u'dr_obs']
        dr_obs_high = data[key_high][u'results'][u'dr_obs']
        shift_obs = dr_obs_low - dr_obs_high

        dr_true_low = data[key_low][u'results'][u'dr_true']
        dr_true_high = data[key_high][u'results'][u'dr_true']
        shift_true = dr_true_low - dr_true_high

        shifts.append([shift_obs, shift_true, shift_obs-shift_true])
        los.append(data[key_high][u'parameters'][u'los'])
    return np.atleast_2d(shifts), np.array(los)


def join_json_histories(result_json_history, overwrite=False, *json_histories):
    """
    Join several json-format history files in one.

    :param result_json_history:
        Location of resulted json history file.
    :param overwrite: (optional)
        Boolean - overwrite repeating keys in resulting history? (default:
        ``False``)
    :param json_histories:
        Iterable of json-format history files to concatenate.
    :return:
        Json-format history file with all info from others included.
    """
    histories = nested_dict()
    for json_history in json_histories:
        with open(json_history, 'r') as fo:
            history = json.load(fo)
        keys = history.keys()
        for key in keys:
            if key in histories.keys():
                if not overwrite:
                    print("Resulted json already has key {}! Skipping".format(key))
                    continue
            histories[key] = history[key]

    with open(result_json_history, 'w') as fo:
        json.dump(histories, fo)


def fit_parameters_to_flux(json_history, flux_min, flux_max, freq):
    """
    :param json_history:
    :param flux_min:
    :param flux_max:
    :param freq:
    :return:

    :note:
        Get parameters that results in fluxes between 0.3 and 0.5:
        ``X[np.logical_and(y>0.3, y<0.5)]``
    """
    from mpl_toolkits.mplot3d import Axes3D
    keys = find_keys_with_flux_between(json_history, flux_min, flux_max)
    keys = select_keys_with_given_frequency(keys, freq)
    parameters = find_params_for_given_keys(json_history, keys)
    fluxes = find_flux_for_given_keys(json_history, keys)
    from sklearn import gaussian_process
    from sklearn.preprocessing import MinMaxScaler
    mms = MinMaxScaler()
    X = mms.fit_transform(np.log10(parameters))
    y = fluxes.copy()
    y -= np.median(fluxes)
    gp = gaussian_process.GaussianProcess(thetaL=(0.1, 0.1, 0.1),
                                          thetaU=(0.5, 0.5, 0.5),
                                          theta0=(0.25, 0.25, 0.25),
                                          nugget=0.05**2,
                                          storage_mode='full')
    gp.fit(X, y)
    XX, YY, ZZ = np.meshgrid(np.linspace(0, 1, 50),
                             np.linspace(0, 1, 50),
                             np.linspace(0, 1, 50))
    X_ = np.vstack((XX.ravel(), YY.ravel(), ZZ.ravel())).T
    y_ = gp.predict(X_)
    y_ = y_.reshape(XX.shape)
    y_ += np.median(fluxes)
    X = np.atleast_2d([np.linspace(0, 1, 50),
                       np.linspace(0, 1, 50),
                       np.linspace(0, 1, 50)]).T
    X = 10**mms.inverse_transform(X)
    XX, YY, ZZ = np.meshgrid(X[:, 0], X[:, 1], X[:, 2])
    plt.matshow(y_[:,:,0])
    plt.contour(XX[:,:,0], YY[:,:,0], y_[:,:,0])
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(XX[:,:,0], YY[:,:,0], y_[:,:,0], c='r', marker='o')

    return X_, y_


def plot_stripe_for_given_key(simulated_image_low, simulated_image_high,
                              json_history, key_low, key_high,
                              difmap_model_low, difmap_model_high,
                              delta=100, savefig=None):
    """
     Plot 1D stripe along jet axis.
     """
    fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
    text = ['{} GHz'.format(key_low.split('_')[1]),
            '{} GHz'.format(key_high.split('_')[1])]
    positions = [(0.95, 0.8), (0.95, 0.10)]
    for i, simulated_image, key, difmap_model in zip((0, 1),
                                           (simulated_image_low, simulated_image_high),
                                           (key_low, key_high),
                                           (difmap_model_low, difmap_model_high)):
        difmap_dir, difmap_fn = os.path.split(difmap_model)
        comps = import_difmap_model(difmap_fn, difmap_dir)
        core = comps[0]

        image = np.loadtxt(simulated_image)
        image[image < 0] = 0
        image[image > 10.0] = 0
        with open(json_history, 'r') as fo:
            simulations_params = json.load(fo)
        imsize = simulations_params[key][u'image'][u'number_of_pixels']
        mas_in_pix = simulations_params[key][u'image'][u'pixel_size_mas']

        y = np.arange(-imsize/2, imsize/2, dtype=float)
        x = np.arange(-imsize/2, imsize/2, dtype=float)
        x *= mas_in_pix
        y *= mas_in_pix
        xx, yy = np.meshgrid(x, y)

        # Fit simulated image in image plane
        p = [1.0, 0.0, 0.0, 0.1, 0.1, 0.0]
        if core is not None:
            p = [core.p[0], core.p[1], core.p[2], core.p[3], core.p[3], 0.0]

        def tiedfunc(g_init):
            y_stddev = g_init.x_stddev
            return y_stddev

        if len(core) == 4:
            tied = {'y_stddev': tiedfunc}
            g_init = models.Gaussian2D(amplitude=p[0], x_mean=p[1], y_mean=p[2],
                                       x_stddev=p[3], y_stddev=p[3], theta=None,
                                       tied=tied)
        elif len(core) == 6:
            g_init = models.Gaussian2D(amplitude=p[0], x_mean=p[1], y_mean=p[2],
                                       x_stddev=p[3], y_stddev=p[3], theta=None)
        else:
            raise Exception

        fit_g = fitting.LevMarLSQFitter()
        g = fit_g(g_init, xx, yy, image)
        g_image = g(xx, yy)

        if len(core) == 6:
            e = core.p[4]
        else:
            e = 1.
        # put 2 instead of 4 before sqrt and remove 0.5 before p[3]
        amp = core.p[0] / (2. * np.pi * e * (core.p[3] / (2. * np.sqrt(2. * np.log(2)) * mas_in_pix)) ** 2)
        p = [amp, core.p[1], core.p[2], core.p[3]] + list(core.p[4:])
        gaus = gaussian(*p)
        gaus_image = gaus(xx, yy)

        x = np.arange(-20, imsize / 2 - delta) * mas_in_pix

        colors = plt.rcParams['axes.color_cycle']
        # Plot simulation results
        axes[i].plot(x, image[imsize / 2, imsize / 2 - 20:-delta], lw=2, label="Simulations",
                  color=colors[0])
        axes[i].axvline(x[np.argmax(image[imsize / 2, imsize / 2 - 20:-delta])], lw=1,
                        color=colors[0])
        axes[i].axvline(core.p[1], lw=2, color=colors[1])
        # Plot difmap model
        axes[i].plot(x, gaus_image[imsize / 2, imsize / 2 - 20:-delta], label="Difmap fit",
                  color=colors[1])
        ratio = image[imsize / 2, imsize / 2 - 20:-delta].max()/gaus_image[imsize / 2, imsize / 2 - 20:-delta].max()
        axes[i].plot(x, ratio*gaus_image[imsize / 2, imsize / 2 - 20:-delta], ':',
                  color=colors[1], lw=1)
        # Optionally plot model fitted in image plane
        axes[i].plot(x, g_image[imsize / 2, imsize / 2 - 20:-delta], label="Sim. image fit",
                  color=colors[2])
        axes[i].axvline(core.p[1], lw=2, color=colors[1])
        axes[i].set_ylabel("Flux, [Jy/pixel]", fontsize=14)
        axes[i].text(positions[i][0], positions[i][1], text[i],
                     transform=axes[i].transAxes, fontsize=14,
                     horizontalalignment='right')
        plt.legend(loc="best", fontsize=14)
    axes[i].set_xlabel("Distance along jet, [mas]", fontsize=14)
    fig.tight_layout()
    if savefig:
        plt.savefig(savefig, bbox_inches='tight', dpi=300)

    return fig


def comoving_transverse_distance(z, H_0=73.0, omega_M=0.3, omega_V=0.7,
                                 format="pc"):
    """
    Given redshift ``z``, Hubble constant ``H_0`` [km/s/Mpc] and
    density parameters ``omega_M`` and ``omega_V``, returns comoving transverse
    distance (see arXiv:astro-ph/9905116v4 formula 14). Angular diameter
    distance is factor (1 + z) lower and luminosity distance is the same factor
    higher.

    """
    from scipy.integrate import quad
    fmt_dict = {"cm": 9.26 * 10.0 ** 27.0, "pc": 3. * 10 ** 9, "Mpc": 3000.0}

    result = (H_0 / 100.0) ** (-1.) * quad(lambda x: (omega_M * (1. + x ** 3) +
                                                      omega_V) ** (-0.5),
                                           0, z)[0]
    try:
        return fmt_dict[format] * result
    except KeyError:
        raise Exception('Format  \"pc\", \"cm\" or \"Mpc\"')


def pc_to_mas(z):
    """
    Return scale factor that convert from parsecs to milliarcseconds .

    """
    # Angular distance in pc
    d_a = comoving_transverse_distance(z, format='pc') / (1. + z)
    # Angle in radians
    angle = 1. / d_a
    return rad_to_mas * angle


def mas_to_pc(z):
    """
    Return scale factor that convert from milliarcseconds to parsecs.

    """
    # Angular distance in pc
    d_a = comoving_transverse_distance(z, format='pc') / (1. + z)
    return mas_to_rad * d_a


def distance_from_SMBH(dr_mas, los_rad, z):
    """
    Returns distance from SMBH [pc] at given projected distance [mas].
    :param dr_mas:
    :param los_rad:
    :param z:
    :return:
    """
    return dr_mas*mas_to_pc(z)/np.sin(los_rad)


def b_field(b1, r_pc, m=1):
    """
    Value of B field at ``r_pc`` distance [pc]f rom apex.

    :param b1:
    :param r_pc:
    :param m:
    :return:
    """
    return b1 * (r_pc)**(-m)


def _find_shifts_from_true_images(freq_true_images_dict, imsize,
                                  pixelsizes_dict):
    """
  Find shift using true images.

  :param freq_true_images_dict:
      Dictionary with keys - frequency [GHz] and values - paths to images of
       true intensity distributions with core.
  :return:
      Dictionary with keys - name of the value used to measure shift and
      values - measured values.
  """
    drs = list()
    freqs = list()
    result_dict = collections.OrderedDict()
    for freq, true_image_path in freq_true_images_dict.items():
        true_image = np.loadtxt(true_image_path)
        # Distance from SMBH (0,0) in mas
        dr = (np.unravel_index(true_image.argmax(),
                               true_image.shape)[1]-imsize/2)*pixelsizes_dict[freq]
        drs.append(dr)
        freqs.append(freq)
    res = curve_fit(shift_model, freqs, drs, p0=[1.0, 1.0])
    result_dict['dr'] = (res[0], drs)
    return result_dict


def _make_plot_i_tau(image_i, image_tau, angle, los_angle, imsize, out_dir=None,
                     name="", title=None):
    if out_dir is None:
        out_dir = os.getcwd()
    image_i = np.loadtxt(image_i)
    image_tau = np.loadtxt(image_tau)
    size = int(0.5*imsize*(angle/los_angle))
    islice = slice(size, imsize-size, None)
    import matplotlib.pyplot as plt
    plt.matshow(image_i[islice])
    plt.contour(np.log10(image_tau[islice]), [-2, -1, -0.5, 0, 0.5, 1, 2],
                cmap='tab10')
    plt.colorbar()
    plt.axis('off')
    if title is not None:
        plt.title(title)
    plt.savefig(os.path.join(out_dir, "image_I_tau_{}.png".format(name)))
    plt.close()


if __name__ == '__main__':
    # Run simulation
    main_dir = '/home/ilya/github/bck/jetshow'
    cfg_file = os.path.join(main_dir, 'config.json')
    with open(cfg_file, "r") as jsonFile:
        simulation_params = json.load(jsonFile)
    path_to_executable = os.path.join(main_dir, 'cmake-build-debug', 'jetshow')
    simulation_params = run_simulations(cfg_file, path_to_executable)

    # Create artificial data set with BK core and jet component
    exe_dir, exe = os.path.split(path_to_executable)
    initial_dfm_model = os.path.join(main_dir, 'initial_eg.mdl')
    # initial_dfm_model = os.path.join(main_dir, 'initial_cg.mdl')
    uv_fits_template = '/home/ilya/github/vlbi_errors/vlbi_errors/'
    uv_fits_template = os.path.join(uv_fits_template,
                                    '0235+164.u.2006_06_15.uvf')
    modelfit_simulation_result(exe_dir, initial_dfm_model, noise_factor=1.0,
                               out_dfm_model_fn="bk_e.mdl", out_dir=exe_dir,
                               params=simulation_params,
                               uv_fits_save_fname="bk.fits",
                               uv_fits_template=uv_fits_template)

    # Find measured and true distance of core to jet component
    # dr_obs = find_core_separation_from_jet_using_difmap_model(os.path.join(exe_dir, "bk.mdl"))
    # Here ``1.`` is because jet component put at 1 mas distance from phase
    # center of the originaldata set.
    # dr_true = 1. - find_core_separation_from_center_using_simulations(os.path.join(exe_dir, "map_i.txt"),
    #                                                                   simulation_params)

    # Plot map with components superimposed
    from spydiff import clean_difmap, import_difmap_model
    from image import plot as iplot
    from image import find_bbox
    from from_fits import create_clean_image_from_fits_file
    from image_ops import rms_image

    path_to_script = '/home/ilya/github/vlbi_errors/difmap/final_clean_nw'
    clean_difmap('bk.fits', 'bk_cc.fits', 'I', (1024, 0.1), path=exe_dir,
                 path_to_script=path_to_script, show_difmap_output=True,
                 outpath=exe_dir)

    ccimage = create_clean_image_from_fits_file(os.path.join(exe_dir,
                                                             'bk_cc.fits'))
    beam = ccimage.beam
    rms = rms_image(ccimage)
    blc, trc = find_bbox(ccimage.image, rms, 10)
    comps = import_difmap_model('bk_e.mdl', exe_dir)
    iplot(ccimage.image, x=ccimage.x, y=ccimage.y, min_abs_level=3*rms,
          beam=beam, show_beam=True, blc=blc, trc=trc, components=comps)

    g = fit_simulations_in_image_plane(os.path.join(exe_dir, "map_i.txt"),
                                       simulation_params, core=comps[0])

    plot_stripe(os.path.join(exe_dir, "map_i.txt"),
                os.path.join(exe_dir, "bk_e.mdl"), simulation_params, g=g)

    # plot_simulations_3d(os.path.join(exe_dir, "map_i.txt"),
    #                     simulation_params, core=comps[0], each=1, delta=130)

    plot_simulations_2d(os.path.join(exe_dir, "map_i.txt"),
                        simulation_params, core=comps[0], each=1,
                        side_delta=110, g=g)

    # initial_dfm_model = os.path.join(main_dir, 'initial_eg.mdl')
    # executable_dir, _ = os.path.split(path_to_executable)
    # uv_fits_template = '/home/ilya/github/vlbi_errors/vlbi_errors/'
    # uv_fits_template_dict = collections.OrderedDict()
    # uv_fits_template_dict[1.7] = os.path.join(uv_fits_template,
    #                                           '0235+164.18cm.2010_06_23.uvf')
    # uv_fits_template_dict[8.1] = os.path.join(uv_fits_template,
    #                                           '0235+164.x.2006_06_15.uvf')
    # uv_fits_template_dict[8.4] = os.path.join(uv_fits_template,
    #                                           '0235+164.y.2006_06_15.uvf')
    # uv_fits_template_dict[12.1] = os.path.join(uv_fits_template,
    #                                            '0235+164.j.2006_06_15.uvf')
    # uv_fits_template_dict[15.4] = os.path.join(uv_fits_template,
    #                                           '0235+164.u.2006_06_15.uvf')
    # update_params_dict = nested_dict()
    # pixsizes_dict = nested_dict()
    # imsizes_dict = nested_dict()
    # freqs = [1.7, 8.1, 12.1, 15.4]
    # # los_angles = [0.035, 0.0525, 0.07, 0.0875]
    # los_angles = [0.0875]
    # # angles = [0.0175, 0.035]
    # angles = [0.0175]
    # bs = [0.1, 1, 10]
    # ns = [50, 500, 5000]
    #
    # # /40
    # pixsizes_dict[0.1][50] = 0.0001
    # imsizes_dict[0.1][50] = 1000
    # # /20
    # pixsizes_dict[0.1][500] = 0.0002
    # imsizes_dict[0.1][500] = 1000
    # # /10
    # pixsizes_dict[0.1][5000] = 0.0004
    # imsizes_dict[0.1][5000] = 1000
    #
    # # /10
    # pixsizes_dict[1][50] = 0.0004
    # imsizes_dict[1][50] = 1000
    # # /4
    # pixsizes_dict[1][500] = 0.001
    # imsizes_dict[1][500] = 1000
    # # /2
    # pixsizes_dict[1][5000] = 0.002
    # imsizes_dict[1][5000] = 1000
    #
    # # /2
    # pixsizes_dict[10][50] = 0.002
    # imsizes_dict[10][50] = 1000
    # # /1
    # pixsizes_dict[10][500] = 0.004
    # imsizes_dict[10][500] = 1000
    # # *2
    # pixsizes_dict[10][5000] = 0.008
    # imsizes_dict[10][5000] = 1000
    #
    #
    #
    # for angle in angles:
    #     for los_angle in los_angles:
    #         for b in bs:
    #             for n in ns:
    #                 if angle == 0.0175:
    #                     if los_angle == 0.035:
    #                         if b == 0.1:
    #                             if n == 50 or n == 500 or n == 5000:
    #                                 continue
    #                         if b == 1:
    #                             if n == 50:
    #                                 continue
    #                 # Here cycle for different values of the parameters
    #                 # los_angle = 0.035
    #                 # angle = 0.0175
    #                 # b = 1.0
    #                 # n = 500.0
    #                 # number_of_pixels = 1000
    #                 number_of_pixels = imsizes_dict[b][n]
    #                 # pixel_size_mas = 0.004
    #                 pixel_size_mas = pixsizes_dict[b][n]
    #
    #                 freq_pixel_size_dict = collections.OrderedDict()
    #                 for freq in freqs:
    #                     if freq == 1.7:
    #                         pixel_size_mas_ = 10*pixel_size_mas
    #                     else:
    #                         pixel_size_mas_ = pixel_size_mas
    #                     freq_pixel_size_dict[freq] = pixel_size_mas_
    #
    #                 print("Pixel sizes : ", freq_pixel_size_dict)
    #                 freq_difmap_models_dict = collections.OrderedDict()
    #                 freq_true_images_dict = collections.OrderedDict()
    #                 out_dir = os.path.join(main_dir,
    #                                        'prod_imsize{}_pix{}_los{}_angle{}_b{}_n{}'.format(number_of_pixels,
    #                                                                     pixel_size_mas, los_angle, angle, b, n))
    #                 if not os.path.exists(out_dir):
    #                     os.mkdir(out_dir)
    #
    #                 # Simulate image for each frequency
    #                 for freq in freqs:
    #                     update_params_dict[u'image'][u'number_of_pixels'] = number_of_pixels
    #                     update_params_dict[u'image'][u'pixel_size_mas'] = freq_pixel_size_dict[freq]
    #                     update_params_dict[u'integration'][u'parameters'][u'log10_tau_min'] = -4.0
    #                     update_params_dict[u'observation'][u'frequency_ghz'] = freq
    #                     update_params_dict[u'observation'][u'los_angle'] = los_angle
    #                     update_params_dict[u'jet'][u'geometry'][u'parameters'][u'angle'] = angle
    #                     update_params_dict[u'jet'][u'bfield'][u'parameters'][u'b_1'] = b
    #                     update_params_dict[u'jet'][u'nfield'][u'parameters'][u'n_1'] = n
    #                     update_params_dict[u'output'][u'file_i'] = 'map_i_{}.txt'.format(freq)
    #                     update_params_dict[u'output'][u'file_tau'] = 'map_tau_{}.txt'.format(freq)
    #                     update_params_dict[u'output'][u'file_length'] = 'map_l_{}.txt'.format(freq)
    #                     run_simulations(cfg_file, out_dir, initial_dfm_model,
    #                                     path_to_executable, uv_fits_template_dict[freq],
    #                                     uv_fits_save_fname='bk_{}.fits'.format(freq),
    #                                     out_dfm_model_fn='bk_{}.mdl'.format(freq),
    #                                     update_params_dict=update_params_dict,
    #                                     noise_factor=0.01)
    #
    #                     # Move simulated images to from directory with executable to ``out_dir``
    #                     os.rename(os.path.join(executable_dir, 'map_i_{}.txt'.format(freq)),
    #                               os.path.join(out_dir, 'map_i_{}.txt'.format(freq)))
    #                     os.rename(os.path.join(executable_dir, 'map_tau_{}.txt'.format(freq)),
    #                               os.path.join(out_dir, 'map_tau_{}.txt'.format(freq)))
    #                     os.rename(os.path.join(executable_dir, 'map_l_{}.txt'.format(freq)),
    #                               os.path.join(out_dir, 'map_l_{}.txt'.format(freq)))
    #
    #                     make_plot_i_tau(os.path.join(out_dir, 'map_i_{}.txt'.format(freq)),
    #                                     os.path.join(out_dir, 'map_tau_{}.txt'.format(freq)),
    #                                     angle, los_angle, number_of_pixels, out_dir=out_dir,
    #                                     name="{}_GHz".format(freq),
    #                                     title="LOS={} Angle={} B={} N={} freq={}".format(los_angle, angle, b, n, freq))
    #                     freq_difmap_models_dict[freq] = os.path.join(out_dir,
    #                                                                  'bk_{}.mdl'.format(freq))
    #                     freq_true_images_dict[freq] = os.path.join(out_dir,
    #                                                                'map_i_{}.txt'.format(freq))
    #
    #                 # Calculate shifts in different ways
    #                 observed_shifts = find_shift_from_difmap_models(freq_difmap_models_dict)
    #                 true_shifts = find_shifts_from_true_images(freq_true_images_dict,
    #                                                            number_of_pixels,
    #                                                            freq_pixel_size_dict)
    #                 bias_dr = observed_shifts['dr'][0][1] - true_shifts['dr'][0][1]
    #                 bias_bmaj = observed_shifts['bmaj'][0][1] - 1.0
    #                 bias_bmin = observed_shifts['bmin'][0][1] - 1.0
    #
    #                 import matplotlib.pyplot as plt
    #                 # Plot all measured values
    #                 plt.plot(freqs, observed_shifts['dr'][1], '.k', ms=10,
    #                          label="k={0:.2f} observed dr".format(observed_shifts['dr'][0][1]))
    #                 plt.plot(freqs, true_shifts['dr'][1], '.r', ms=10,
    #                          label="k={0:.2f} true dr".format(true_shifts['dr'][0][1]))
    #                 plt.plot(freqs, observed_shifts['bmaj'][1], '.b', ms=10,
    #                          label="k={0:.2f} observed bmaj".format(observed_shifts['bmaj'][0][1]))
    #                 plt.plot(freqs, observed_shifts['bmin'][1], '.g', ms=10,
    #                          label="k={0:.2f} observed bmin".format(observed_shifts['bmin'][0][1]))
    #
    #                 freqs_grid = np.linspace(freqs[0], freqs[-1], 100)
    #                 drs_observed_fit = shift_model(freqs_grid, observed_shifts['dr'][0][0],
    #                                                observed_shifts['dr'][0][1])
    #                 bmajs_observed_fit = shift_model(freqs_grid, observed_shifts['bmaj'][0][0],
    #                                                  observed_shifts['bmaj'][0][1])
    #                 bmins_observed_fit = shift_model(freqs_grid, observed_shifts['bmin'][0][0],
    #                                                  observed_shifts['bmin'][0][1])
    #                 drs_true_fit = shift_model(freqs_grid, true_shifts['dr'][0][0],
    #                                            true_shifts['dr'][0][1])
    #                 plt.plot(freqs_grid, drs_observed_fit, 'k')
    #                 plt.plot(freqs_grid, bmajs_observed_fit, 'b')
    #                 plt.plot(freqs_grid, bmins_observed_fit, 'g')
    #                 plt.plot(freqs_grid, drs_true_fit, 'r')
    #                 plt.xlabel("Frequency, GHz")
    #                 plt.ylabel("Shift/Size, mas")
    #                 plt.legend(loc='best')
    #                 plt.title("LOS={} Angle={} B={} N={}".format(los_angle, angle, b, n))
    #                 plt.savefig(os.path.join(out_dir, 'fits.png'), bbox_inches='tight')
    #                 plt.close()
    #
    #                 print("True k (dr) = {}".format(true_shifts['dr'][0][1]))
    #                 print("Bias of observed k (dr) = {}".format(bias_dr))
    #                 print("Bias of observed k (bmaj) = {}".format(bias_bmaj))
    #                 print("Bias of observed k (bmin) = {}".format(bias_bmin))