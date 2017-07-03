import os
import sys
sys.path.insert(0, '/home/ilya/github/vlbi_errors/vlbi_errors')
import subprocess
import collections
import json
import numpy as np
from scipy.optimize import curve_fit
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
                    update_params_dict=None, noise_factor=1.0):
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
    for key, value in noise.items():
        noise[key] = noise_factor*value
    model = Model(stokes='I')
    model.add_component(icomp)
    uvdata.substitute([model])
    uvdata.noise_add(noise)
    uvdata.save(os.path.join(out_dir, uv_fits_save_fname))

    initial_dfm_dir, initial_dfm_model_fn = os.path.split(initial_dfm_model)
    # Fit model in difmap
    modelfit_difmap(uv_fits_save_fname, initial_dfm_model_fn, out_dfm_model_fn,
                    niter=300, path=out_dir,
                    mdl_path=initial_dfm_dir,
                    out_path=out_dir)


def nested_dict():
    return collections.defaultdict(nested_dict)


def shift_model(nu, a, k):
    return a * nu**(-1./k)


def find_shift_from_difmap_models(freq_difmap_models_dict):
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


def find_shifts_from_true_images(freq_true_images_dict, imsize, pixelsizes_dict):
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


def make_plot_i_tau(image_i, image_tau, angle, los_angle, imsize, out_dir=None,
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
    main_dir = '/home/ilya/github/bck/jetshow'

    cfg_file = os.path.join(main_dir, 'config.json')
    initial_dfm_model = os.path.join(main_dir, 'initial_eg.mdl')
    path_to_executable = os.path.join(main_dir, 'cmake-build-debug', 'jetshow')
    executable_dir, _ = os.path.split(path_to_executable)
    uv_fits_template = '/home/ilya/github/vlbi_errors/vlbi_errors/'
    uv_fits_template_dict = collections.OrderedDict()
    uv_fits_template_dict[1.7] = os.path.join(uv_fits_template,
                                              '0235+164.18cm.2010_06_23.uvf')
    uv_fits_template_dict[8.1] = os.path.join(uv_fits_template,
                                              '0235+164.x.2006_06_15.uvf')
    uv_fits_template_dict[8.4] = os.path.join(uv_fits_template,
                                              '0235+164.y.2006_06_15.uvf')
    uv_fits_template_dict[12.1] = os.path.join(uv_fits_template,
                                               '0235+164.j.2006_06_15.uvf')
    uv_fits_template_dict[15.4] = os.path.join(uv_fits_template,
                                              '0235+164.u.2006_06_15.uvf')
    update_params_dict = nested_dict()
    pixsizes_dict = nested_dict()
    imsizes_dict = nested_dict()
    freqs = [1.7, 8.1, 12.1, 15.4]
    # los_angles = [0.035, 0.0525, 0.07, 0.0875]
    los_angles = [0.0875]
    # angles = [0.0175, 0.035]
    angles = [0.0175]
    bs = [0.1, 1, 10]
    ns = [50, 500, 5000]

    # /40
    pixsizes_dict[0.1][50] = 0.0001
    imsizes_dict[0.1][50] = 1000
    # /20
    pixsizes_dict[0.1][500] = 0.0002
    imsizes_dict[0.1][500] = 1000
    # /10
    pixsizes_dict[0.1][5000] = 0.0004
    imsizes_dict[0.1][5000] = 1000

    # /10
    pixsizes_dict[1][50] = 0.0004
    imsizes_dict[1][50] = 1000
    # /4
    pixsizes_dict[1][500] = 0.001
    imsizes_dict[1][500] = 1000
    # /2
    pixsizes_dict[1][5000] = 0.002
    imsizes_dict[1][5000] = 1000

    # /2
    pixsizes_dict[10][50] = 0.002
    imsizes_dict[10][50] = 1000
    # /1
    pixsizes_dict[10][500] = 0.004
    imsizes_dict[10][500] = 1000
    # *2
    pixsizes_dict[10][5000] = 0.008
    imsizes_dict[10][5000] = 1000



    for angle in angles:
        for los_angle in los_angles:
            for b in bs:
                for n in ns:
                    if angle == 0.0175:
                        if los_angle == 0.035:
                            if b == 0.1:
                                if n == 50 or n == 500 or n == 5000:
                                    continue
                            if b == 1:
                                if n == 50:
                                    continue
                    # Here cycle for different values of the parameters
                    # los_angle = 0.035
                    # angle = 0.0175
                    # b = 1.0
                    # n = 500.0
                    # number_of_pixels = 1000
                    number_of_pixels = imsizes_dict[b][n]
                    # pixel_size_mas = 0.004
                    pixel_size_mas = pixsizes_dict[b][n]

                    freq_pixel_size_dict = collections.OrderedDict()
                    for freq in freqs:
                        if freq == 1.7:
                            pixel_size_mas_ = 10*pixel_size_mas
                        else:
                            pixel_size_mas_ = pixel_size_mas
                        freq_pixel_size_dict[freq] = pixel_size_mas_

                    print("Pixel sizes : ", freq_pixel_size_dict)
                    freq_difmap_models_dict = collections.OrderedDict()
                    freq_true_images_dict = collections.OrderedDict()
                    out_dir = os.path.join(main_dir,
                                           'prod_imsize{}_pix{}_los{}_angle{}_b{}_n{}'.format(number_of_pixels,
                                                                        pixel_size_mas, los_angle, angle, b, n))
                    if not os.path.exists(out_dir):
                        os.mkdir(out_dir)

                    # Simulate image for each frequency
                    for freq in freqs:
                        update_params_dict[u'image'][u'number_of_pixels'] = number_of_pixels
                        update_params_dict[u'image'][u'pixel_size_mas'] = freq_pixel_size_dict[freq]
                        update_params_dict[u'integration'][u'parameters'][u'log10_tau_min'] = -4.0
                        update_params_dict[u'observation'][u'frequency_ghz'] = freq
                        update_params_dict[u'observation'][u'los_angle'] = los_angle
                        update_params_dict[u'jet'][u'geometry'][u'parameters'][u'angle'] = angle
                        update_params_dict[u'jet'][u'bfield'][u'parameters'][u'b_1'] = b
                        update_params_dict[u'jet'][u'nfield'][u'parameters'][u'n_1'] = n
                        update_params_dict[u'output'][u'file_i'] = 'map_i_{}.txt'.format(freq)
                        update_params_dict[u'output'][u'file_tau'] = 'map_tau_{}.txt'.format(freq)
                        update_params_dict[u'output'][u'file_length'] = 'map_l_{}.txt'.format(freq)
                        run_simulations(cfg_file, out_dir, initial_dfm_model,
                                        path_to_executable, uv_fits_template_dict[freq],
                                        uv_fits_save_fname='bk_{}.fits'.format(freq),
                                        out_dfm_model_fn='bk_{}.mdl'.format(freq),
                                        update_params_dict=update_params_dict,
                                        noise_factor=0.01)

                        # Move simulated images to from directory with executable to ``out_dir``
                        os.rename(os.path.join(executable_dir, 'map_i_{}.txt'.format(freq)),
                                  os.path.join(out_dir, 'map_i_{}.txt'.format(freq)))
                        os.rename(os.path.join(executable_dir, 'map_tau_{}.txt'.format(freq)),
                                  os.path.join(out_dir, 'map_tau_{}.txt'.format(freq)))
                        os.rename(os.path.join(executable_dir, 'map_l_{}.txt'.format(freq)),
                                  os.path.join(out_dir, 'map_l_{}.txt'.format(freq)))

                        make_plot_i_tau(os.path.join(out_dir, 'map_i_{}.txt'.format(freq)),
                                        os.path.join(out_dir, 'map_tau_{}.txt'.format(freq)),
                                        angle, los_angle, number_of_pixels, out_dir=out_dir,
                                        name="{}_GHz".format(freq),
                                        title="LOS={} Angle={} B={} N={} freq={}".format(los_angle, angle, b, n, freq))
                        freq_difmap_models_dict[freq] = os.path.join(out_dir,
                                                                     'bk_{}.mdl'.format(freq))
                        freq_true_images_dict[freq] = os.path.join(out_dir,
                                                                   'map_i_{}.txt'.format(freq))

                    # Calculate shifts in different ways
                    observed_shifts = find_shift_from_difmap_models(freq_difmap_models_dict)
                    true_shifts = find_shifts_from_true_images(freq_true_images_dict,
                                                               number_of_pixels,
                                                               freq_pixel_size_dict)
                    bias_dr = observed_shifts['dr'][0][1] - true_shifts['dr'][0][1]
                    bias_bmaj = observed_shifts['bmaj'][0][1] - 1.0
                    bias_bmin = observed_shifts['bmin'][0][1] - 1.0

                    import matplotlib.pyplot as plt
                    # Plot all measured values
                    plt.plot(freqs, observed_shifts['dr'][1], '.k', ms=10,
                             label="k={0:.2f} observed dr".format(observed_shifts['dr'][0][1]))
                    plt.plot(freqs, true_shifts['dr'][1], '.r', ms=10,
                             label="k={0:.2f} true dr".format(true_shifts['dr'][0][1]))
                    plt.plot(freqs, observed_shifts['bmaj'][1], '.b', ms=10,
                             label="k={0:.2f} observed bmaj".format(observed_shifts['bmaj'][0][1]))
                    plt.plot(freqs, observed_shifts['bmin'][1], '.g', ms=10,
                             label="k={0:.2f} observed bmin".format(observed_shifts['bmin'][0][1]))

                    freqs_grid = np.linspace(freqs[0], freqs[-1], 100)
                    drs_observed_fit = shift_model(freqs_grid, observed_shifts['dr'][0][0],
                                                   observed_shifts['dr'][0][1])
                    bmajs_observed_fit = shift_model(freqs_grid, observed_shifts['bmaj'][0][0],
                                                     observed_shifts['bmaj'][0][1])
                    bmins_observed_fit = shift_model(freqs_grid, observed_shifts['bmin'][0][0],
                                                     observed_shifts['bmin'][0][1])
                    drs_true_fit = shift_model(freqs_grid, true_shifts['dr'][0][0],
                                               true_shifts['dr'][0][1])
                    plt.plot(freqs_grid, drs_observed_fit, 'k')
                    plt.plot(freqs_grid, bmajs_observed_fit, 'b')
                    plt.plot(freqs_grid, bmins_observed_fit, 'g')
                    plt.plot(freqs_grid, drs_true_fit, 'r')
                    plt.xlabel("Frequency, GHz")
                    plt.ylabel("Shift/Size, mas")
                    plt.legend(loc='best')
                    plt.title("LOS={} Angle={} B={} N={}".format(los_angle, angle, b, n))
                    plt.savefig(os.path.join(out_dir, 'fits.png'), bbox_inches='tight')
                    plt.close()

                    print("True k (dr) = {}".format(true_shifts['dr'][0][1]))
                    print("Bias of observed k (dr) = {}".format(bias_dr))
                    print("Bias of observed k (bmaj) = {}".format(bias_bmaj))
                    print("Bias of observed k (bmin) = {}".format(bias_bmin))