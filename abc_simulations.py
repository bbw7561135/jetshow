import os
import shutil
import glob
import json
import pickle
import numpy as np
from example_simple import (generate_sample_points, update_config,
                            find_core_separation_from_jet_using_difmap_model,
                            find_core_separation_from_center_using_simulations,
                            run_simulations,
                            fit_simulations_in_image_plane,
                            plot_stripe, plot_simulations_2d,
                            modelfit_simulation_result,
                            nested_dict,
                            FailedFindBestImageParamsException,
                            b_to_n_energy_ratio,
                            tb, tb_comp, t_syn,
                            plot_fitted_model, n1_equipartition,
                            log_base,
                            distance_from_SMBH, b_field)
from spydiff import clean_difmap, import_difmap_model
from image import plot as iplot
from image import find_bbox
from from_fits import create_clean_image_from_fits_file
from image_ops import rms_image
from scipy.ndimage.measurements import label
from scipy.ndimage.morphology import generate_binary_structure
from skimage.measure import regionprops
from skimage.morphology import opening
from collections import namedtuple


class TotalDisasterException(Exception):
    pass


Parameters = namedtuple("Parameters", ["b", "n", "los"])


def create_summary_from_result_dict(result_dict, dr_obs_true):
    """
    Create vector of summary statistics from results of simulations and some
    observed values.

    :param result_dict:
        Dictionary with results.
    :param dr_obs_true:
        Iterable of true values of core shift for all pairs of frequencies,
        e.g. ``(1.6GHz-15GHz, 8GHz-15GHz, 12GHz-15GHz)``
    :return:
        Vector of summary statistics. As last element is ``0`` for simulations
        that reproduce data completely, the observed value of summary statistic
        should have ``0`` as the last element.
    """
    freqs = result_dict.keys()
    freqs = sorted(freqs)
    dr_obs_true = np.array(dr_obs_true)
    flux_high = result_dict[freqs[-1]["result"]["flux_obs"]]
    size_high = result_dict[freqs[-1]["result"]["size_obs"]]
    dr_obs = np.array([result_dict[freq_low]["results"]["dr_obs"] -
                       result_dict[freqs[-1]]["results"]["dr_obs"] for
                       freq_low in freqs])
    dr = np.hypot(dr_obs, dr_obs_true)
    return flux_high, size_high, dr


def dist_matric(d, x):
    return np.sum((d - x)**2)


def cut_image(txt_file):
    """
    Function that cut image symmetrically removing zeros.
    :param txt_file:
        Files with image.
    """
    image = np.loadtxt(txt_file)
    s = generate_binary_structure(2, 2)
    labeled_array, num_features = label(image, structure=s)
    prop = regionprops(labeled_array, intensity_image=image)[0]
    np.savetxt(txt_file, prop.intensity_image)


def suggest_image_parameters(p, freq, pickle_history):
    """
    Suggest values of image size and pixel size using history.
    :param p:
        Instance of ``Parameter`` namedtuple.
    :param freq:
        Frequency [GHz].
    :param pickle_history:
        Path to pickle file with simulations history.
    :return:
        Tuple (image size, pixel size [mas])
    """
    with open(pickle_history, 'r') as fo:
        history = pickle.load(fo)
    parameters = history.keys()
    parametes = np.atleast_2d([[p.b, p.n, p.los] for p in parameters])


def abc_simulations(param,
                    z,
                    freqs,
                    uv_fits_templates,
                    cc_images,
                    pickle_history,
                    out_dir=None):
    """
    Simulation function for ABC. Simulates data given parameters (``b, n, los``)
    and returns the vector of the summary statistics.

    :param b:
        Value of the magnetic field at 1 pc [G].
    :param n:
        Value of particle density at 1 pc [cm^(-3)].
    :param los:
        Jet angle to the line of site [rad].
    :param z:
        Redshift of the source.
    :param freqs:
        Iterable of frequencies.
    :param uv_fits_templates:
        Iterable of paths to FITS files with self-calibrated uv-data. In order
        of ``freqs``.
    :param cc_images:
        Iterable of paths to FITS files with CLEAN maps. In order of ``freqs``.
    :param out_dir: (optional)
        Directory to store files. If ``None`` then use CWD. (default: ``None``)
    :param pickle_history:
        Path to pickle file with dictionary of simulations history.

    :return:
        Summary statistics. E.g. "observed" core flux at highest frequency,
        "observed" core size at highest frequency, "observed" ellipticity of the
        core at highest frequency, core shift between frequencies (distance
        between "core shift - frequency" curves).
    """
    if out_dir is None:
        out_dir = os.getcwd()
    else:
        # FIXME: Create directory if it doesn't exist
        if not os.path.exists(out_dir):
            pass

    result_dict = dict()
    b, n, los = np.exp(param)
    p = Parameters(b, n, los)
    with open(pickle_history, 'r') as fo:
        history = pickle.load(fo)
    print("Running simulations with b={}, n={}, los={}".format(b, n, los))
    for freq, uv_fits_template, cc_image in zip(freqs,
                                                uv_fits_templates,
                                                cc_images):

        # Cleaning old results if any
        simulated_maps_old = glob.glob(os.path.join(exe_dir, "map*.txt"))
        for to_remove in simulated_maps_old:
            os.unlink(to_remove)

        # Checking if simulations on highest frequency are done. If so then use
        # scaled image parameters
        if freqs[0] in history[p].keys():
            pixel_size_mas = history[p][freqs[0]]["pixel_size_mas"]*freqs[0]/freq
            number_of_pixels = history[p][freqs[0]]["number_of_pixels"]*freq/freqs[0]
            map_size = (number_of_pixels, pixel_size_mas)
            print("Using scaled image parameters: {}, {}".format(number_of_pixels,
                                                                 pixel_size_mas))
        else:
            pixel_size_mas = 0.01
            number_of_pixels = 400
            map_size = None

        update_dict = {"jet": {"bfield": {"parameters": {"b_1": b}},
                               "nfield": {"parameters": {"n_1": n}}},
                       "observation": {"los_angle": los,
                                       "frequency_ghz": freq,
                                       "redshift": z},
                       "image": {"pixel_size_mas": pixel_size_mas,
                                 "number_of_pixels": number_of_pixels}}
        update_config(cfg_file, update_dict)

        # FIXME: Handle ``FailedFindBestImageParamsException`` during ABC run
        simulation_params = run_simulations(cfg_file, path_to_executable,
                                            map_size=map_size)

        # Find total flux on simulated image
        image = os.path.join(exe_dir, "map_i.txt")
        image = np.loadtxt(image)

        # Rare case of strange fluxes
        image[image < 0] = 0
        image[image > 10.0] = 0

        # Total model flux at current frequency
        total_flux = image.sum()
        # Maximum model pixel flux at current frequency
        max_flux = image.max()
        cc_image = create_clean_image_from_fits_file(cc_image)
        noise_factor = 1.0

        initial_dfm_model = os.path.join(main_dir, 'initial_cg.mdl')
        out_dfm_model_fn = "bk_{}.mdl".format(freq)
        uv_fits_save_fname = "bk_{}.fits".format(freq)
        modelfit_simulation_result(exe_dir, initial_dfm_model,
                                   noise_factor=noise_factor,
                                   out_dfm_model_fn=out_dfm_model_fn,
                                   out_dir=out_dir,
                                   params=simulation_params,
                                   uv_fits_save_fname=uv_fits_save_fname,
                                   uv_fits_template=uv_fits_template)

        # Find measured and true distance of core to jet component
        dr_obs = find_core_separation_from_jet_using_difmap_model(os.path.join(out_dir,
                                                                               out_dfm_model_fn))
        dr_true = find_core_separation_from_center_using_simulations(os.path.join(exe_dir, "map_i.txt"),
                                                                     simulation_params)

        # This means something wrong with jetshow
        if total_flux < 0:
            print("b = {}, n = {}, los = {}".format(b, n, los))
            open(os.path.join(out_dir,
                              "total_disaster_{}_{}_{}_{}.txt".format(b, n, los, freq)), 'a').close()
            raise TotalDisasterException

        # Plot map with components superimposed
        cc_fits_save_fname = "bk_cc_{}.fits".format(freq)

        # For 18 cm we need large pixel size
        cellsize = 0.1
        if freq == 1.665:
            cellsize = 0.5
        elif freq == 8.1:
            cellsize = 0.2

        clean_difmap(uv_fits_save_fname, cc_fits_save_fname, 'I',
                     (1024, cellsize), path=out_dir,
                     path_to_script=path_to_script, show_difmap_output=False,
                     outpath=out_dir)

        ccimage = create_clean_image_from_fits_file(os.path.join(out_dir,
                                                                 cc_fits_save_fname))
        beam = ccimage.beam
        rms = rms_image(ccimage)
        blc, trc = find_bbox(ccimage.image, rms, 10)
        comps = import_difmap_model(out_dfm_model_fn, out_dir)

        plot_fitted_model(os.path.join(out_dir, uv_fits_save_fname), comps,
                          savefig=os.path.join(out_dir,
                                               "difmap_model_uvplot_{}.png".format(freq)))

        fig = iplot(ccimage.image, x=ccimage.x, y=ccimage.y, min_abs_level=3*rms,
                    beam=beam, show_beam=True, blc=blc, trc=trc, components=comps,
                    close=True, colorbar_label="Jy/beam")
        fig.savefig(os.path.join(out_dir, "cc_{}.png".format(freq)))

        # Move simulated images to data directory
        cut_image(os.path.join(exe_dir, "map_i.txt"))
        cut_image(os.path.join(exe_dir, "map_l.txt"))
        cut_image(os.path.join(exe_dir, "map_q.txt"))
        cut_image(os.path.join(exe_dir, "map_u.txt"))
        cut_image(os.path.join(exe_dir, "map_v.txt"))
        cut_image(os.path.join(exe_dir, "map_tau.txt"))
        for name in ('i', 'q', 'u', 'v', 'tau', 'l'):
            shutil.move(os.path.join(exe_dir, "map_{}.txt".format(name)),
                        os.path.join(out_dir, "map_{}_{}.txt".format(name,
                                                                     freq)))

        # Calculate some info
        dr_pc = distance_from_SMBH(dr_true, los, z=z)
        b_core = b_field(b, dr_pc)
        t_syn_years = t_syn(b_core, freq)/(np.pi*10**7)
        result_dict[freq] = dict()
        to_results = {"dr_obs": dr_obs,
                      "dr_true": dr_true,
                      "flux": total_flux,
                      "flux_obs": comps[0].p[0],
                      "bmaj_obs": comps[0].p[3],
                      "tb_difmap": np.log10(tb_comp(comps[0].p[0], comps[0].p[3], freq, z=z)),
                      "tb_pix": np.log10(tb(max_flux, freq, simulation_params[u'image'][u'pixel_size_mas'], z=z)),
                      "b_core": b_core,
                      "dr_core_pc": dr_pc,
                      "t_syn_core": t_syn_years,
                      "pixel_size_mas": simulation_params[u'image'][u'pixel_size_mas'],
                      "number_of_pixels": simulation_params[u'image'][u'number_of_pixels']}
        result_dict[freq] = to_results

        history[p] = result_dict
        with open(pickle_history, 'w') as fo:
            pickle.dump(history, fo)

    return create_summary_from_result_dict(result_dict, (0.3, 0.2, 0.1))


if __name__ == '__main__':
    main_dir = '/home/ilya/github/bck/jetshow'
    out_dir = '/home/ilya/github/bck/jetshow/test_abc'
    path_to_executable = os.path.join(main_dir, 'cmake-build-debug', 'jetshow')
    exe_dir, exe = os.path.split(path_to_executable)
    cfg_file = os.path.join(main_dir, 'config.json')
    uv_fits_templates_dir = '/home/ilya/github/bck/jetshow/uvf'
    uv_fits_template_u = os.path.join(uv_fits_templates_dir,
                                      '0235+164.u.2006_06_15.uvf')
    uv_fits_template_x = os.path.join(uv_fits_templates_dir,
                                      '0235+164.x.2006_06_15.uvf')
    uv_fits_template_j = os.path.join(uv_fits_templates_dir,
                                      '0235+164.j.2006_06_15.uvf')
    uv_fits_template_18 = os.path.join(uv_fits_templates_dir,
                                       '0235+164.18cm.2010_06_23.uvf')
    cc_image_u = os.path.join(uv_fits_templates_dir, '0235+164.u.2006_06_15_cc.fits')
    cc_image_x = os.path.join(uv_fits_templates_dir, '0235+164.x.2006_06_15_cc.fits')
    cc_image_j = os.path.join(uv_fits_templates_dir, '0235+164.j.2006_06_15_cc.fits')
    cc_image_18 = os.path.join(uv_fits_templates_dir, '0235+164.18cm.2010_06_23_cc.fits')

    path_to_script = '/home/ilya/github/vlbi_errors/difmap/final_clean_nw'
    pickle_history = os.path.join(out_dir, 'abc_history.pkl')

    if not os.path.exists(pickle_history):
        with open(pickle_history, 'w') as fo:
            pickle.dump(nested_dict(), fo)

    b = np.log(1.0)
    n = np.log(100.0)
    los = np.log(0.1)
    z = 0.3
    param = (b, n, los)
    summary = abc_simulations(param,
                              z,
                              freqs=(15.4, 12.1, 8.1),
                              uv_fits_templates=(uv_fits_template_u,
                                                 uv_fits_template_j,
                                                 uv_fits_template_x),
                              cc_images=(cc_image_u,
                                         cc_image_j,
                                         cc_image_x),
                              pickle_history=pickle_history,
                              out_dir=out_dir)
    # Observed values
    data = np.array([0.5, 0.1, 0.0])
    priors = [('uniform', [-3, 1]), ('uniform', [-2, 5])]
    prop = {'dfunc':dist_metric, 'outfile':"gaussian_example.txt", 'verbose':1, 'adapt_t': True, 'pert_kernel':2}
    sampler = astroabc.ABC_class(3, 100, data,[0.7,0.05], 20,priors,**prop)
    sampler.sample(abc_simulations)