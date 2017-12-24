import os
import shutil
import glob
import json
import numpy as np
from example_simple import (generate_sample_points, update_config,
                            find_core_separation_from_jet_using_difmap_model,
                            find_core_separation_from_center_using_simulations,
                            run_simulations,
                            fit_simulations_in_image_plane,
                            plot_stripe, plot_simulations_2d,
                            modelfit_simulation_result,
                            automodelfit_simulation_result,
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


class TotalDisasterException(Exception):
    pass


use_equipartition = False
# base = 2.0
# B = 0.01-10 Gs
# jet.bfield.parameters.b1
# b_uniform_borders = [log_base(0.01, base), log_base(10, base)]
# N = 50-5000 cm**(-3)
# jet.nfield.parameters.n1
# n_uniform_borders = [log_base(1, base), log_base(1000, base)]
# LOS = 1/(2\Gamma) - 2/(/Gamma)
# observation.los_angle
# los_uniform_borders = [log_base(2.865*np.pi/180, base),
#                        log_base(11.459*np.pi/180, base)]
# from itertools import product
# b_values = [0.01, 0.1, 1, 10]
# n_values = [10, 100, 1000, 5000]
# los_values = np.linspace(2.865*np.pi/180, 3*2.865*np.pi/180, 3)

# sample = list()
# for b_ in b_values:
#     for n_ in n_values:
#         if use_equipartition:
#             n_ = n1_equipartition(b_)
#         for los_ in los_values:
#             sample.append([b_, n_, los_])

# That is master sample
# sample= [[0.01, 10000, 0.050003683069637546],
#          [0.01, 10000, 0.10000736613927511],
#          [0.01, 10000, 0.15001104920891264],
#          [0.1, 10, 0.050003683069637546],
#          [0.1, 10, 0.10000736613927511],
#          [0.1, 10, 0.15001104920891264],
#          [0.1, 100, 0.050003683069637546],
#          [0.1, 100, 0.10000736613927511],
#          [0.1, 100, 0.15001104920891264],
#          [0.1, 1000, 0.050003683069637546],
#          [0.1, 1000, 0.10000736613927511],
#          [0.1, 1000, 0.15001104920891264],
#          [1, 1, 0.050003683069637546],
#          [1, 1, 0.10000736613927511],
#          [1, 1, 0.15001104920891264],
#          [1, 10, 0.050003683069637546],
#          [1, 10, 0.10000736613927511],
#          [1, 10, 0.15001104920891264],
#          [1, 100, 0.050003683069637546],
#          [1, 100, 0.10000736613927511],
#          [1, 100, 0.15001104920891264],
#          [1, 1000, 0.050003683069637546],
#          [1, 1000, 0.10000736613927511],
#          [1, 1000, 0.15001104920891264],
#          [10, 0.1, 0.050003683069637546],
#          [10, 0.1, 0.10000736613927511],
#          [10, 0.1, 0.15001104920891264],
#          [10, 1, 0.050003683069637546],
#          [10, 1, 0.10000736613927511],
#          [10, 1, 0.15001104920891264],
#          [10, 10, 0.050003683069637546],
#          [10, 10, 0.10000736613927511],
#          [10, 10, 0.15001104920891264]]


sample = [[1, 1000, 0.10000736613927511],
          [1, 10000, 0.10000736613927511],
          [1, 10000, 0.15001104920891264],
          [1, 1000, 0.050003683069637546],
          [0.1, 1000, 0.050003683069637546],
          [0.1, 10000, 0.050003683069637546],
          [2, 1000, 0.10000736613927511],
          [0.5, 1000, 0.050003683069637546],
          [1, 5000, 0.10000736613927511],
          [2, 5000, 0.15001104920891264],
          [10, 1000, 0.10000736613927511],
          [5, 1000, 0.10000736613927511],
          [5, 5000, 0.15001104920891264],
          [0.5, 5000, 0.050003683069637546],
          [0.5, 5000, 0.10000736613927511],
          [0.25, 10000, 0.10000736613927511],
          [0.25, 1000, 0.050003683069637546],
          [0.5, 10000, 0.10000736613927511],
          [0.75, 10000, 0.10000736613927511],
          [1.25, 5000, 0.10000736613927511],
          [1.25, 5000, 0.15001104920891264],
          [0.75, 5000, 0.10000736613927511],
          [1.75, 5000, 0.10000736613927511],
          [1.75, 5000, 0.050003683069637546],
          [0.25, 50000, 0.10000736613927511],
          [2, 5000, 0.10000736613927511],
          [1.5, 1000, 0.10000736613927511],
          [0.1, 50000, 0.10000736613927511],
          [0.5, 50000, 0.15001104920891264],
          [0.125, 50000, 0.10000736613927511]]

sample = np.atleast_2d(sample)
b_values = sample[:, 0]
n_values = sample[:, 1]
los_values = sample[:, 2]

# if not use_equipartition:
#     uniform_borders = [b_uniform_borders, los_uniform_borders, n_uniform_borders]
# else:
#     uniform_borders = [b_uniform_borders, los_uniform_borders]
#
# sample = generate_sample_points(uniform_borders, n_samples=100)
#
# b_values = base**sample[:, 0]
# los_values = base**sample[:, 1]
# if not use_equipartition:
#     n_values = base**sample[:, 2]
# else:
#     n_values = [n1_equipartition(b) for b in b_values]

# Run simulation
# Only simulations result with flux between these two are analyzed
total_flux_min = 0.1
total_flux_max = 3.0
# We adding noise that is equal to
# ``noise_scale * (model flux / original flux) * original_noise``
noise_scale = 0.1
main_dir = '/home/ilya/github/bck/jetshow'
out_dir = '/home/ilya/github/bck/jetshow/uvf_mf_adds'
path_to_executable = os.path.join(main_dir, 'cmake-build-debug', 'jetshow')
exe_dir, exe = os.path.split(path_to_executable)
cfg_file = os.path.join(main_dir, 'config.json')
uv_fits_template = '/home/ilya/github/bck/jetshow/uvf'
uv_fits_template_u = os.path.join(uv_fits_template,
                                '0235+164.u.2006_06_15.uvf')
uv_fits_template_x = os.path.join(uv_fits_template,
                                '0235+164.x.2006_06_15.uvf')
uv_fits_template_j = os.path.join(uv_fits_template,
                                '0235+164.j.2006_06_15.uvf')
uv_fits_template_18 = os.path.join(uv_fits_template,
                                   '0235+164.18cm.2010_06_23.uvf')
cc_image_u = os.path.join(uv_fits_template, '0235+164.u.2006_06_15_cc.fits')
cc_image_x = os.path.join(uv_fits_template, '0235+164.x.2006_06_15_cc.fits')
cc_image_j = os.path.join(uv_fits_template, '0235+164.j.2006_06_15_cc.fits')
cc_image_18 = os.path.join(uv_fits_template, '0235+164.18cm.2010_06_23_cc.fits')


path_to_script = '/home/ilya/github/vlbi_errors/difmap/final_clean_nw'
json_out = os.path.join(out_dir, 'history_mf.json')

if not os.path.exists(json_out):
    with open(json_out, 'w') as fo:
        json.dump(nested_dict(), fo)

i = 1
for b, n, los in zip(b_values, n_values, los_values):
    print("Running simulations with b={}, n={}, los={}".format(b, n, los))
    for freq, uv_fits_template, cc_image in zip((15.4, 12.1, 8.1, 1.665),
                                                (uv_fits_template_u,
                                                 uv_fits_template_j,
                                                 uv_fits_template_x,
                                                 uv_fits_template_18),
                                                (cc_image_u,
                                                 cc_image_j,
                                                 cc_image_x,
                                                 cc_image_18)):

        # Cleaning old results if any
        simulated_maps_old = glob.glob(os.path.join(exe_dir, "map*.txt"))
        for to_remove in simulated_maps_old:
            os.unlink(to_remove)

        update_dict = {"jet": {"bfield": {"parameters": {"b_1": b}},
                               "nfield": {"parameters": {"n_1": n}}},
                       "observation": {"los_angle": los,
                                       "frequency_ghz": freq},
                       "image": {"pixel_size_mas": 0.01, "number_of_pixels": 400}}
        update_config(cfg_file, update_dict)

        # If we failed to find best image params - just continue
        try:
            simulation_params = run_simulations(cfg_file, path_to_executable)
        except FailedFindBestImageParamsException:
            open(os.path.join(out_dir,
                              "failed_find_best_image_params_{}_{}_{}_{}.txt".format(b, n, los, freq)), 'a').close()
            break

        # Find total flux on simulated image
        image = os.path.join(exe_dir, "map_i.txt")
        image = np.loadtxt(image)

        # Rare case of strange fluxes
        image[image < 0] = 0
        image[image > 10.0] = 0

        total_flux = image.sum()
        print("TOTAL FLUX = {}".format(total_flux))
        if freq == 15.4:
            if total_flux > total_flux_max or total_flux < total_flux_min:
                print("SKIPPING...")
                open(os.path.join(out_dir, "bad_total_flux_{}_{}_{}_{}.txt".format(total_flux, b, n, los)),'a').close()
                break
        max_flux = image.max()
        cc_image = create_clean_image_from_fits_file(cc_image)
        noise_factor = noise_scale*image.sum()/cc_image.cc.sum()

        initial_dfm_model = os.path.join(main_dir, 'initial_cg.mdl')
        out_dfm_model_fn = "bk_{}_{}.mdl".format(str(i).zfill(2), freq)
        uv_fits_save_fname = "bk_{}_{}.fits".format(str(i).zfill(2), freq)
        modelfit_simulation_result(exe_dir, initial_dfm_model,
                                   noise_factor=noise_factor,
                                   out_dfm_model_fn=out_dfm_model_fn,
                                   out_dir=out_dir,
                                   params=simulation_params,
                                   uv_fits_save_fname=uv_fits_save_fname,
                                   uv_fits_template=uv_fits_template)

        # Modelfit results of simulations by first FT model to template uv-fits,
        # adding specified fraction of the observed noise and automodel
        # automodelfit_simulation_result(exe_dir,
        #                                noise_factor=noise_factor,
        #                                out_dfm_model_fn=out_dfm_model_fn,
        #                                out_dir=out_dir,
        #                                params=simulation_params,
        #                                uv_fits_save_fname=uv_fits_save_fname,
        #                                uv_fits_template=uv_fits_template,
        #                                core_elliptic=False,
        #                                n_max_components=20,
        #                                mapsize_clean=(512, 0.1))

        # Find measured and true distance of core to jet component
        dr_obs = find_core_separation_from_jet_using_difmap_model(os.path.join(out_dir,
                                                                               out_dfm_model_fn))
        dr_true = find_core_separation_from_center_using_simulations(os.path.join(exe_dir, "map_i.txt"),
                                                                     simulation_params)

        # This means smth wrong with jetshow
        if total_flux < 0:
            print("b = {}, n = {}, los = {}".format(b, n, los))
            open(os.path.join(out_dir,
                              "total_disaster_{}_{}_{}_{}.txt".format(b, n, los, freq)), 'a').close()
            # raise TotalDisasterException
            continue

        # Plot map with components superimposed
        cc_fits_save_fname = "bk_cc_{}_{}.fits".format(str(i).zfill(2),
                                                       freq)

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
                                               "difmap_model_uvplot_{}_{}.png".format(str(i).zfill(2), freq)))

        fig = iplot(ccimage.image, x=ccimage.x, y=ccimage.y, min_abs_level=3*rms,
                    beam=beam, show_beam=True, blc=blc, trc=trc, components=comps,
                    close=True, colorbar_label="Jy/beam")
        fig.savefig(os.path.join(out_dir, "cc_{}_{}.png".format(str(i).zfill(2),
                                                                freq)))

        g = fit_simulations_in_image_plane(os.path.join(exe_dir, "map_i.txt"),
                                           simulation_params, core=comps[0])
        # Find gaussian parameters
        g_flux = 2.0*np.pi*g.amplitude.value*g.x_stddev.value*g.y_stddev.value
        g_bmaj = 2.0*np.sqrt(2.0*np.log(2)) * g.x_stddev

        plot_stripe(os.path.join(exe_dir, "map_i.txt"),
                    os.path.join(out_dir, out_dfm_model_fn), simulation_params,
                    g=g, close=True, savefig=os.path.join(out_dir, "stripe_{}_{}.png".format(str(i).zfill(2),
                                                                                             freq)))

        plot_simulations_2d(os.path.join(exe_dir, "map_i.txt"),
                            simulation_params, core=comps[0], g=g, close=True,
                            savefig=os.path.join(out_dir, "sim2d_{}_{}.png".format(str(i).zfill(2),
                                                                                   freq)))

        # Move simulated images to data directory
        for name in ('i', 'q', 'u', 'v', 'tau', 'l'):
            shutil.move(os.path.join(exe_dir, "map_{}.txt".format(name)),
                        os.path.join(out_dir, "map_{}_{}_{}.txt".format(name,
                                                                        str(i).zfill(2),
                                                                        freq)))

        # Calculate some info
        dr_pc = distance_from_SMBH(dr_true, los, z=0.5)
        b_core = b_field(b, dr_pc)
        t_syn_years = t_syn(b_core, freq)/(np.pi*10**7)

        with open(json_out, 'r') as fo:
            history = json.load(fo)
        history["{}_{}".format(str(i).zfill(2), freq)] = {"parameters": {"b": b,
                                                                         "n": n,
                                                                         "los": los,
                                                                         "sigma": b_to_n_energy_ratio(b, n)},
                                                          "results": {"dr_obs": dr_obs,
                                                                      "dr_true": dr_true,
                                                                      "flux": total_flux,
                                                                      "tb_difmap": np.log10(tb_comp(comps[0].p[0], comps[0].p[3], freq, z=0.5)),
                                                                      "tb_pix": np.log10(tb(max_flux, freq, simulation_params[u'image'][u'pixel_size_mas'], z=0.5)),
                                                                      "tb_gsim": np.log10(tb_comp(g_flux, g_bmaj, freq, z=0.5)),
                                                                      "b_core": b_core,
                                                                      "dr_core_pc": dr_pc,
                                                                      "t_syn_core": t_syn_years},
                                                          "image": {"pixel_size_mas": simulation_params[u'image'][u'pixel_size_mas'],
                                                                    "number_of_pixels": simulation_params[u'image'][u'number_of_pixels']}}
        with open(json_out, 'w') as fo:
            json.dump(history, fo)

     # Update counter
    i += 1