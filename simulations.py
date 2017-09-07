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
                            nested_dict,
                            FailedFindBestImageParamsException)
from spydiff import clean_difmap, import_difmap_model
from image import plot as iplot
from image import find_bbox
from from_fits import create_clean_image_from_fits_file
from image_ops import rms_image


class TotalDisasterException(Exception):
    pass


# B = 0.01-10 Gs
# jet.bfield.parameters.b1
b_uniform_borders = [np.log10(0.1), np.log10(10)]
# N = 50-5000 cm**(-3)
# jet.nfield.parameters.n1
n_uniform_borders = [np.log10(50), np.log10(5000)]
# LOS = 1/2G - 2/G
# observation.los_angle
los_uniform_borders = [2.865*np.pi/180, 11.459*np.pi/180]

uniform_borders = [b_uniform_borders, n_uniform_borders, los_uniform_borders]

sample = generate_sample_points(uniform_borders, n_samples=30)

b_values = 10.**sample[:, 0]
n_values = 10.**sample[:, 1]
los_values = sample[:, 2]

# Run simulation
main_dir = '/home/ilya/github/bck/jetshow'
out_dir = '/home/ilya/github/bck/jetshow/uvf'
path_to_executable = os.path.join(main_dir, 'cmake-build-debug', 'jetshow')
exe_dir, exe = os.path.split(path_to_executable)
cfg_file = os.path.join(main_dir, 'config.json')
uv_fits_template = '/home/ilya/github/bck/jetshow/uvf'
uv_fits_template = os.path.join(uv_fits_template,
                                '0235+164.u.2006_06_15.uvf')
path_to_script = '/home/ilya/github/vlbi_errors/difmap/final_clean_nw'
json_out = os.path.join(out_dir, 'history.json')
with open(json_out, 'w') as fo:
    json.dump(nested_dict(), fo)

i = 1
for b, n, los in zip(b_values, n_values, los_values):
    print("Running simulations with b={}, n={}, los={}".format(b, n, los))

    # Cleaning old results if any
    simulated_maps_old = glob.glob(os.path.join("exe_dir", "map*.txt"))
    for to_remove in simulated_maps_old:
        os.unlink(to_remove)

    update_dict = {"jet": {"bfield": {"parameters": {"b_1": b}},
                   "nfield": {"parameters": {"n_1": n}}},
                   "observation": {"los_angle": los},
                   "image": {"pixel_size_mas": 0.001, "number_of_pixels": 200}}
    update_config(cfg_file, update_dict)

    # If we failed to find best image params - just continue
    try:
        simulation_params = run_simulations(cfg_file, path_to_executable)
    except FailedFindBestImageParamsException:
        continue

    initial_dfm_model = os.path.join(main_dir, 'initial_cg.mdl')
    out_dfm_model_fn = "bk_{}.mdl".format(str(i).zfill(2))
    uv_fits_save_fname = "bk_{}.fits".format(str(i).zfill(2))
    modelfit_simulation_result(exe_dir, initial_dfm_model, noise_factor=0.01,
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

    # Find total flux on simulated image
    image = os.path.join(exe_dir, "map_i.txt")
    image = np.loadtxt(image)
    total_flux = image.sum()

    # This means smth wrong with jetshow
    if total_flux < 0:
        print("b = {}, n = {}, los = {}".format(b, n, los))
        # raise TotalDisasterException
        open(os.path.join(out_dir, "total_disaster_{}_{}_{}.txt".format(b, n, los)), 'a').close()
        continue

    # Plot map with components superimposed
    cc_fits_save_fname = "bk_cc_{}.fits".format(str(i).zfill(2))
    clean_difmap(uv_fits_save_fname, cc_fits_save_fname, 'I', (1024, 0.1),
                 path=out_dir, path_to_script=path_to_script,
                 show_difmap_output=False, outpath=out_dir)

    ccimage = create_clean_image_from_fits_file(os.path.join(out_dir,
                                                             cc_fits_save_fname))
    beam = ccimage.beam
    rms = rms_image(ccimage)
    blc, trc = find_bbox(ccimage.image, rms, 10)
    comps = import_difmap_model(out_dfm_model_fn, out_dir)
    fig = iplot(ccimage.image, x=ccimage.x, y=ccimage.y, min_abs_level=3*rms,
                beam=beam, show_beam=True, blc=blc, trc=trc, components=comps,
                close=True)
    fig.savefig(os.path.join(out_dir, "cc_{}.png".format(str(i).zfill(2))))

    g = fit_simulations_in_image_plane(os.path.join(exe_dir, "map_i.txt"),
                                       simulation_params, core=comps[0])

    plot_stripe(os.path.join(exe_dir, "map_i.txt"),
                os.path.join(out_dir, out_dfm_model_fn), simulation_params, g=g,
                close=True, savefig=os.p    ath.join(out_dir, "stripe_{}.png".format(str(i).zfill(2))))

    plot_simulations_2d(os.path.join(exe_dir, "map_i.txt"),
                        simulation_params, core=comps[0], g=g, close=True,
                        savefig=os.path.join(out_dir,
                                             "sim2d_{}.png".format(str(i).zfill(2))))

    # Move simulated images to data directory
    for name in ('i', 'q', 'u', 'v', 'tau', 'l'):
        shutil.move(os.path.join(exe_dir, "map_{}.txt".format(name)),
                    os.path.join(out_dir, "map_{}_{}.txt".format(name,
                                                                 str(i).zfill(2))))

    # # Save results to history
    # to_save = nested_dict()
    # to_save.update({"{}".format(str(i).zfill(2)): {"parameters": {"b": b,
    #                                                               "n": n,
    #                                                               "los": los},
    #                                                "results": {"dr_obs": dr_obs,
    #                                                            "dr_true": dr_true,
    #                                                            "flux": total_flux}}})
    with open(json_out, 'r') as fo:
        history = json.load(fo)
    history["{}".format(str(i).zfill(2))] = {"parameters": {"b": b,
                                                            "n": n,
                                                            "los": los},
                                             "results": {"dr_obs": dr_obs,
                                                         "dr_true": dr_true,
                                                         "flux": total_flux}}
    with open(json_out, 'w') as fo:
        json.dump(history, fo)

    # Update counter
    i += 1