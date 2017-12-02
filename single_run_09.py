from example_simple import run_simulations, update_config
import os


main_dir = '/home/ilya/github/bck/jetshow'
path_to_executable = os.path.join(main_dir, 'cmake-build-debug', 'jetshow')
cfg_file = os.path.join(main_dir, 'config.json')
los = 0.1
b = 1.0
n = 5000.0
freq = 22.4

update_dict = {"jet": {"bfield": {"parameters": {"b_1": b}},
                       "nfield": {"parameters": {"n_1": n}}},
               "observation": {"los_angle": los,
                               "frequency_ghz": freq},
               "image": {"pixel_size_mas": 0.00253, "number_of_pixels": 756}}
update_config(cfg_file, update_dict)

simulation_params = run_simulations(cfg_file, path_to_executable)
