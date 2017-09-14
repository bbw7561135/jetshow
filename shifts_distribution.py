import os
import json
from example_simple import (find_shifts_and_fluxes_for_given_keys_and_freqs,
                            strip_frequency_from_keys)


out_dir = '/home/ilya/github/bck/jetshow/uvf_mf_adds'
json_history = os.path.join(out_dir, 'history_mf.json')

with open(json_history, "r") as fo:
    data = json.load(fo)

keys = data.keys()
keys = strip_frequency_from_keys(keys)
# keys = ['09']

shifts, fluxes = find_shifts_and_fluxes_for_given_keys_and_freqs(json_history,
                                                                 keys, 8.1,
                                                                 15.4)
shifts_obs = shifts[:, 0]
shifts_true = shifts[:, 1]
shifts_bias = shifts[:, 2]

# Histograms of shifts
# fig, axes = plt.subplots(1, 1)
# axes.hist(shifts_obs, range=[0, 2.5], alpha=0.5, label='observed')
# axes.hist(shifts_true, range=[0, 2.5], alpha=0.5, label='true')
# axes.set_xlabel("1.7 - 8.1 GHz core shift, [mas]", fontsize=14)
# axes.set_ylabel("N", fontsize=14)
# axes.tick_params(labelsize=12)
# plt.legend(loc='best', fontsize=14)
# fig.savefig('/home/ilya/github/bck/jetshow/uvf_mf_adds/shifts_1.7_8.1.png',
#             bbox_inches='tight', dpi=300)

# # Dependence of bias on flux
# fig, axes = plt.subplots(1, 1)
# axes.scatter(fluxes, shifts_bias, marker="o")
# axes.set_ylabel("8.1 - 15.4 GHz core shift bias, [mas]", fontsize=14)
# axes.set_xlabel("15.4 GHz flux, [Jy]", fontsize=14)
# axes.tick_params(labelsize=12)
# fig.savefig('/home/ilya/github/bck/jetshow/uvf_mf_adds/shifts_bias_vs_flux_8.1_15.4.png',
#             bbox_inches='tight', dpi=300)

# # Dependence of observed shift on flux
# fig, axes = plt.subplots(1, 1)
# axes.scatter(fluxes, shifts_obs, marker="o")
# axes.set_ylabel("8.1 - 15.4 GHz observed core shift, [mas]", fontsize=14)
# axes.set_xlabel("15.4 GHz flux, [Jy]", fontsize=14)
# axes.tick_params(labelsize=12)
# fig.savefig('/home/ilya/github/bck/jetshow/uvf_mf_adds/shifts_obs_vs_flux_8.1_15.4.png',
#             bbox_inches='tight', dpi=300)

# Dependence of true shift on flux
fig, axes = plt.subplots(1, 1)
axes.scatter(fluxes, shifts_true, marker="o")
axes.set_ylabel("8.1 - 15.4 GHz true core shift, [mas]", fontsize=14)
axes.set_xlabel("15.4 GHz flux, [Jy]", fontsize=14)
axes.tick_params(labelsize=12)
fig.savefig('/home/ilya/github/bck/jetshow/uvf_mf_adds/shifts_true_vs_flux_8.1_15.4.png',
            bbox_inches='tight', dpi=300)