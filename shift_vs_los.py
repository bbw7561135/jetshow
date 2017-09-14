import os
import json
from example_simple import (find_shifts_and_los_for_given_keys_and_freqs,
                            strip_frequency_from_keys)


out_dir = '/home/ilya/github/bck/jetshow/uvf_mf_adds'
json_history = os.path.join(out_dir, 'history_mf.json')

with open(json_history, "r") as fo:
    data = json.load(fo)

keys = data.keys()
keys = strip_frequency_from_keys(keys)
# keys = ['09']

shifts, loses = find_shifts_and_los_for_given_keys_and_freqs(json_history,
                                                                 keys, 8.1,
                                                                 15.4)
shifts_obs = shifts[:, 0]
shifts_true = shifts[:, 1]
shifts_bias = shifts[:, 2]

loses += np.random.normal(0, 0.0003, size=len(loses))

# # Dependence of observed shift on los
# fig, axes = plt.subplots(1, 1)
# axes.scatter(loses, shifts_obs, marker="o")
# axes.set_ylabel("8.1 - 15.4 GHz observed core shift, [mas]", fontsize=14)
# axes.set_xlabel("LOS", fontsize=14)
# axes.tick_params(labelsize=12)
# fig.savefig('/home/ilya/github/bck/jetshow/uvf_mf_adds/shifts_obs_vs_los_8.1_15.4.png',
#             bbox_inches='tight', dpi=300)
#
# # Dependence of observed shift on los
# fig, axes = plt.subplots(1, 1)
# axes.scatter(loses, shifts_true, marker="o")
# axes.set_ylabel("8.1 - 15.4 GHz true core shift, [mas]", fontsize=14)
# axes.set_xlabel("LOS", fontsize=14)
# axes.tick_params(labelsize=12)
# fig.savefig('/home/ilya/github/bck/jetshow/uvf_mf_adds/shifts_true_vs_los_8.1_15.4.png',
#             bbox_inches='tight', dpi=300)

# Dependence of bias on los
fig, axes = plt.subplots(1, 1)
axes.scatter(loses, shifts_bias, marker="o")
axes.set_ylabel("8.1 - 15.4 GHz core shift bias, [mas]", fontsize=14)
axes.set_xlabel("LOS", fontsize=14)
axes.tick_params(labelsize=12)
fig.savefig('/home/ilya/github/bck/jetshow/uvf_mf_adds/shifts_bias_vs_los_8.1_15.4.png',
            bbox_inches='tight', dpi=300)