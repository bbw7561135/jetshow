import os
import json
from example_simple import (find_tb_for_given_keys,
                            select_keys_with_given_frequency)


freq = 15.4
out_dir = '/home/ilya/github/bck/jetshow/uvf_mf_adds'
json_history = os.path.join(out_dir, 'history_mf.json')

with open(json_history, "r") as fo:
    data = json.load(fo)

keys = data.keys()
keys = select_keys_with_given_frequency(keys, 15.4)
tb = find_tb_for_given_keys(json_history, keys)

tb_true = tb[:, 0]
tb_difmap = tb[:, 1]

# Histograms of Tb
fig, axes = plt.subplots(1, 1)
axes.hist(tb_difmap, alpha=0.5, range=[10, 12.5], label='difmap')
axes.hist(tb_true, alpha=0.5, range=[10, 12.5], label='simulation pixel')
axes.set_xlabel("lg of z-corrected Tb", fontsize=14)
axes.set_ylabel("N", fontsize=14)
axes.tick_params(labelsize=12)
plt.legend(loc='best', fontsize=14)
fig.savefig('/home/ilya/github/bck/jetshow/uvf_mf_adds/Tb_15.4.png',
            bbox_inches='tight', dpi=300)

# # Histograms of bias Tb
# fig, axes = plt.subplots(1, 1)
# axes.hist(tb_true-tb_difmap)
# axes.set_xlabel("bias of lg of z-corrected Tb", fontsize=14)
# axes.set_ylabel("N", fontsize=14)
# axes.tick_params(labelsize=12)
# plt.legend(loc='best', fontsize=14)
# fig.savefig('/home/ilya/github/bck/jetshow/uvf_mf_adds/Tb_15.4_bias.png',
#             bbox_inches='tight', dpi=300)