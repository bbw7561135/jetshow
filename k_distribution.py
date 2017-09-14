import os
import json
from example_simple import parse_history_mf__


out_dir = '/home/ilya/github/bck/jetshow/uvf_mf_adds'
json_history = os.path.join(out_dir, 'history_mf.json')
X, y = parse_history_mf__(json_history)

k_obs = y[:, 0]
k_true = y[:, 1]
k_bias = y[:, 2]

# # Histograms of k
# fig, axes = plt.subplots(1, 1)
# axes.hist(k_obs, alpha=0.5, label='observed')
# axes.hist(k_true, alpha=0.5, label='true')
# axes.set_xlabel("k", fontsize=14)
# axes.set_ylabel("N", fontsize=14)
# axes.tick_params(labelsize=12)
# plt.legend(loc='best', fontsize=14)
# fig.savefig('/home/ilya/github/bck/jetshow/uvf_mf_adds/k_obs_true.png',
#             bbox_inches='tight', dpi=300)

# Histograms of k bias
fig, axes = plt.subplots(1, 1)
axes.hist(k_bias)
axes.set_xlabel("bias of k", fontsize=14)
axes.set_ylabel("N", fontsize=14)
axes.tick_params(labelsize=12)
plt.legend(loc='best', fontsize=14)
fig.savefig('/home/ilya/github/bck/jetshow/uvf_mf_adds/k_bias.png',
            bbox_inches='tight', dpi=300)