import json
import glob
import os
import matplotlib.pyplot as plt
import matplotlib
label_size = 14
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['xtick.labelsize'] = label_size
matplotlib.rcParams['ytick.labelsize'] = label_size
matplotlib.rcParams['axes.titlesize'] = label_size
matplotlib.rcParams['axes.labelsize'] = label_size
matplotlib.rcParams['font.size'] = label_size
matplotlib.rcParams['legend.fontsize'] = label_size
from example_simple import (find_shifts_and_fluxes_for_given_keys_and_freqs,
                            find_shifts_and_los_for_given_keys_and_freqs,
                            strip_frequency_from_keys,
                            find_params_for_given_keys,
                            convert_json_history_to_params_dictionary,
                            find_core_separation_from_jet_using_difmap_model_,
                            find_flux_for_given_keys)


models_dir = '/home/ilya/github/bck/jetshow/automodel/'
json_history = '/home/ilya/github/bck/jetshow/uvf_mf_adds/history_mf.json'

with open(json_history, "r") as fo:
    data = json.load(fo)

keys = data.keys()
params_dict = convert_json_history_to_params_dictionary(json_history)
keys = strip_frequency_from_keys(keys)

shifts, fluxes = find_shifts_and_fluxes_for_given_keys_and_freqs(json_history,
                                                                 keys, 8.1,
                                                                 15.4)
# For all cases
shifts_true = shifts[:, 1]
# These ones are for 1 GC
shifts_obs = shifts[:, 0]
shift_bias_1cg = shifts_obs - shifts_true

# Now parse external model and get shifts from them
freq_low = '8_1'
freq_high = '15_4'
shifts_params_dict = {}
for key in keys:
    params = params_dict[key]
    key_low = u'{}_{}'.format(key, freq_low)
    key_high = u'{}_{}'.format(key, freq_high)
    model_low = glob.glob(os.path.join(models_dir, "bk_{}_fits_cg_fitted_*.mdl".format(key_low)))[0]
    model_high = glob.glob(os.path.join(models_dir, "bk_{}_fits_cg_fitted_*.mdl".format(key_high)))[0]
    dr_obs_low = find_core_separation_from_jet_using_difmap_model_(model_low)
    dr_obs_high = find_core_separation_from_jet_using_difmap_model_(model_high)
    shift_obs = dr_obs_high - dr_obs_low
    shifts, fluxes = find_shifts_and_fluxes_for_given_keys_and_freqs(json_history, [key], 8.1, 15.4)
    shift_true = shifts[:, 1][0]
    flux = fluxes[0]
    shifts_params_dict.update({key: [shift_true, shift_obs, flux]})


shifts_obs = np.atleast_2d(shifts_params_dict.values())[:, 1]
shifts_true = np.atleast_2d(shifts_params_dict.values())[:, 0]
shifts_bias = shifts_obs - shifts_true

# bias 1cg = 0.0968 (median=0.096, max=0.2, min=0.027) +/- 0.04
# bias 1eg = 0.0864 (median=0.079, max=0.205, min=0.027) +/- 0.04
# bias cgauto = 0.0706 (median=0.068, max=0.117, min=0.027) +/- 0.019
# bias egauto = 0.0889 (median=0.087, max=0.166, min=0.027) +/- 0.036

# Histogram of shifts
fig, axes = plt.subplots(1, 1)
# axes.hist(shifts_obs, range=[0, 0.3], alpha=0.5, label='observed')
axes.hist(shift_bias_1cg, range=[0, 0.15], alpha=0.25, label='1 CG')
axes.hist(shift_bias_cgauto, range=[0, 0.15], alpha=0.25, label='auto CG')
axes.hist(shift_bias_1eg, range=[0, 0.15], alpha=0.25, label='1 EG')
axes.hist(shift_bias_egauto, range=[0, 0.15], alpha=0.25, label='EG + auto CG')

# axes.hist(shifts_true, range=[0, 0.3], alpha=0.5, label='true')
axes.set_xlabel("8.1 - 15.4 GHz core shift bias, [mas]", fontsize=14)
axes.set_ylabel("N", fontsize=14)
axes.tick_params(labelsize=12)
plt.legend(loc='best', fontsize=14)
fig.savefig('/home/ilya/github/bck/jetshow/uvf_mf_adds/shifts_8.1_15.4_auto.png',
            bbox_inches='tight', dpi=300)


# Plot the bias violin plot
labels = ('CG', 'EG', 'auto CG', 'EG + auto CG')
colors = [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728']
double_colors = [color for color in colors for _ in (0, 1)]
# Create a figure instance
fig = plt.figure()

# Create an axes instance
ax = fig.add_subplot(111)
boxes = [shift_bias_1cg, shift_bias_1eg, shift_bias_cgauto, shift_bias_egauto]
bp = ax.boxplot(boxes, notch=True, vert=True, bootstrap=1000, whis=1.5)
for i, box in enumerate(bp['boxes']):
    # change outline color
    box.set(color=colors[i], linewidth=2)
    # change fill color
    # box.set(facecolor=colors[i], alpha=0.5)

## change color and linewidth of the whiskers
for i, whisker in enumerate(bp['whiskers']):
    whisker.set(linewidth=2, color=double_colors[i])

## change color and linewidth of the caps
for i, cap in enumerate(bp['caps']):
    cap.set(linewidth=2, color=double_colors[i])

## change color and linewidth of the medians
for i, median in enumerate(bp['medians']):
    median.set(linewidth=2, color=colors[i])

## change the style of fliers and their fill
for flier in bp['fliers']:
    flier.set(marker='o', color=colors[i])

## Custom x-axis labels
ax.set_xticklabels(labels, fontsize=14)
# ax.set_yticks([0.7, 0.8])
# ax.set_yticklabels([0.025, 0.21], fontsize=14)
## Remove top axes and right axes ticks
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

ax.set_ylabel("Core shift bias, [mas]", fontsize=16)
fig.show()


# Dependence of bias on flux
shifts, fluxes = find_shifts_and_fluxes_for_given_keys_and_freqs(json_history,
                                                                 keys, 8.1,
                                                                 15.4)
fig, axes = plt.subplots(1, 1)
axes.scatter(fluxes, shifts_bias, marker="o")
axes.set_ylabel("8.1 - 15.4 GHz core shift bias, [mas]", fontsize=14)
axes.set_xlabel("15.4 GHz flux, [Jy]", fontsize=14)
axes.tick_params(labelsize=12)
fig.savefig('/home/ilya/github/bck/jetshow/uvf_mf_adds/shifts_bias_vs_flux_8.1_15.4_auto.png',
            bbox_inches='tight', dpi=300)

# Dependence of observed shift on flux
fig, axes = plt.subplots(1, 1)
axes.scatter(fluxes, shifts_obs, marker="o")
axes.set_ylabel("8.1 - 15.4 GHz observed core shift, [mas]", fontsize=14)
axes.set_xlabel("15.4 GHz flux, [Jy]", fontsize=14)
axes.tick_params(labelsize=12)
fig.savefig('/home/ilya/github/bck/jetshow/uvf_mf_adds/shifts_obs_vs_flux_8.1_15.4_auto.png',
            bbox_inches='tight', dpi=300)

# Dependence of true shift on flux
fig, axes = plt.subplots(1, 1)
axes.scatter(fluxes, shifts_true, marker="o")
axes.set_ylabel("8.1 - 15.4 GHz true core shift, [mas]", fontsize=14)
axes.set_xlabel("15.4 GHz flux, [Jy]", fontsize=14)
axes.tick_params(labelsize=12)
fig.savefig('/home/ilya/github/bck/jetshow/uvf_mf_adds/shifts_true_vs_flux_8.1_15.4_auto.png',
            bbox_inches='tight', dpi=300)

# Dependense on the parameters: logbias = 0.32*logb + 0.13*logn - 3.59 (1/gamma)
#  0.62,  0.27, -5.05 for 1/2gamma
#  0.18,  0.07, -2.88 for 2/gamma
keys_1gamma = [key for key, value in params_dict.items()
               if value[-1] == 0.1000073661392751]
keys_halfgamma = [key for key, value in params_dict.items()
                  if value[-1] == 0.050003683069637546]
keys_dgamma = [key for key, value in params_dict.items()
                  if value[-1] == 0.15001104920891264]


biases_5000 = []
bs_5000 = []
biases_1000 = []
bs_1000 = []

for key in keys_1gamma:
    if params_dict[key][1] == 5000.0:
        # print(shifts_params_dict[key][1]-shifts_params_dict[key][0], params_dict[key][0], params_dict[key][1])
        biases_5000.append(shifts_params_dict[key][1]-shifts_params_dict[key][0])
        bs_5000.append(params_dict[key][0])
    elif params_dict[key][1] == 1000.0:
        # print(shifts_params_dict[key][1]-shifts_params_dict[key][0], params_dict[key][0], params_dict[key][1])
        biases_1000.append(shifts_params_dict[key][1]-shifts_params_dict[key][0])
        bs_1000.append(params_dict[key][0])

biases_dgamma = []
fluxes_dgamma = []
bs = []
ns = []
for key in keys_dgamma:
    print(shifts_params_dict[key][1]-shifts_params_dict[key][0], params_dict[key][0], params_dict[key][1])
    biases_dgamma.append(shifts_params_dict[key][1]-shifts_params_dict[key][0])
    fluxes_dgamma.append(shifts_params_dict[key][2])
    bs.append(params_dict[key][0])
    ns.append(params_dict[key][1])

logbiases = np.log(biases)
logbiases
logbs = log(bs)
logns = log(ns)
import scipy


def f(x, logb, logn, logbias):
    alpha = x[0]
    beta = x[1]
    interc = x[2]
    modelled = alpha*logb + beta*logn + interc
    diffs = modelled - logbias
    return diffs.flatten()
result = scipy.optimize.leastsq(f, [1.0, 1.0, 0.0],
                                args=(logbs , logns, logbiases))


fluxes_halfgamma = find_flux_for_given_keys(json_history, [key+"_15.4" for key in keys_halfgamma])
fluxes_1gamma = find_flux_for_given_keys(json_history, [key+"_15.4" for key in keys_1gamma])
fluxes_dgamma = find_flux_for_given_keys(json_history, [key+"_15.4" for key in keys_dgamma])
plt.plot(fluxes_halfgamma, biases_halfgamma, 'o', label=r"$1/2\Gamma$")
plt.plot(fluxes_1gamma, biases_1gamma, 'o', label=r"$1/\Gamma$")
plt.plot(fluxes_dgamma, biases_dgamma, 'o', label=r"$2/\Gamma$")
plt.legend(loc="best")
plt.xlabel(r"Flux, [Jy]")
plt.ylabel(r"Bias, [mas]")