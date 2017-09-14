import os
from example_simple import plot_stripe_for_given_key


out_dir = '/home/ilya/github/bck/jetshow/uvf_mf_adds'
savefig = os.path.join(out_dir, 'bk_3d_1.7_8.1.png')
key_low = '09_1.665'
key_high = '09_8.1'
simulated_image_low = os.path.join(out_dir, "map_i_{}.txt".format(key_low))
simulated_image_high = os.path.join(out_dir, "map_i_{}.txt".format(key_high))
difmap_model_low = os.path.join(out_dir, "bk_{}.mdl".format(key_low))
difmap_model_high = os.path.join(out_dir, "bk_{}.mdl".format(key_high))
json_history = os.path.join(out_dir, 'history_mf.json')

fig = plot_stripe_for_given_key(simulated_image_low, simulated_image_high,
                                json_history, key_low, key_high,
                                difmap_model_low, difmap_model_high,
                                savefig=savefig, delta=200)
