import os
import json
from example_simple import nested_dict


master_dir = '/home/ilya/github/bck/jetshow/pixel_all'
out_dir1 = '/home/ilya/github/bck/jetshow/pixel_25_33'
out_dir2 = '/home/ilya/github/bck/jetshow/pixel_19_24'
out_dir3 = '/home/ilya/github/bck/jetshow/pixel'

out_dirs = [out_dir1, out_dir2, out_dir3]

histories = nested_dict()
for out_dir in out_dirs:
    json_out = os.path.join(out_dir, 'history_mf.json')
    with open(json_out, 'r') as fo:
        history = json.load(fo)
        keys = history.keys()
        for key in keys:
            histories[key] = history[key]

json_out = os.path.join(master_dir, 'history_mf.json')
with open(json_out, 'w') as fo:
    json.dump(histories, fo)
