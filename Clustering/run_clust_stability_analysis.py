import argparse
from os import path 
import json
import pandas as pd 
from sklearn.model_selection import train_test_split
import numpy as np
from plotting import create_sim_subdirectories
import stability_spectral

with open(path.expanduser("~") + "/sys_specific.json") as file: 
        params = json.load(file)
_HOME_ = params["home"]

parser = argparse.ArgumentParser()
parser.add_argument("params", help="Path to the params json file")
args = parser.parse_args()
params_path = args.params
with open(params_path, "r") as file:
    params = json.load(file)
method = params["method"]
output_dir = params["output_dir"]
path_to_output = _HOME_ + params["path_to_output"]
path_to_data = _HOME_ + params["path_to_data"]
min_clusts = params["min_clusts"]
max_clusts = params["max_clusts"]
sigmas = params["sigmas"]
nrands = params["nrands"]

description = f"stability analysis"

output_dir = create_sim_subdirectories(path_to_output, output_dir, description)
sites_taxa = pd.read_csv(path_to_data, index_col = [0,1])

#hellinger transformation
abund_hellinger = sites_taxa.apply(lambda x: np.sqrt(x / sum(x)), axis = 1)
X = abund_hellinger.to_numpy()
X_tr, X_ts, = train_test_split(X,
                                test_size=0.25,
                                random_state=42)


#SELECTING METHOD
if method == "spectral": 
    stability_spectral.compute(X_tr, min_clusts, max_clusts, sigmas, output_dir, nrands = nrands)
elif method in ["ward", "complete"]: 
    pass 
