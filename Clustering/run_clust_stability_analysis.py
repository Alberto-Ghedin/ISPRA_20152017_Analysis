import argparse
from os import path 
import json
import pandas as pd 
from sklearn.model_selection import train_test_split
import numpy as np
from plotting import create_sim_subdirectories
import stability_spectral, stability_hierarchical
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

if rank == 0: 
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
    path_to_relevant_species = _HOME_ + params.get("path_to_relevant_species", "")
    min_clusts = params["min_clusts"]
    max_clusts = params["max_clusts"]
    sigmas = params.get("sigmas", None)
    nrands = params["nrands"]

    metadata = {"title" : "stability analysis", 
                "additional_info" : "stability analysis of spectral clustering on hellinger transformed data, only relevant species included",
                "method" : method,
                } | (params)


    output_dir = create_sim_subdirectories(path_to_output, output_dir, metadata = metadata)
    sites_taxa = pd.read_csv(path_to_data, index_col=[0, 1])

    with open(path_to_relevant_species, 'r') as file:
        relevant_species = [line.strip() for line in file]

    #hellinger transformation
    abund_hellinger = sites_taxa.loc[:, relevant_species].apply(lambda x: np.sqrt(x / sum(x)), axis=1)
    X = abund_hellinger.to_numpy()
    X_tr, X_ts, = train_test_split(X,
                                    test_size=0.25,
                                    random_state=42)
else: 
    output_dir = None
    X_tr = None
    min_clusts = None
    max_clusts = None
    sigmas = None
    nrands = None
    method = None

output_dir = comm.bcast(output_dir, root = 0)
X_tr = comm.bcast(X_tr, root = 0)
min_clusts = comm.bcast(min_clusts, root = 0)
max_clusts = comm.bcast(max_clusts, root = 0)
sigmas = comm.bcast(sigmas, root = 0)
nrands = comm.bcast(nrands, root = 0)
method = comm.bcast(method, root = 0)



#SELECTING METHOD
if method == "spectral": 
    stability_spectral.compute(X_tr, min_clusts, max_clusts, sigmas, output_dir, nrands = nrands)
elif method in ["ward", "complete"]: 
     stability_hierarchical.compute(X_tr, min_clusts, max_clusts, method, output_dir, nrands = nrands)
else: 
    raise ValueError(f"Method {method} not recognized")
