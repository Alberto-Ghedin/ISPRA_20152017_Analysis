import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering
from tqdm.auto import tqdm
from pandas.tseries.offsets import MonthEnd
import datetime
from pandas.tseries.offsets import MonthEnd
import json
from sklearn import metrics
from scipy.cluster.hierarchy import fcluster
from sklearn.cluster import SpectralClustering
from itertools import combinations
from numba import njit, prange
from os import path
from mpi4py import MPI
from sklearn.metrics.pairwise import euclidean_distances
import argparse
from plotting import make_sim_directory, create_sim_subdirectories


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

with open(path.expanduser("~") + "/sys_specific.json") as file: 
    params = json.load(file)
_HOME_ = params["home"]


if rank == 0: 
    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("params", help="Path to the params json file")
    args = parser.parse_args()
    params_path = args.params
    with open(params_path, "r") as file:
        params = json.load(file)
    output_dir = params["output_dir"]
    path_to_outout = _HOME_ + params["path_to_output"]
    path_to_data = _HOME_ + params["path_to_data"]
    min_clusts = params["min_clusts"]
    max_clusts = params["max_clusts"]
    sigmas = params["sigmas"]
   
    description = f"Consensus among spectral clustering with varying sigma \n \
                    min_clusts = {min_clusts} \n \
                    max_clusts = {max_clusts} \n \
                    sigmas = {sigmas} \n \
                    "

    output_dir = create_sim_subdirectories(path_to_outout, output_dir, description)
    sites_taxa = pd.read_csv(path_to_data, index_col = [0, 1])
    #hellinger transformation
    abund_hellinger = sites_taxa.apply(lambda x: np.sqrt(x / sum(x)), axis = 1)
    X = abund_hellinger.to_numpy()
    hell_dist = euclidean_distances(X)
    n_clusters = np.arange(min_clusts,max_clusts,1)

else: 
    output_dir = None
    hell_dist = None
    n_clusters = None
    sigmas = None

output_dir = comm.bcast(output_dir, root = 0)
hell_dist = comm.bcast(hell_dist, root = 0)
n_clusters = comm.bcast(n_clusters, root = 0)
sigmas = comm.bcast(sigmas, root = 0)

method_labels = {}

n_times = len(n_clusters)

method = "spectral"

for sigma in sigmas[rank::size]:
    similarity = np.exp(-hell_dist**2 / (2 * sigma))
    spectral_labels = np.empty((n_times, hell_dist.shape[0]))
    for i, n in enumerate(n_clusters):
        #print(f"Rank {rank} commencing spectral clustering with sigma = {sigma} and n_clusters = {n}")
        spectral_clustering = SpectralClustering(n_clusters = n, affinity = "precomputed", assign_labels="cluster_qr")
        spectral_labels[i,:] = spectral_clustering.fit_predict(similarity)
    method_labels[sigma] = spectral_labels

method_labels = comm.gather(method_labels, root = 0)

if rank == 0:
    print(f"Computing diagnostics")
    
    method_labels = {k: v for d in method_labels for k, v in d.items()}
    #consensun among different methods
    df_consensus = pd.DataFrame(
        {},
        columns = pd.MultiIndex(levels = [[], []], codes = [[], []], names = ["Index", "Methods"])
    )


    for method_1, method_2 in combinations(method_labels.keys(), 2): 
        labels_1 = method_labels[method_1]
        labels_2 = method_labels[method_2]
        scores_1 = np.zeros(labels_1.shape[0])
        scores_2 = np.zeros(labels_2.shape[0])
        for i in range(labels_1.shape[0]): 
            scores_1[i] = metrics.adjusted_rand_score(labels_1[i,:], labels_2[i,:])
            scores_2[i] = metrics.normalized_mutual_info_score(labels_1[i,:], labels_2[i,:])
        df_consensus[("ARI", f"{method_1}_{method_2}")] = scores_1
        df_consensus[("NMI", f"{method_1}_{method_2}")] = scores_2
    df_consensus = df_consensus.loc[:, ["ARI", "NMI"]]
    df_consensus.index = n_clusters
    df_consensus.sort_index().to_csv(f"{output_dir}/clust_consensus_spectral_sigma_{min_clusts}_{max_clusts}.csv")