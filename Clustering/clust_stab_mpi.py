import datetime
import json
import numpy as np
import pandas as pd
import seaborn as sns
from mpi4py import MPI
from os import path
from pandas.tseries.offsets import MonthEnd
from reval.best_nclust_cv import FindBestClustCV
from reval.visualization import plot_metrics
from scipy import stats
from scipy.cluster.hierarchy import cophenet, fcluster
from scipy.spatial.distance import pdist
from sklearn import cluster, datasets, decomposition, metrics, mixture, preprocessing
from sklearn.cluster import AgglomerativeClustering, SpectralClustering
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier, kneighbors_graph
import time
import matplotlib.pyplot as plt

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
method = "spectral"

if rank == 0: 
    with open(path.expanduser("~") + "/sys_specific.json") as file: 
        params = json.load(file)
    _HOME_ = params["home"]

    sites_taxa = pd.read_csv(_HOME_ + "/PHD/ISPRA_20152017_Analysis/site_taxa_matrix.csv", index_col = [0,1])
    #hellinger transformation
    abund_hellinger = sites_taxa.apply(lambda x: np.sqrt(x / sum(x)), axis = 1)
    X = abund_hellinger.to_numpy()
    X_tr, X_ts, = train_test_split(X,
                                    test_size=0.25,
                                    random_state=42)
else: 
    X_tr = None
    X_ts = None

X_tr = comm.bcast(X_tr, root = 0)
X_ts = comm.bcast(X_ts, root = 0)



#CLUSTERING 
min_clusts = 3 
max_clusts = 7
iter_cv = 10

clust_list = list(range(min_clusts, max_clusts))

print(
    f"Rank {rank} commencing stability analysis of {clust_list[rank::size]}"
)
init = time.time()
classifier = KNeighborsClassifier()

if method == "spectral": 
    if rank == 0:
        hell_dist = euclidean_distances(X)
        sigma_2 = np.var(hell_dist)
    else: 
        sigma_2 = None
    
    sigma_2 = comm.bcast(sigma_2, root = 0)
    clustering = SpectralClustering(affinity = "rbf", gamma = 1 / (2 * sigma_2), assign_labels="cluster_qr")
else: 
    clustering = AgglomerativeClustering(metric = "euclidean", linkage = method)

findbestclust = FindBestClustCV(nfold=2,
                                nclust_range=clust_list[rank::size],
                                s=classifier,
                                c=clustering,
                                n_jobs=1,
                                nrand=100)

stab_analysis, nbest = findbestclust.best_nclust(X_tr, iter_cv=iter_cv)

print(f"Rank {rank} finished in {time.time() - init} seconds")

results = {}
results["mean"] = [value[1][0] for key, value in stab_analysis["val"].items()]
results["std"] = [value[1][1] for key, value in stab_analysis["val"].items()]

results = pd.DataFrame(results, index = clust_list[rank::size])
results = comm.gather(results, root = 0)

if rank == 0:
    results = pd.concat(results)
    results.to_csv(_HOME_ + f"/PHD/ISPRA_20152017_Analysis/Tests/{method}_stability_results_{min_clusts}_{max_clusts}.csv")
    print("end")
