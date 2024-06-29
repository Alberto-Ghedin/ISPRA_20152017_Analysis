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
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier, kneighbors_graph
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def compute(train : np.ndarray[np.float64], 
         #test : np.array[float],
         min_clusts : int, 
         max_clusts : int, 
         method : str, 
         output_dir : str, 
         iter_cv : int = 10, 
         nfolds : int = 2,
         classifier = KNeighborsClassifier(), 
         nrands : int = 100
         ):


    clust_list = list(range(min_clusts, max_clusts))

    print(
    f"Rank {rank} commencing stability analysis of {clust_list[rank::size]}"
    )
    init = time.time()
    
    clustering = AgglomerativeClustering(metric = "euclidean", linkage = method)


    findbestclust = FindBestClustCV(nfold=nfolds,
                                    nclust_range=clust_list[rank::size],
                                    s=classifier,
                                    c=clustering,
                                    n_jobs=1,
                                    nrand=nrands)

    stab_analysis, nbest = findbestclust.best_nclust(train, iter_cv=iter_cv)

    print(f"Rank {rank} finished in {time.time() - init} seconds")

    results = {}
    results["mean"] = [value[1][0] for key, value in stab_analysis["val"].items()]
    results["std"] = [value[1][1] for key, value in stab_analysis["val"].items()]

    results = pd.DataFrame(results, index = clust_list[rank::size])
    results = comm.gather(results, root = 0)

    if rank == 0:
        results = pd.concat(results)
        results.to_csv(f"{output_dir}/{method}_stability_results_{min_clusts}_{max_clusts}.csv")
        print("end")
