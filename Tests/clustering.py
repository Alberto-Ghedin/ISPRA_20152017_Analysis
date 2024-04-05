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

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

with open(path.expanduser("~") + "/sys_specific.json") as file: 
    params = json.load(file)
_HOME_ = params["home"]
from plotting import * 

sites_taxa = pd.read_csv(_HOME_ + "/PHD/ISPRA_20152017_Analysis/site_taxa_matrix.csv", index_col = 0)

## Clustering functions utility
def build_linkage_matrix(model) -> np.ndarray : 
    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    return np.column_stack(
        [model.children_, model.distances_, counts]
    ).astype(float)

## Taxa subset
with open(_HOME_ + "/PHD/ISPRA_20152017_Analysis/selected_species.txt", 'r') as file:
    relevant_species = [line.strip() for line in file]
sites_taxa = sites_taxa.loc[:, relevant_species]

#hellinger transformation
abund_hellinger = sites_taxa.apply(lambda x: np.sqrt(x / sum(x)), axis = 1)
X = abund_hellinger.to_numpy()

#CLUSTERING 
min_clusts = 3
max_clusts = 20

silhouette_scores = {}
calisnki_scores = {}
davies_scores = {}
method_labels = {}

n_clusters = np.arange(min_clusts,max_clusts,1)
n_times = len(n_clusters[rank::size])
silhouette = np.empty(n_times)
calisnki = np.empty(n_times)
davies = np.empty(n_times)




#SINGLE LINKAGE CLUSTERING
method = "single"

if rank == 0:
    print(f"Computing {method} clustering")
if rank == 0: 
    single_linkage_clustering = AgglomerativeClustering(n_clusters = None, metric = "euclidean", linkage = "single", compute_distances = True, compute_full_tree=True, distance_threshold=0.0)
    linkage_clusters = single_linkage_clustering.fit(X)
    linkage_matrix = build_linkage_matrix(linkage_clusters)
    del single_linkage_clustering, linkage_clusters
else: 
    linkage_matrix = None
linkage_matrix = comm.bcast(linkage_matrix, root = 0)


single_link_labels = np.empty((n_times, X.shape[0]))
for i, n in enumerate(n_clusters[rank::size]):
    single_link_labels[i,:] = fcluster(linkage_matrix, t = n, criterion = "maxclust")
    silhouette[i] = metrics.silhouette_score(X, single_link_labels[i,:])
    calisnki[i] = metrics.calinski_harabasz_score(X, single_link_labels[i,:])
    davies[i] = metrics.davies_bouldin_score(X, single_link_labels[i,:])

silhouette_scores[method] = silhouette
calisnki_scores[method] = calisnki
davies_scores[method] = davies
method_labels[method] = single_link_labels
del linkage_matrix

#COMPLETE LINKAGE CLUSTERING
method = "complete"
if rank == 0:
    print(f"Computing {method} clustering")

complete_link_lables = np.empty((n_times, X.shape[0]))
for i, n in enumerate(n_clusters[rank::size]):
    complete_linkage_clustering = AgglomerativeClustering(n_clusters = n, metric = "euclidean", linkage = method, compute_distances = True, compute_full_tree=True)
    linkage_clusters = complete_linkage_clustering.fit(X)
    complete_link_lables[i,:] = linkage_clusters.labels_
    silhouette[i] = metrics.silhouette_score(X, complete_link_lables[i,:])
    calisnki[i] = metrics.calinski_harabasz_score(X, complete_link_lables[i,:])
    davies[i] = metrics.davies_bouldin_score(X, complete_link_lables[i,:])

silhouette_scores[method] = silhouette
calisnki_scores[method] = calisnki
davies_scores[method] = davies
method_labels[method] = complete_link_lables
del complete_linkage_clustering, linkage_clusters

#WARD LINKAGE CLUSTERING
method = "ward"
if rank == 0:
    print(f"Computing {method} clustering")

if rank == 0:
    ward_clustering = AgglomerativeClustering(n_clusters = None, metric = "euclidean", linkage = "ward", compute_distances = True, compute_full_tree=True, distance_threshold=0.0)
    ward_clusters = ward_clustering.fit(X)
    linkage_matrix = build_linkage_matrix(ward_clusters)
    del ward_clustering, ward_clusters
else:
    linkage_matrix = None
linkage_matrix = comm.bcast(linkage_matrix, root = 0)

ward_labels = np.empty((n_times, X.shape[0]))
for i, n in enumerate(n_clusters[rank::size]):
    ward_labels[i,:] = fcluster(linkage_matrix, t = n, criterion = "maxclust")
    silhouette[i] = metrics.silhouette_score(X, ward_labels[i,:])
    calisnki[i] = metrics.calinski_harabasz_score(X, ward_labels[i,:])
    davies[i] = metrics.davies_bouldin_score(X, ward_labels[i,:])

silhouette_scores[method] = silhouette
calisnki_scores[method] = calisnki
davies_scores[method] = davies
method_labels[method] = ward_labels
del linkage_matrix

#SPECTRAL CLUSTERING correlation
method = "spectral_corr"
if rank == 0:
    print(f"Computing {method} clustering")
if rank == 0:
    eco_matrix = sites_taxa.to_numpy()
    D_c = np.diag(np.sum(eco_matrix, axis = 0))
    #compute diagonal inverse of D_c 
    inv = np.diag(1 / np.diag(D_c))
    similarity = eco_matrix @ inv @ eco_matrix.T
    del D_c, inv
else: 
    similarity = None
    eco_matrix = None
    
similarity = comm.bcast(similarity, root = 0)
eco_matrix = comm.bcast(eco_matrix, root = 0)
spectral_labels = np.empty((n_times, X.shape[0]))
for i, n in enumerate(n_clusters[rank::size]):
    spectral_clustering = SpectralClustering(n_clusters = n, affinity = "precomputed", assign_labels="cluster_qr")
    spectral_labels[i,:] = spectral_clustering.fit_predict(similarity)
    silhouette[i] = metrics.silhouette_score(eco_matrix, spectral_labels[i,:])
    calisnki[i] = metrics.calinski_harabasz_score(eco_matrix, spectral_labels[i,:])
    davies[i] = metrics.davies_bouldin_score(eco_matrix, spectral_labels[i,:])

silhouette_scores[method] = silhouette
calisnki_scores[method] = calisnki
davies_scores[method] = davies
method_labels[method] = spectral_labels

#SPECTRAL CLUSTERING euclidean
method = "spectral_dist"
if rank == 0:
    print(f"Computing {method} clustering")
if rank == 0:
    hell_dist = euclidean_distances(X)
    delta_2 = np.var(hell_dist)
    similarity = np.exp(-hell_dist**2 / (2 * delta_2))
    del hell_dist, delta_2
else: 
    similarity = None

similarity = comm.bcast(similarity, root = 0)


for i, n in enumerate(n_clusters[rank::size]):
    spectral_clustering = SpectralClustering(n_clusters = n, affinity = "precomputed", assign_labels="cluster_qr")
    spectral_labels[i,:] = spectral_clustering.fit_predict(similarity)
    silhouette[i] = metrics.silhouette_score(X, spectral_labels[i,:])
    calisnki[i] = metrics.calinski_harabasz_score(X, spectral_labels[i,:])
    davies[i] = metrics.davies_bouldin_score(X, spectral_labels[i,:])
silhouette_scores[method] = silhouette
calisnki_scores[method] = calisnki
davies_scores[method] = davies
method_labels[method] = spectral_labels


if rank == 0:
    print(f"Computing diagnostics")

# number of items per cluster
values = {}
for method, labels in method_labels.items(): 
    values[method] = []
    for i in range(labels.shape[0]): 
        values[method] += [np.unique(labels[i,:], return_counts = True)[1]]
    values[method] = np.concatenate(values[method])
df_n_cluster_members = pd.DataFrame(
    values,
    index = pd.MultiIndex.from_tuples([(n_tot, ith) for n_tot in n_clusters[rank::size] for ith in range(n_tot)], names = ["N_clusters", "Index"])
)

df_n_cluster_members = comm.gather(df_n_cluster_members, root = 0)

if rank == 0: 
    df_n_cluster_members = pd.concat(df_n_cluster_members, axis = 0)
    df_n_cluster_members.sort_index().to_csv(_HOME_ + "/PHD/ISPRA_20152017_Analysis/n_cluster_members.csv")
 
#consensun amog different methods
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
    df_consensus[("ARI", method_1 + "_" + method_2)] = scores_1
    df_consensus[("NMI", method_1 + "_" + method_2)] = scores_2
df_consensus = df_consensus.loc[:, ["ARI", "NMI"]]
df_consensus.index = n_clusters[rank::size]

df_consensus = comm.gather(df_consensus, root = 0)
if rank == 0:
    df_consensus = pd.concat(df_consensus, axis = 0)
    df_consensus.sort_index().to_csv(_HOME_ + "/PHD/ISPRA_20152017_Analysis/clust_consensus.csv")

# #pooling everything together
# silhouette_scores = comm.gather(silhouette_scores, root = 0)
# calisnki_scores = comm.gather(calisnki_scores, root = 0)
# davies_scores = comm.gather(davies_scores, root = 0)
# method_labels = comm.gather(method_labels, root = 0)
# index_used = comm.gather(n_clusters[rank::size], root = 0)

# if rank == 0: 
#     index_used = np.concatenate(index_used, axis = None)
#     scores = {}
#     pooled_dict = {}
#     for (score, dics) in zip([silhouette_scores, calisnki_scores, davies_scores], ["silhouette", "calinski", "davies"]): 
#         pooled_dict = {}
#         for dic in score: 
#             for method in dic: 
#                 if method not in pooled_dict: 
#                     pooled_dict[method] = dic[method]
#                 else:
#                     pooled_dict[method] = np.concatenate([pooled_dict[method], dic[method]], axis = None)
#         scores[dics] = pooled_dict

#     cluster_labels = {}
#     for dic in method_labels: 
#         for method in dic: 
#             if method not in cluster_labels: 
#                 cluster_labels[method] = dic[method]
#             else:
#                 cluster_labels[method] = np.concatenate([cluster_labels[method], dic[method]], axis = 0)
    
#     print(
#          f"{index_used}\n"
#          f"Silhouette scores: {scores['silhouette']}\n"
#          f"Calinski scores: {scores['calinski']}\n"
#          f"Davies scores: {scores['davies']}\n"
#          f"Labels: {cluster_labels}"
#     )

