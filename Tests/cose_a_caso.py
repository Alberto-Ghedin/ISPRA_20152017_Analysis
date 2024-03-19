## Correlation among species
#focused on Sample taken May 2017 (case where all regions sampled)
#Kendall tau and Spearman r produced 90% of pairs nto statistically significant
class Kendall_w(): 
    def __init__(self, a_posteriori : bool = True): 
        self.statistic = None
        self.p_value = None
        if a_posteriori: 
            self.spearman = None
            
    def compute_statistic(self, data: pd.DataFrame): 
        #each column is a judge
        n_objects, n_judges  = data.shape
        ranks = data.rank()
        r_i = ranks.sum(axis = 1)
        S = np.sum((r_i - np.mean(r_i))**2)
        #computing T 
        t_k = np.array([])
        for names, col in ranks.items(): 
            _, n_ties = np.unique(col.values, return_counts=True)
            t_k = np.append(
                t_k,
                n_ties[n_ties > 1]
            )
        T = np.sum(t_k**3 - t_k)
        W = 12 * S / (n_judges**2 * (n_objects**3 - n_objects) - n_judges * T)
        self.statistic = W
        F = (n_judges - 1) * W / (1 - W)
        d1 = n_objects - 1 - 2 / n_judges
        d2 = d1 * (n_judges - 1)
        self.p_value = 1 - stats.f.cdf(F, d1, d2)
    
    def spear2W(self, data : pd.DataFrame): 
        n_objects, n_judges  = data.shape
        r_mean = data.corr("spearman").where(
            np.triu(
                np.ones((n_judges,n_judges), dtype = bool), 
                k = 1
            )
        ).mean()
        return ((n_judges - 1) * r_mean + 1 ) / n_judges
    
    def compute_permutation_p_value(self, data: pd.DataFrame, col2perm: int, r_0 : float, n_perm : int = 999):
        n_objects, n_judges = data.shape
        spear_means = np.zeros(n_perm)
        df_shuffled = data.copy()
        for i in range(n_perm): 
            df_shuffled[col2perm] = np.random.permutation(df_shuffled[col2perm])
            spear_means[i] = np.delete(df_shuffled.corr("spearman")[col2perm], col2perm).mean()
        return np.sum(spear_means[spear_means > r_0]) / n_perm
    
    def a_posteriori(self, data: pd.DataFrame, n_perm : int): 
        n_objects, n_judges  = data.shape
        r_0 = data.corr("spearman").where(~np.identity(n_judges, dtype = bool)).mean().values
        print(r_0)
        p_values = np.vectorize(
            lambda col, r_0: self.compute_permutation_p_value(data, col, r_0, n_perm)
        )(
            range(n_judges),
            r_0
        )
        return pd.DataFrame({"p_vales" : p_values}, index = list(data.columns))
relative_abund = None
dates = phyto_abund_simplified["Date"]#.between(datetime.datetime(2017,5,1), datetime.datetime(2017,5,31))
for (region_name, df_region) in phyto_abund_simplified.loc[:, ["Region", "id", "Date", "Taxon", "Num_cell_l"]].groupby("Region"): 
    relative_abund_reg = find_most_representative_species(df_region.drop(columns="Region"), groupby_columns= ["id", "Date"], taxon_column="Taxon",  numeric_column = "Num_cell_l", threshold=1.1, relative_abundance=False)
    relative_abund_reg.reset_index(level = "Taxon", inplace=True)
    relative_abund_reg = relative_abund_reg[~relative_abund_reg["Taxon"].str.contains("Other")]
    relative_abund_reg = pd.pivot_table(relative_abund_reg, index= "id", columns = "Taxon", values="Num_cell_l").fillna(0).T
    relative_abund_reg.columns = pd.MultiIndex.from_tuples((region_name, station_id) for station_id in relative_abund_reg.columns)
    if relative_abund is None: 
        relative_abund = relative_abund_reg
    else: 
        relative_abund = pd.merge(left=relative_abund, right=relative_abund_reg, how = "outer", left_index=True, right_index=True)
relative_abund.fillna(0,inplace=True)
del relative_abund_reg
relative_abund = relative_abund.T
relative_abund.index.names = ["Region", "id"]
IndVal(relative_abund)
relative_abund = None
for (region_name, df_region) in phyto_abund_simplified.loc[:, ["Region", "id", "Date", "Taxon", "Num_cell_l"]].groupby("Region"): 
    relative_abund_reg = find_most_representative_species(df_region.drop(columns="Region"), groupby_columns= ["id", "Date"], taxon_column="Taxon",  numeric_column = "Num_cell_l", threshold=1.1, relative_abundance=False)
    relative_abund_reg.reset_index(level = "Taxon", inplace=True)
    relative_abund_reg = relative_abund_reg[~relative_abund_reg["Taxon"].str.contains("Other")]
    relative_abund_reg = pd.pivot_table(relative_abund_reg, index= "id", columns = "Taxon", values="Num_cell_l").fillna(0).T
    relative_abund_reg.columns = pd.MultiIndex.from_tuples((region_name, station_id) for station_id in relative_abund_reg.columns)
    if relative_abund is None: 
        relative_abund = relative_abund_reg
    else: 
        relative_abund = pd.merge(left=relative_abund, right=relative_abund_reg, how = "outer", left_index=True, right_index=True)
relative_abund.fillna(0,inplace=True)
del relative_abund_reg
df = pd.pivot_table(phyto_abund_simplified.loc[:,["Region", "id", "Date", "Taxon", "Num_cell_l"]], index = ["Region", "id", "Date"], columns = "Taxon", values = "Num_cell_l").fillna(0)
df.reset_index().to_csv("/mnt/d/PHD/real_abund_phyto_tot.csv", index = False)
IndVal(relative_abund.T)
relative_abund.astype(bool).sum(axis=1).describe()
relative_abund = relative_abund.T.apply(lambda x: np.sqrt(x / sum(x)), axis = 1)
from scipy.stats import spearmanr
from scipy import stats
values = relative_abund.to_numpy().T
N, _ = values.shape
index_p_values = {}
for i in range(N): 
    for j in range(i+1,N):
        p_val = stats.kendalltau(values[i, :], values[j,:], method="asymptotic").pvalue
        if p_val > 0.05: 
            index_p_values[(i,j)] = p_val
2 * 75529 / (410 * (410-1))
len(index_p_values)
p_values = relative_abund.corr(method= lambda a,b: stats.kendalltau(a,b, method="asymptotic").pvalue)
len(p_values.to_numpy().flatten())
410**2
sum(p_values.to_numpy().flatten() < 0.05)
X = relative_abund.apply(lambda x: (x - np.mean(x)) / np.std(x), axis = 1).to_numpy()

X = relative_abund.apply(lambda x: (x - np.mean(x)) / np.std(x), axis = 0).to_numpy()
pca = decomposition.PCA(n_components = 2)
pca.fit(X)
df = pd.DataFrame(data = pca.transform(X)[:, [0,1]], columns=["PC1", "PC2"])
H_clusters = AgglomerativeClustering(n_clusters = 4, metric = "euclidean", linkage = "ward", compute_distances = True)
clusters = H_clusters.fit(X)
fig, ax = plt.subplots(1,1, figsize=(20,10))
plot_dendrogram(clusters, ax = ax)
plt.xlabel("Number of points in node (or index of point if no parenthesis).")
plt.show()
pca.explained_variance_
df["clusters"] = clusters.labels_
fig, ax = plt.subplots(1,1, figsize=(15,10))
sns.scatterplot(df.loc[:, :], x = "PC1", y = "PC2", s = 120, hue = "clusters", palette = "deep")
spear = relative_abund.corr("spearman")
spear
fig, ax = plt.subplots(1,1, figsize=(13, 11))
sns.heatmap(spear, ax = ax)
H_clusters = AgglomerativeClustering(n_clusters = 4, metric = "precomputed", linkage = "average", compute_distances = True)
clusters = H_clusters.fit((1 - spear)**2)
fig, ax = plt.subplots(1,1, figsize=(20,10))
plot_dendrogram(clusters, ax = ax , truncate_mode = "level", p = 3
               )
plt.xlabel("Number of points in node (or index of point if no parenthesis).")
plt.show()
X = ((1 - spear)**2).to_numpy()
pca = decomposition.KernelPCA(n_components = 200, kernel="precomputed")
pca.fit(X)
df = pd.DataFrame(data = pca.transform(X)[:, [0,1]], columns=["PC1", "PC2"])
pca.eigenvalues_[[1,2]] / sum(pca.eigenvalues_)
df["cluster"] = clusters.labels_
fig, ax = plt.subplots(1,1, figsize=(15,10))
sns.scatterplot(df.loc[:, :], x = "PC1", y = "PC2", s = 120, hue = "cluster", palette = "deep")
relative_abund = None
dates = phyto_abund_simplified["Date"].between(datetime.datetime(2017,5,1), datetime.datetime(2017,5,31))
for (region_name, df_region) in phyto_abund_simplified.loc[dates, ["Region", "id", "Date", "Taxon", "Num_cell_l"]].groupby("Region"): 
    relative_abund_reg = find_most_representative_species(df_region.drop(columns="Region"), groupby_columns= ["id", "Date"], taxon_column="Taxon",  numeric_column = "Num_cell_l", threshold=1.1)
    relative_abund_reg.reset_index(level = "Taxon", inplace=True)
    relative_abund_reg = relative_abund_reg[~relative_abund_reg["Taxon"].str.contains("Other")]
    relative_abund_reg = pd.pivot_table(relative_abund_reg, index= "id", columns = "Taxon", values="Num_cell_l").fillna(0).T
    relative_abund_reg.columns = pd.MultiIndex.from_tuples((region_name, station_id) for station_id in relative_abund_reg.columns)
    if relative_abund is None: 
        relative_abund = relative_abund_reg
    else: 
        relative_abund = pd.merge(left=relative_abund, right=relative_abund_reg, how = "outer", left_index=True, right_index=True)
relative_abund.fillna(0,inplace=True)
del relative_abund_reg
### abundances
phyto_abundances_per_species = phyto_abund_simplified[["Taxon", "Num_cell_l"]].groupby("Taxon").median().sort_values("Num_cell_l", ascending=False).reset_index()
n_species = 20
abundance = phyto_abundances_per_species.loc[0:n_species-1, "Num_cell_l"].copy().to_numpy()  
species = phyto_abundances_per_species.loc[0:n_species-1, "Taxon"].copy().to_numpy()
cm, bounds, norm = discrete_colormap(np.arange(n_species),n_species)
fig, ax = plt.subplots(1,1, figsize=(9,6))
ax.set_yscale("log")
for i, (spec, taxon) in enumerate(zip(abundance, species)): 
    ax.bar(i + 1, height = spec, label = taxon, color = cm(i))
ax.legend(loc=(1.04, 0))
ax.set_ylabel("Abundance [cell / l]")
ax.set_xlabel("")
ax.set_title("Abundance of species across Italy, surface level")
ax.set_xticks([i for i in range(1, n_species, 2)]);
#plt.savefig("/mnt/d/PHD/Most_abundant_species.png")
df_most_abund_species = phyto_abundances_simplified.loc[phyto_abundances_simplified["Taxon"].isin(most_abund_species), ["id", "Longitude", "Latitude", "Taxon", "Num_cell_l"]].groupby(["id", "Taxon"]).agg({"Longitude" : "first", "Latitude" : "first", "Num_cell_l" : "mean"}).reset_index()
df_most_abund_species["Num_cell_log"] = df_most_abund_species["Num_cell_l"].apply(np.log10)
#df_most_abund_species = df_most_abund_species.query("Num_cell_l >= 3.0")
n_rows = len(most_abund_species) // 2 
figsize_dims = (12,10)
fig, axs = plt.subplots(n_rows,2, figsize=figsize_dims);
axs = axs.flat[:len(most_abund_species)+1]

for ax, (name, df) in zip(axs, df_most_abund_species.groupby(by = "Taxon")): 
    cm, bounds, norm = discrete_colormap(df["Num_cell_log"].to_numpy(), 15);
    plot_italian_coast(ax, (x_pen, x_sic, x_sard), (y_pen, y_sic, y_sard), alpha = 0.7);
    sc = ax.scatter(*df[["Longitude", "Latitude"]].to_numpy().T, c = df["Num_cell_log"].to_numpy(), cmap = cm, norm = norm, s = 20);
    ax.set_xlim(df["Longitude"].min() * 0.99, df["Longitude"].max() * 1.01)
    ax.set_ylim(df["Latitude"].min() * 0.99, df["Latitude"].max() * 1.01)
    ax.set_xlabel("Longitude");
    ax.set_ylabel("Latitude");
    ax.set_title(name);
    plt.colorbar(sc, cmap=cm, norm=norm,
        spacing='proportional', ticks=bounds, boundaries=bounds, format='%2.2f', label = "log((# individuals) / l)");
plt.tight_layout()
plt.show()
phyto_abundances_per_genus = phyto_abundances_simplified[["Genus", "Num_cell_l"]].groupby("Genus").median().sort_values("Num_cell_l", ascending=False).reset_index()
n_species = 20
abundance = phyto_abundances_per_genus.loc[0:n_species-1, "Num_cell_l"].copy().to_numpy()  
order_of_magnitude = floor(abs(np.log10(abundance.max())))
abundance /= 10 ** order_of_magnitude
species = phyto_abundances_per_genus.loc[0:n_species-1, "Genus"].copy().to_numpy()
cm, bounds, norm = discrete_colormap(np.arange(n_species),n_species)
fig, ax = plt.subplots(1,1, figsize=(9,6))
for i, (spec, taxon) in enumerate(zip(abundance, species)): 
    ax.bar(i + 1, height = spec, label = taxon, color = cm(i))
ax.legend(loc=(1.04, 0))
ax.set_ylabel(r"Abundance (ind. / l) $\times 10^{}$".format(order_of_magnitude))
ax.set_xlabel("")
ax.set_title("Abundance of species across Italy (averages over time and " + r"depth$ = 0 \; m$)")
ax.set_xticks([i for i in range(1, n_species, 2)]);
#plt.savefig("/mnt/d/PHD/Most_abundant_species.png")
perform PCA on IndVal
X = IndVal_italy.to_numpy().T
X = np.sqrt(X / np.sum(X, axis = 1)[:, None])
kpca = decomposition.KernelPCA(n_components=2)
kpca.fit(X)
df = pd.DataFrame(data = kpca.transform(X), columns=["PC1", "PC2"])
df["Region"] = IndVal_italy.columns.get_level_values(0)

fig, ax = plt.subplots(1,1, figsize=(15,10))
sns.scatterplot(df, x = "PC1", y = "PC2", hue = "Region", palette = "deep", ax = ax, style = "Region", s = 120)
Hireranchical clustering
X = IndVal_italy.to_numpy().T
X = np.sqrt(X / np.sum(X, axis = 1)[:, None])
model = AgglomerativeClustering(linkage="ward", compute_distances=True, n_clusters = 8).fit(X)
fig, ax = plt.subplots(1,1, figsize=(20,10))
plot_dendrogram(model, ax = ax)
plt.xlabel("Number of points in node (or index of point if no parenthesis).")
plt.show()
fig, ax = plt.subplots(1,1, figsize=(15,10))
colors = make_colors_from_labels(model.labels_)
cm, bounds, norm = discrete_colormap(model.labels_)
plot_italian_coast(ax, (x_pen, x_sic, x_sard), (y_pen, y_sic, y_sard), alpha = 0.5)
sc = ax.scatter(*df_lat_long_stations.set_index("id").loc[IndVal_italy.columns.get_level_values(1).to_list(), ["Longitude", "Latitude"]].to_numpy().T, c = model.labels_, cmap = cm, norm = norm)
ax.set_xlabel("Longitude")
ax.set_ylabel("Latitude")
plt.colorbar(sc, cmap=cm, norm=norm,
    spacing='proportional', ticks=bounds, boundaries=bounds, format='%1i', label = "# of species")
plt.show()