import pandas as pd 
import numpy as np

def IndVal(data : pd.DataFrame): 
    specificity =  data.groupby(level = data.index.names[0]).mean().apply(lambda x: x / x.sum())
    fidelity = data.groupby(level = data.index.names[0]).apply(lambda x: x.astype(bool).sum() / len(x))
    return specificity * fidelity

def group_and_order_abundance_species(df : pd.DataFrame, groupby_columns : list[object], taxon_column: str, numeric_column : str): 
    df = df.groupby(groupby_columns + [taxon_column]).mean()
    df = df.groupby(level = groupby_columns, group_keys=False).apply(lambda x: x.sort_values([numeric_column], ascending=False))
    return df

def select_threshold(group, q : float ): 
    cumsum = group.cumsum() / group.sum()
    idx = np.sum(cumsum.to_numpy() < q) + 1
    return group.iloc[:idx]

def find_most_representative_species(df : pd.DataFrame, groupby_columns : list[object], taxon_column : str, numeric_column : str, threshold : float = 0.9, n_greates = None, relative_abundance : bool = True): 
    df_ordered = group_and_order_abundance_species(df, groupby_columns, taxon_column, numeric_column)
    
    if relative_abundance: 
        df_ordered = df_ordered.groupby(groupby_columns, group_keys=False).apply(lambda x: x / x.sum())
        
    if n_greates: 
        if threshold: 
            raise ValueError("threshold and n_greates can not be simultaneously specified!")
        return df_ordered.groupby(groupby_columns, group_keys=False).head(n_greates)
    elif threshold: 
        return df_ordered.groupby(groupby_columns, group_keys=False).apply(select_threshold, threshold)
    else: 
        raise ValueError("Either n_greates or threshold must be set!")

def find_outliers(df : pd.DataFrame, column : str = None, index : str = None):
    if column: 
        df_slice = df.loc[:, column]
    elif index: 
        df_slice = df.loc[index]
    else: 
        df_slice = df
        
    Q1 = df_slice.quantile(0.25)
    Q3 = df_slice.quantile(0.75)
    IQR = Q3 - Q1
    return df_slice[~df_slice.between(Q1 - 1.5 * IQR, Q3 + 1.5 * IQR, inclusive="both")].dropna()

def make_simplified_dataset(phyto_abundances : pd.DataFrame, MAX_DEPTH : float == 0.7): 
    phyto_abundances_simplified = phyto_abundances.loc[phyto_abundances["Sample_depth"] <= MAX_DEPTH, ["Region", "id", "Longitude", "Latitude", "Closest_coast", "Date", "Sample_depth", "Taxon", "Num_cell_l"]].copy()
    phyto_abundances_simplified["Genus"] = phyto_abundances_simplified["Taxon"].apply(lambda x: x.split(" ", maxsplit = 1)[0])
    phyto_abundances_simplified = phyto_abundances_simplified[["Region", "id", "Longitude",	"Latitude", "Closest_coast", "Date", "Sample_depth", "Genus", "Taxon", "Num_cell_l"]]
    #sometimes a reagion sampled more than once at the same depth
    return phyto_abundances_simplified.groupby(["Region", "id", "Longitude", "Latitude", "Date", "Sample_depth", "Genus", "Taxon"]).mean().reset_index()

def find_season(month, seasons): 
    for name, months in seasons.items(): 
        if month in months: 
            return name, month 

def make_string_season(dates, seasons): 
    season_list = np.empty(len(dates), dtype = "<U11")
    for i, date in enumerate(dates): 
        season, month = find_season(date.month, seasons) 
        season_list[i] = f"{season}-{date.year}"
    return season_list

def add_det_level_column(phyto_abund_simplified):
    phyto_abund_simplified["Det_level"] = "Species"
    phyto_abund_simplified.loc[phyto_abund_simplified["Genus"] == "Other", "Det_level"] = "Unknown"
    phyto_abund_simplified.loc[phyto_abund_simplified["Taxon"].str.contains("indet") & ~phyto_abund_simplified["Taxon"].str.contains("Other"), "Det_level"] = "Higher cat."
    phyto_abund_simplified.loc[phyto_abund_simplified["Taxon"].str.contains("spp."), "Det_level"] = "Genus"
    return phyto_abund_simplified

def add_season_column(phyto_abund_simplified, seasons): 
    phyto_abund_simplified["Season_year"] = make_string_season(phyto_abund_simplified["Date"], seasons)
    phyto_abund_simplified["Season"] = list(map(lambda x: x.split("-")[0], phyto_abund_simplified["Season_year"]))
    return phyto_abund_simplified

def add_coast_dist_column(df): 
    df["Coast_dist"] = "far"
    near_sta = df["Closest_coast"].between(3, 6.5, inclusive = "both")
    middle_sta = df["Closest_coast"].between(9, 12.5, inclusive = "both")
    df.loc[near_sta, "Coast_dist"] = "near"
    df.loc[middle_sta, "Coast_dist"] = "middle"
    return df