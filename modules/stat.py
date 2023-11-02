import numpy as np
import pandas as pd

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
    

def IndVal(data : pd.DataFrame): 
    specificity =  data.groupby(level = data.index.names[0]).mean().apply(lambda x: x / x.sum())
    fidelity = data.groupby(level = data.index.names[0]).apply(lambda x: x.astype(bool).sum() / len(x))
    return specificity * fidelity


def Correspondance_analysis(data : np.ndarray, scaling : int = 1, n_rows : int = None, n_cols : int = None, n_axis : int = 2) -> np.ndarray: 
    """
    Compute correspondace analysis of a non-normalizied matrix.
    (row vectors, columns vectors)
    n_rows & n_cols are the number of rows and column to return 
    if scaling == 1: data rows are at the centroids, if scaling == 2 otherwise
    """
    Q = data / data.sum()
    weight_row = Q.sum(axis=1)
    weight_col = Q.sum(axis=0)
    Q_bar = (Q- np.outer(weight_row, weight_col)) / np.sqrt(np.outer(weight_row, weight_col))

    if scaling == 1:
        _, eigenvects = np.linalg.eigh(Q_bar.T @ Q_bar)
        if n_cols is None: 
            _, n_cols = Q.shape
        if n_rows is None: 
            n_rows, _ = Q.shape
    
        eigenvects = np.flip(eigenvects, axis = 1)[:n_cols,:n_axis]
        V = np.diag(1 / np.sqrt(weight_col[:n_cols])) @ eigenvects 
        F = np.diag(1 / weight_row[:n_rows]) @ Q[:n_rows, :n_cols] @ V
        return F, V 

    elif scaling == 2: 
        _, eigenvects = np.linalg.eigh(Q_bar @ Q_bar.T)
        if n_cols is None: 
            _, n_cols = Q.shape
        if n_rows is None: 
            n_rows, _ = Q.shape
        eigenvects = np.flip(eigenvects, axis = 1)[:n_rows,:n_axis]
        V = np.diag(1 / np.sqrt(weight_row[:n_rows])) @ eigenvects 
        F = np.diag(1 / weight_col[:n_cols]) @ Q.T[:n_cols, :n_rows] @ V
        return F, V
    
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