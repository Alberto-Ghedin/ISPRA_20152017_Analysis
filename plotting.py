import numpy as np
import seaborn as sns 
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram
import pandas as pd
from datetime import datetime
import os

def make_color_map(labels, palette_name : str = "Spectral"): 
    color_labels = np.unique(labels)
    color_values = sns.color_palette(palette_name, len(color_labels))
    return dict(zip(color_labels, color_values))

def make_colors_from_labels(labels, palette_name : str = "Spectral"):     
    color_map = make_color_map(labels, palette_name)
    return list(map(lambda i : color_map[i], labels))

def discrete_colormap(z : list, n_colors : int = None, palette_name : str = "jet"):
    
    z_min, z_max = min(z), max(z)
    if n_colors is None: 
        n_colors = len(np.unique(z))
    cm = plt.cm.get_cmap(palette_name, n_colors)
    bounds = np.linspace(z_min, z_max, n_colors)
    norm = mpl.colors.BoundaryNorm(bounds, n_colors)
    
    return cm, bounds, norm

def plot_italian_coast(ax, xs, ys, **kargs): 

    
    x_max = -np.inf
    x_min = np.inf
    y_max = -np.inf
    y_min = np.inf
    for x, y in zip(xs, ys): 
        x_max = max(x) if max(x) > x_max else x_max
        x_min = min(x) if min(x) < x_min else x_min
        y_max = max(y) if max(y) > y_max else y_max
        y_min = min(y) if min(y) < y_min else y_min
        ax.plot(x, y, 2, c = 'k', **kargs)
    ax.set_xlim(x_min * 0.99, x_max * 1.01)
    ax.set_ylim(y_min * 0.99, y_max * 1.01)

    return ax 


def compare_histograms(
    dfs,  
    variables : list[str], 
    sup_title : str=None, 
    titles : list[str]=None,
    legend_labels : list[str]=None,
    x_labels : list[str]=None,
    y_labels : list[str]="Relative frequency",
    n_bins : int = 50,
    **kwargs
): 
    if "figsize" in kwargs: 
        figsize_dims = kwargs["figsize"]
    else: 
        figsize_dims = (8,8) if len(variables) == 1 else (8,13)
        
    if len(variables) == 1: 
        fig, axs = plt.subplots(1,1, figsize=figsize_dims)
    else : 
        n_rows = len(variables) // 2 
        fig, axs = plt.subplots(n_rows,2, figsize=figsize_dims)
        axs = axs.flat
    
    if sup_title: 
        plt.suptitle(sup_title)
        
    for ax, var, title, x_label, y_label in zip(fig.axes, variables, titles, x_labels, y_labels):
        right_edge = -np.inf
        left_edge = np.inf
        datas = [df[var].dropna().to_numpy() for df in dfs]
        for data in datas:
            right_edge = data.max() if right_edge < data.max() else right_edge
            left_edge = data.min() if left_edge > data.min() else left_edge
            
        ax.set_title(title)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        
        for data, label in zip(datas, legend_labels): 
            ax.hist(data, 
                    bins=np.linspace(left_edge,right_edge, n_bins, endpoint=True), 
                    alpha=0.5, 
                    weights=np.ones_like(data) / len(data),
                    label=label)
            ax.legend()
        
        fig.tight_layout()


def plot_dendrogram(model, **kwargs):
    # Create linkage matrix and then plot the dendrogram

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

    linkage_matrix = np.column_stack(
        [model.children_, model.distances_, counts]
    ).astype(float)

    # Plot the corresponding dendrogram
    dendrogram(linkage_matrix, **kwargs)


# create directory to insert output of script
def make_sim_directory(description : str):
    def wrap(script):  
        def wrap_script(file_path : str, output_dir: str): 
            current_datetime = datetime.now()
            current_date = current_datetime.date().strftime("%Y_%m_%d")
            current_time = current_datetime.time().strftime("%H_%M_%S")
            output_dir += f"/{current_date}/{current_time}"
            os.makedirs(file_path + f"/{output_dir}", exist_ok=True)
            with open(file_path + f"/{output_dir}" + "/descr.txt", "a+") as description_file:
                description_file.write(description)
            script(file_path, output_dir)
        return wrap_script
    return wrap
