### EXAMPLE OF PLOTS AND MAPS WITH GEOPANDAS AND MATPLOTLIB


abund_basin = phyto_abund_simplified.groupby("Basin", observed=True).sum(numeric_only=True)[["Num_cell_l"]]
abund_taxon_basin = phyto_abund_simplified.groupby(["Basin", "Taxon"], observed=True).sum(numeric_only=True)[["Num_cell_l"]].sort_values("Num_cell_l")
top_abund_taxa_basin = abund_taxon_basin.div(abund_basin, axis=0).reset_index().sort_values(["Basin", "Num_cell_l"], ascending=[True, False]).groupby("Basin", observed=True, sort=False).head(5)
fig, ax = plt.subplots(figsize=(15, 15))

# Plot the map of Italy
italy.plot(ax=ax, color='lightgrey')

# Define the positions for each pie chart
# Adjust these values based on the number of basins and desired layout
positions = [
    [0.1, 0.6, 0.15, 0.15],  # [left, bottom, width, height]
    [0.3, 0.6, 0.15, 0.15],
    [0.5, 0.6, 0.15, 0.15],
    [0.7, 0.6, 0.15, 0.15],
    [0.1, 0.3, 0.15, 0.15],
    [0.3, 0.3, 0.15, 0.15],
    [0.5, 0.3, 0.15, 0.15],
    [0.7, 0.3, 0.15, 0.15]
]

# Plot each pie chart in its respective position
for pos, (basin, group) in zip(positions, top_abund_taxa_basin.groupby('Basin', observed=True)):
    inset_ax = fig.add_axes(pos)
    inset_ax.pie(group['Num_cell_l'], labels=group['Taxon'], autopct='%1.1f%%')
    inset_ax.set_title(f'Basin: {basin}', fontsize=10)

# Adjust layout
plt.tight_layout()
plt.show()




fig, ax = plt.subplots(figsize=(15, 15))

# Plot the map of Italy
italy.plot(ax=ax, color='lightgrey')


# Plot each bar plot at the centroid of its respective basin
for basin_name, group in top_abund_taxa_basin.groupby('Basin', observed=True):
    # Get the centroid coordinates
    centroid = basins.loc[basins["Basin"] == basin_name].dissolve("Basin", aggfunc = "mean").geometry.centroid;
    
    # Transform the centroid coordinates to the axes coordinate system
    x_fig, y_fig = ax.transData.transform((centroid.x[0], centroid.y[0]))
    x_fig, y_fig = fig.transFigure.inverted().transform((x_fig, y_fig))
    
    # Define the position for the bar plot
    pos = [x_fig - 0.075, y_fig - 0.075, 0.15, 0.15]  # Adjust the size and position as needed
    
    # Create an inset axis for the bar plot
    inset_ax = fig.add_axes(pos)
    inset_ax.bar(group['Taxon'], group['Num_cell_l'])
    inset_ax.set_title(f'Basin: {basin}', fontsize=10)
    inset_ax.set_xticklabels(group['Taxon'], rotation=90, fontsize=8)
    inset_ax.set_aspect('auto')  # Ensure the bar plot is properly scaled

# Adjust layout
plt.tight_layout()
plt.show()