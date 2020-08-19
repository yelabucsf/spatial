from .imports import *
from pdb import set_trace


def subselect_data(data, celltype, grouping_fn):
    '''
    Returns a subselection of the data with a particular celltype and grouping
    function
    '''
    data_cell = data[data.obs.loc[data.obs.celltype == celltype].index].copy()
    grouping = grouping_fn(data_cell)
    data_cell.obs["grouping"] = grouping
    return data_cell


def plot_grouping(cell_spec_data, data, layer, perc):
    '''
    Plots the hexagonal grid colored by grouping labels
    '''

    perc_label = int(100*perc)
    blue_label = f"Percent Tumor Within {layer} layers > {perc_label}%"
    red_label = f"Percent Tumor Within {layer} layers < {perc_label}%"



    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(20,6))
    ax1 = axes[0]
    ax2 = axes[1]

    sub_df = cell_spec_data.obs.loc[~cell_spec_data.obs.grouping]
    x = sub_df.X 
    y = sub_df.Y - (sub_df.X % 2) / 2
    red_sc = ax1.scatter(x,y, c="Red", cmap='bwr', marker='H',s=20*10, label=red_label)
    sub_df = cell_spec_data.obs.loc[cell_spec_data.obs.grouping]
    x = sub_df.X 
    y = sub_df.Y - (sub_df.X % 2) / 2
    blue_sc = ax1.scatter(x,y, c="Blue", cmap='bwr', marker='H',s=20*10, label=blue_label)
    ax1.set(xlim=(0,42), ylim=(0,20))
    ax1.invert_yaxis()
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles[::-1], labels[::-1], fontsize=15)
    #ax1.legend(, fontsize=15)



    sub_df = data.obs.loc[data.obs.celltype=='Hepatocyte']
    x = sub_df.X 
    y = sub_df.Y - (sub_df.X % 2) / 2
    ax2.scatter(x,y, c="Red", marker='H',s=20*10, label="Hepatocyte")

    sub_df = data.obs.loc[data.obs.celltype=='mc38']
    x = sub_df.X 
    y = sub_df.Y - (sub_df.X % 2) / 2
    ax2.scatter(x,y, c="Blue", marker='H',s=20*10, label='Cancer')

    ax2.set(xlim=(0,42), ylim=(0,20))
    ax2.invert_yaxis()
    handles, labels = ax2.get_legend_handles_labels()
    ax2.legend(handles[::-1], labels[::-1], fontsize=15)
    return fig


def perform_mann_whitney(hepatocytes, all_genes):
    '''
    Runs and plots results of mann_whitney test
    '''
    mann_whitney = de.test.rank_test(hepatocytes, grouping="grouping")
    idxs = np.argsort(mann_whitney.pval)
    pvals_sorted = mann_whitney.pval[idxs]
    qvals_sorted = mann_whitney.qval[idxs]
    genes_sorted = mann_whitney._gene_names[idxs]
    fc_sorted = mann_whitney.log2_fold_change()[idxs]
    genes_to_write = [x[5:].upper() for x in genes_sorted]
    genes_to_write_set = set(genes_to_write)
    other_genes = [x[5:].upper() for x in all_genes if x[5:].upper() not in genes_to_write_set]
    df = pd.DataFrame(list(zip(genes_sorted, pvals_sorted, qvals_sorted)), columns=["Gene", "Pval", "Qval"])

    #mann_whitney.plot_volcano(corrected_pval=True, size=40, save=save)

    # idxs = np.argsort(mann_whitney.pval)
    # pvals_sorted = mann_whitney.pval[idxs]
    # genes_sorted = mann_whitney._gene_names[idxs]
    # fc_sorted = mann_whitney.log2_fold_change()[idxs]


    # mean_blue_label = np.mean(np.asarray(hepatocytes[hepatocytes.obs.grouping].X.todense()[:, idxs[0]]).reshape(-1))
    # mean_red_label = np.mean(np.asarray(hepatocytes[~hepatocytes.obs.grouping].X.todense()[:, idxs[0]]).reshape(-1))


    # median_blue_label = np.median(np.asarray(hepatocytes[hepatocytes.obs.grouping].X.todense()[:, idxs[0]]).reshape(-1))
    # median_red_label = np.median(np.asarray(hepatocytes[~hepatocytes.obs.grouping].X.todense()[:, idxs[0]]).reshape(-1))

    # y = np.asarray(hepatocytes.X.todense()[:, idxs[0]]).reshape(-1)
    # x = hepatocytes.obs.grouping
    #ax = sns.boxplot(x,y)

    # print(f"1.  {genes_sorted[0][5:]}, Fold Change: {fc_sorted[0]:.2f}")
    # print(f"2.  {genes_sorted[1][5:]}, Fold Change: {fc_sorted[1]:.2f}")
    # print(f"3.  {genes_sorted[2][5:]}, Fold Change: {fc_sorted[2]:.2f}")
    # print(f"4.  {genes_sorted[3][5:]}, Fold Change: {fc_sorted[3]:.2f}")
    # print(f"5.  {genes_sorted[4][5:]}, Fold Change: {fc_sorted[4]:.2f}")
    # print(f"5.  {genes_sorted[5][5:]}, Fold Change: {fc_sorted[5]:.2f}")
    # print(f"6.  {genes_sorted[6][5:]}, Fold Change: {fc_sorted[6]:.2f}")
    # print(f"7.  {genes_sorted[7][5:]}, Fold Change: {fc_sorted[7]:.2f}")
    # print(f"8.  {genes_sorted[8][5:]}, Fold Change: {fc_sorted[8]:.2f}")
    # print(f"9.  {genes_sorted[9][5:]}, Fold Change: {fc_sorted[9]:.2f}")
    # print(f"10. {genes_sorted[10][5:]}, Fold Change: {fc_sorted[10]:.2f} ")
    # print(f"11. {genes_sorted[11][5:]}, Fold Change: {fc_sorted[11]:.2f}")

    # print()
    # print(f"Mean Expression {blue_label}: {mean_blue_label:.2f}")
    # print(f"Mean Expression {red_label}: {mean_red_label:.2f}")
    # print()
    # print(f"Median Expression {blue_label}: {median_blue_label:.2f}")
    # print(f"Median Expression {red_label}: {median_red_label:.2f}")

    return mann_whitney, df



def get_nbr_coords(center):
    '''
    Get coordinates of the first layer around a cell
    '''
    x, y = center
    nbrs = [(x, y+1), (x, y-1), (x+1, y+0.5), (x+1, y-0.5), (x-1, y+0.5), (x-1, y-0.5)]
    return set(nbrs)

def get_layer_coords(center, max_layer):
    '''
    Get all points within a particular layer
    '''
    if max_layer == 0:
        return [center]
    total = get_nbr_coords(center) 
    layer_points = total.copy()
    if max_layer==1:
        return layer_points
    for layer in range(2, max_layer+1):
        next_layer = set([nbr for point in layer_points for nbr in get_nbr_coords(point)])
        total.update(x for x in next_layer)
        layer_points = next_layer
    return total



def get_average_stats(data, nbrs, celltypes):
    num_cells = 0
    num_celltypes = 0
    for nbr in nbrs:
        sub_df = data.obs.loc[data.obs.Grid_Point == nbr]
        num_cells += sub_df.shape[0]
        celltype_series = [sub_df.celltype == celltype for celltype in celltypes]
        num_celltypes += pd.DataFrame(celltype_series).sum().sum()
    if num_cells == 0:
        return 0.0
    return num_celltypes / num_cells


def get_in_well_stats(data, nbrs, celltypes):
    num_wells = 0
    num_has_celltypes = 0
    for nbr in nbrs:
        sub_df = data.obs.loc[data.obs.Grid_Point == nbr]
        num_wells += 1
        celltype_series = [sub_df.celltype == celltype for celltype in celltypes]
        num_celltypes = pd.DataFrame(celltype_series).sum().sum()
        num_has_celltypes += 1 if num_celltypes > 0 else 0
    return num_has_celltypes / num_wells


def add_tumor_spatial_data(data):
    data.obs["Grid_Point"] = data.obs.loc[:, ["X_hex", "Y_hex"]].apply(lambda row: (row[0], row[1]), axis=1)
    max_layers = 2
    batch_dfs = []
    for batch in tqdm(range(5), total=5, desc="Processing Tumor Spatial Information"):
        sub_data = data[data.obs[f'batch_{batch}'] == 1].copy()
        sub_df = data.obs.loc[data.obs[f'batch_{batch}'] == 1].copy()

        grid_points = sub_df.Grid_Point.unique()
        grid_to_layer_to_percent_tumor = defaultdict(dict)
        for layer in range(max_layers):
            for grid_point in grid_points:
                nbrs = get_layer_coords(grid_point, layer)
                percent_tumor = get_average_stats(sub_data, nbrs, celltypes=['mc38'])
                grid_to_layer_to_percent_tumor[grid_point][layer] = percent_tumor

        for layer in range(max_layers):
            sub_df[f"Percent_Tumor_{layer}"] = sub_df.Grid_Point.map(lambda x: grid_to_layer_to_percent_tumor[x][layer])

        grid_to_layer_to_wells_with_tumor = defaultdict(dict)
        for layer in range(max_layers):
            for grid_point in grid_points:
                nbrs = get_layer_coords(grid_point, layer)
                percent_tumor = get_in_well_stats(sub_data, nbrs, celltypes=["mc38"])
                grid_to_layer_to_wells_with_tumor[grid_point][layer] = percent_tumor

        for layer in range(max_layers):
            sub_df[f"Percent_Wells_With_Tumor_{layer}"] = sub_df.Grid_Point.map(lambda x: grid_to_layer_to_wells_with_tumor[x][layer])
        batch_dfs.append(sub_df)

    df = pd.concat(batch_dfs)
    new_columns = [c for c in df.columns if c not in data.obs.columns]
    data.obs = data.obs.merge(df.loc[:,new_columns], how='outer', left_index=True, right_index=True)
    return data
