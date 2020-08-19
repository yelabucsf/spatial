from .imports import *

def get_grid_to_celltype(data, celltype):
    '''
    Return a dictionary that maps the gridpoints to whether or not 
    they have a celltype in that gridpoint.
    '''
    grid_to_celltype = defaultdict(bool)
    for k, v in zip(zip(data.obs.X, data.obs.Y), (data.obs.celltype == celltype)):
        if v:
            grid_to_celltype[k] = True 
    return grid_to_celltype


def add_celltype_columns(data):
    '''
    Adds columns to data that says whether or not a particular celltype exists
    at a gridpoint
    '''
    immune_cells = set(['Marcrophage', 'Lymphocyte'])
    data.obs["Is_Immune"] = data.obs.celltype.map(lambda x: x in immune_cells)
    grid_to_has_immune = (data.obs.groupby(["X", "Y"]).Is_Immune.sum() > 0).to_dict()
    data.obs["Has_Immune"] = [grid_to_has_immune[k] for k in zip(data.obs.X, data.obs.Y)]

    grid_to_has_macrophage = get_grid_to_celltype(data, 'Macrophage')
    data.obs["Has_Macrophage"] = [grid_to_has_macrophage[k] for k in zip(data.obs.X, data.obs.Y)]

    grid_to_has_lymphocyte = get_grid_to_celltype(data, 'Lymphocyte')
    data.obs["Has_Lymphocyte"] = [grid_to_has_lymphocyte[k] for k in zip(data.obs.X, data.obs.Y)]

    grid_to_has_kupffer = get_grid_to_celltype(data, 'Kupffer')
    data.obs["Has_Kupffer"] = [grid_to_has_kupffer[k] for k in zip(data.obs.X, data.obs.Y)]
    return data

def add_hex_coords(data):
    '''
    Adds the hexagonal grid coordinates from the (X,Y) coordinates provided
    '''
    x = data.obs.X 
    y = data.obs.Y - (data.obs.X % 2) / 2
    data.obs["X_hex"] = x
    data.obs["Y_hex"] = y
    return data 


def load_data_for_notebook():
    data_path = Path('../data')
    data = sc.read(data_path / 'concat_20200311.loom')
    data.obs.rename({"CellType": "celltype"}, axis=1, inplace=True)
    
    # keep only mouse genes
    genes_to_keep = data.var.index[data.var.index.map(lambda x: x.startswith('mm'))]
    data = data[:, genes_to_keep].copy()
    # Remove log transformation 
    data.X = np.exp(data.X.todense()) - 1
    data.X = scipy.sparse.csr_matrix(data.X)

    data = add_grid_data(data)
    data = add_celltype_columns(data)

    return data

def add_grid_data(data):
    data_path = Path('../data')
    df = pd.read_csv(data_path / 'plate23_map.csv', header=None)
    df.columns=["barcode", "X", "Y"]
    barcode_to_xy = {}
    for b, x, y in df.to_numpy():
        assert(b not in barcode_to_xy)
        barcode_to_xy[b] = (x,y)
    data.obs['X']= data.obs.barcode.map(lambda b: barcode_to_xy[b][0])
    data.obs['Y']= data.obs.barcode.map(lambda b: barcode_to_xy[b][1])
    data = add_hex_coords(data)
    return data


def filter_genes_mincount_hepatocytes(data, count_min=100):
    '''
    Filter gene if it does not have `count_min` counts among hepatocytes
    '''
    row_inds = data.obs[data.obs.celltype == "Hepatocyte"].index
    gene_has_counts = np.asarray(data[row_inds,:].X.todense().sum(axis=0) > count_min).reshape(-1)
    genes_to_keep = data.var.index[gene_has_counts]
    data = data[:, genes_to_keep].copy()
    return data


def normalize_data(data):
    median_reads = data.obs.n_counts.median()
    sc.pp.normalize_total(data, target_sum=median_reads)
    return data
