"""
"""
from __init__ import *

import subprocess as sp
import os
from scipy import sparse
import anndata


def scf_rna_format_to_h5ad(meta, gc_mat):
    """
    input:
        - meta (cell metadata)
        - gc_mat
    
    output:
        - anndata
    """
    X = gc_mat.data.T # cell by gene [scipy sparse matrix]
    obs = meta # cell annotation [pandas dataframe]
    var = pd.DataFrame(index=gc_mat.gene) # gene annotation [pandas dataframe]
    
    h5ad_mat = anndata.AnnData(X, obs, var,)
    
    return h5ad_mat

def scf_mc_format_to_h5ad(meta, mat):
    """
    input:
        - meta (cell metadata)
        - gc_mat
    
    output:
        - anndata
    """
    X = mat.T.values # cell by gene [numpy array]
    obs = meta # cell annotation [pandas dataframe]
    var = pd.DataFrame(index=mat.index) # gene annotation [pandas dataframe]
    
    h5ad_mat = anndata.AnnData(X, obs, var,)
    
    return h5ad_mat

def h5ad_to_scf_rna_format(h5ad_mat, gid_col='', cid_col=''):
    """
    input:
        - anndata
    output:
        - meta (cell metadata)
        - gc_mat

    """
    meta = h5ad_mat.obs
    if gid_col:
        genes = h5ad_mat.var[gid_col].values
    else:
        genes = h5ad_mat.var.index.values
    if cid_col:
        cells = h5ad_mat.obs[cid_col].values
    else:
        cells = h5ad_mat.obs.index.values 
    gc_mat = GC_matrix(genes,
                       cells,
                       h5ad_mat.X.T,
                      )
    return meta, gc_mat

def h5ad_to_scf_mc_format(h5ad_mat, gid_col='', cid_col=''):
    """
    input:
        - anndata
    output:
        - meta (cell metadata)
        - gc_mat
    """

    meta = h5ad_mat.obs
    if gid_col:
        genes = h5ad_mat.var[gid_col].values
    else:
        genes = h5ad_mat.var.index.values
    if cid_col:
        cells = h5ad_mat.obs[cid_col].values
    else:
        cells = h5ad_mat.obs.index.values 
    mat = pd.DataFrame(h5ad_mat.X.T,
                       index=genes,
                       columns=cells,
                      )
    return meta, mat

def diag_matrix(X, rows=np.array([]), cols=np.array([]), threshold=None):
    """Diagonalize a matrix as much as possible
    """
    di, dj = X.shape
    transposed = 0
    
    if di > dj:
        di, dj = dj, di
        X = X.T.copy()
        rows, cols = cols.copy(), rows.copy()
        transposed = 1
        
    # start (di <= dj)
    new_X = X.copy()
    new_rows = rows.copy() 
    new_cols = cols.copy() 
    if new_rows.size == 0:
        new_rows = np.arange(di)
    if new_cols.size == 0:
        new_cols = np.arange(dj)
        
    # bring the greatest values in the lower right matrix to diagnal position 
    for idx in range(min(di, dj)):

        T = new_X[idx: , idx: ]
        i, j = np.unravel_index(T.argmax(), T.shape) # get the coords of the max element of T
        
        if threshold and T[i, j] < threshold:
            dm = idx # new_X[:dm, :dm] is done (0, 1, ..., dm-1) excluding dm
            break
        else:
            dm = idx+1 # new_X[:dm, :dm] will be done

        # swap row idx, idx+i
        tmp = new_X[idx, :].copy()
        new_X[idx, :] = new_X[idx+i, :].copy() 
        new_X[idx+i, :] = tmp 
        
        tmp = new_rows[idx]
        new_rows[idx] = new_rows[idx+i]
        new_rows[idx+i] = tmp

        # swap col idx, idx+j
        tmp = new_X[:, idx].copy()
        new_X[:, idx] = new_X[:, idx+j].copy() 
        new_X[:, idx+j] = tmp 
        
        tmp = new_cols[idx]
        new_cols[idx] = new_cols[idx+j]
        new_cols[idx+j] = tmp
        
    # 
    if dm == dj:
        pass
    elif dm < dj: # free columns

        col_dict = {}
        sorted_col_idx = np.arange(dm)
        free_col_idx = np.arange(dm, dj)
        linked_rowcol_idx = new_X[:, dm:].argmax(axis=0)
        
        for col in sorted_col_idx:
            col_dict[col] = [col]
        for col, key in zip(free_col_idx, linked_rowcol_idx): 
            if key < dm:
                col_dict[key] = col_dict[key] + [col]
            else:
                col_dict[key] = [col]
                
            
        new_col_order = np.hstack([col_dict[key] for key in sorted(col_dict.keys())])
        
        # update new_X new_cols
        new_X = new_X[:, new_col_order].copy()
        new_cols = new_cols[new_col_order]
    else:
        raise ValueError("Unexpected situation: dm > dj")
    
    if transposed:
        new_X = new_X.T
        new_rows, new_cols = new_cols, new_rows
    return new_X, new_rows, new_cols 

def diag_matrix_rows(X, rows=np.array([]), cols=np.array([]),):
    """Diagonalize a matrix as much as possible by only rearrange rows
    """
    di, dj = X.shape
    
    new_X = X.copy()
    new_rows = rows.copy() 
    new_cols = cols.copy() 
    
    # free to move rows
    row_dict = {}
    free_row_idx = np.arange(di)
    linked_rowcol_idx = new_X.argmax(axis=1) # the column with max value for each row
    
    for row, key in zip(free_row_idx, linked_rowcol_idx): 
        if key in row_dict.keys():
            row_dict[key] = row_dict[key] + [row]
        else:
            row_dict[key] = [row]
            
    new_row_order = np.hstack([row_dict[key] for key in sorted(row_dict.keys())])
    # update new_X new_cols
    new_X = new_X[new_row_order, :].copy()
    new_rows = new_rows[new_row_order]
    
    return new_X, new_rows, new_cols 

def get_grad_colors(n, cmap='copper'):
    """Generate n colors from a given colormap (a matplotlib.cm)
    """
    from matplotlib import cm
    cmap = cm.get_cmap(cmap)
    return [cmap(int(i)) for i in np.linspace(0, 255, n)] 

def logcpm(counts):
    """
    Args:
        - gene-cell matrix
    """
    cov = counts.sum(axis=0)
    logcpm = np.log10(counts.divide(cov, axis=1)*1000000 + 1)
    return logcpm

def logtpm(counts, gene_lengths):
    """
    Args:
        - gene-cell matrix
        - gene_lengths: a series indexed by gene_id
    """
    tpm = counts.divide(gene_lengths.loc[counts.index], axis=0)
    cov = tpm.sum(axis=0)
    logtpm = np.log10((tpm.divide(cov, axis=1))*1000000 + 1)
    return logtpm

def sparse_logcpm(gc_matrix, mode='logcpm', lib_size=[]):
    """
    """
    lib_size = np.array(lib_size)
    if np.size(lib_size) == 0:
        lib_size = gc_matrix.data.sum(axis=0)

    lib_size_inv = sparse.diags(np.ravel(1.0/(1e-7+lib_size)))
    cpm = (gc_matrix.data).dot(lib_size_inv*1e6).tocoo()

    if mode == 'logcpm':
        cpm.data = np.log10(cpm.data + 1)
    elif mode == 'cpm':
        pass

    gc_cpm = GC_matrix(
        gc_matrix.gene, 
        gc_matrix.cell, 
        cpm,
    )
    
    return gc_cpm

def sparse_logtpm(gc_matrix, gene_lengths):
    """
    gene_lengths: array like 
    
    """
    gene_lengths = np.array(gene_lengths)
    gene_length_inv = sparse.diags(np.ravel(1.0/gene_lengths))
    tmp = (gene_length_inv).dot(gc_matrix.data).tocoo()
    lib_size_inv = sparse.diags(np.ravel(1.0/tmp.sum(axis=0)))
    
    logtpm = tmp.dot(lib_size_inv*1e6).tocoo()
    logtpm.data = np.log10(logtpm.data + 1)

    gc_logtpm = GC_matrix(
        gc_matrix.gene, 
        gc_matrix.cell, 
        logtpm,
    )
    
    return gc_logtpm

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

def create_logger(name='log'):
    """
    args: logger name

    return: a logger object
    """
    logging.basicConfig(
        format='%(asctime)s %(message)s', 
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logging.INFO)
    return logging.getLogger(name)

def set_value_by_percentile(this, lo, hi):
    """set `this` below or above percentiles to given values
    this (float)
    lo(float)
    hi(float)
    """
    if this < lo:
        return lo
    elif this > hi:
        return hi
    else:
        return this

def mcc_percentile_norm(mcc, low_p=5, hi_p=95):
    """
    set values above and below specific percentiles to be at the value of percentiles 

    args: mcc, low_p, hi_p  

    return: normalized mcc levels
    """
#   mcc_norm = [np.isnan(mcc) for mcc_i in list(mcc)]
    mcc_norm = np.copy(mcc)
    mcc_norm = mcc_norm[~np.isnan(mcc_norm)]

    lo = np.percentile(mcc_norm, low_p)
    hi = np.percentile(mcc_norm, hi_p)

    mcc_norm = [set_value_by_percentile(mcc_i, lo, hi) for mcc_i in list(mcc)]
    mcc_norm = np.array(mcc_norm)

    return mcc_norm

def plot_tsne_values(df, tx='tsne_x', ty='tsne_y', tc='mCH',
                    low_p=5, hi_p=95,
                    s=2,
                    cbar_label=None,
                    output=None, show=True, close=False, 
                    t_xlim='auto', t_ylim='auto', title=None, figsize=(8,6), **kwargs):
    """
    tSNE plot

    xlim, ylim is set to facilitate displaying glial clusters only

    """
    import matplotlib.pyplot as plt
    import seaborn as sns

    fig, ax = plt.subplots(figsize=figsize)

    im = ax.scatter(df[tx], df[ty], s=s, 
        c=mcc_percentile_norm(df[tc].values, low_p=low_p, hi_p=hi_p), **kwargs)
    if title:
        ax.set_title(title)
    else:
        ax.set_title(tc)
    ax.set_xlabel(tx)
    ax.set_ylabel(ty)
    # ax.set_aspect('auto')


    clb = plt.colorbar(im, ax=ax)
    if cbar_label:
        clb.set_label(cbar_label, rotation=270, labelpad=10)

    if t_xlim == 'auto':
        t_xlim = [np.nanpercentile(df[tx].values, 0.1), np.nanpercentile(df[tx].values, 99.9)]
        t_xlim[0] = t_xlim[0] - 0.1*(t_xlim[1] - t_xlim[0])
        t_xlim[1] = t_xlim[1] + 0.1*(t_xlim[1] - t_xlim[0])
        ax.set_xlim(t_xlim)
    elif t_xlim:
        ax.set_xlim(t_xlim)
    else:
        pass  

    if t_ylim == 'auto':
        t_ylim = [np.nanpercentile(df[ty].values, 0.1), np.nanpercentile(df[ty].values, 99.9)]
        t_ylim[0] = t_ylim[0] - 0.1*(t_ylim[1] - t_ylim[0])
        t_ylim[1] = t_ylim[1] + 0.1*(t_ylim[1] - t_ylim[0])
        ax.set_ylim(t_ylim)
    elif t_ylim:
        ax.set_ylim(t_ylim)
    else:
        pass

    fig.tight_layout()
    if output:
        fig.savefig(output)
        print('Saved to ' + output) 
    if show:
        plt.show()
    if close:
        plt.close(fig)

def get_kwcolors(labels, colors):
    """Generate a dictinary of {label: color} using unique labels and a list of availabel colors
    """
    nc = len(colors)
    nl = len(labels)
    n_repeats = int((nl + nc - 1)/nc)
    colors = list(colors)*n_repeats
    
    kw_colors = {l:c for (l,c) in zip(labels, colors)}
    return kw_colors

def rgb2hex(r,g,b):
    """From rgb (255, 255, 255) to hex
    """
    hex = "#{:02x}{:02x}{:02x}".format(int(r),int(g),int(b))
    return hex

def gen_colors(n, l=0.6, s=0.6, colors=None):
    """Generate compatible and distinct hex colors
    """
    if not colors:
        import colorsys
        hs = np.linspace(0, 1, n, endpoint=False)
        rgbs = [rgb2hex(*(256*np.array(colorsys.hls_to_rgb(h, l, s))))
                 for h in hs]
        return rgbs
    else:
        clrs = [colors[i%len(colors)] for i in range(n)] 
        return clrs 

def myScatter(ax, df, x, y, l, 
              s=20,
              sample_frac=None,
              sample_n=None,
              legend_size=None,
              legend_kws=None,
              grey_label='unlabeled',
              shuffle=True,
              random_state=None,
              legend_mode=0,
              kw_colors=False,
              colors=['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C8', 'C9'], **kwargs):
    """
    take an axis object and make a scatter plot

    - kw_colors is a dictinary {label: color}
    """

    import matplotlib.pyplot as plt
    import seaborn as sns
    df = df.copy()
    # shuffle (and copy) data
    if sample_n:
        df = (df.groupby(l).apply(lambda x: x.sample(min(len(x), sample_n), random_state=random_state))
                            .reset_index(level=0, drop=True)
            )
    if sample_frac:
        df = (df.groupby(l).apply(lambda x: x.sample(frac=sample_frac, random_state=random_state))
                            .reset_index(level=0, drop=True)
            )
    if shuffle:
        df = df.sample(frac=1, random_state=random_state)

    if not kw_colors:
        # add a color column
        inds, catgs = pd.factorize(df[l])
        df['c'] = [colors[i%len(colors)] if catgs[i]!=grey_label else 'grey' 
                    for i in inds]
    else:
        df['c'] = [kw_colors[i] if i!=grey_label else 'grey' for i in df[l]]
    
    # take care of legend
    if legend_mode != -1:
        for ind, row in df.groupby(l).first().iterrows():
            ax.scatter(row[x], row[y], c=row['c'], label=ind, s=s, **kwargs)
        
    if legend_mode == -1:
        pass
    elif legend_mode == 0:
        lgnd = ax.legend()
    elif legend_mode == 1:
        # Shrink current axis's height by 10% on the bottom
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.1,
                         box.width, box.height * 0.9])
        lgnd = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.07),
              ncol=6, fancybox=False, shadow=False) 
    elif legend_mode == 2:
        # Shrink current axis's width by 10% on the bottom
        box = ax.get_position()
        ax.set_position([box.x0 + box.width*0.1, box.y0,
                                 box.width*0.8, box.height])

    if legend_kws:
        lgnd = ax.legend(**legend_kws)

    if legend_mode != -1 and legend_size:
        for handle in lgnd.legendHandles:
            handle._sizes = [legend_size] 

    # backgroud (grey)
    df_grey = df.loc[df['c']=='grey']
    if not df_grey.empty:
        ax.scatter(df_grey[x], 
                   df_grey[y],
                   c=df_grey['c'], s=s, **kwargs)
    # actual plot
    df_tmp = df.loc[df['c']!='grey']
    ax.scatter(df_tmp[x], 
               df_tmp[y],
               c=df_tmp['c'], s=s, **kwargs)
    
    return

def plot_tsne_labels_ax(df, ax, tx='tsne_x', ty='tsne_y', tc='cluster_ID', 
                    sample_frac=None,
                    sample_n=None,
                    legend_size=None,
                    legend_kws=None,
                    grey_label='unlabeled',
                    legend_mode=0,
                    s=1,
                    shuffle=True,
                    random_state=None,
                    t_xlim='auto', t_ylim='auto', title=None, 
                    legend_loc='lower right',
                    colors=['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C8', 'C9'], **kwargs):
    """
    tSNE plot

    xlim, ylim is set to facilitate displaying glial clusters only

    # avoid gray-like 'C7' in colors
    # color orders are arranged for exci-inhi-glia plot 11/1/2017
    """
    import matplotlib.pyplot as plt

    myScatter(ax, df, tx, ty, tc,
             s=s,
             sample_frac=sample_frac,
             sample_n=sample_n,
             legend_size=legend_size,
             legend_kws=legend_kws,
             shuffle=shuffle,
             grey_label=grey_label,
             random_state=random_state, 
             legend_mode=legend_mode, 
             colors=colors, **kwargs)

    if title:
        ax.set_title(title)
    else:
        ax.set_title(tc)
    ax.set_xlabel(tx)
    ax.set_ylabel(ty)
    # ax.set_aspect('auto')

    if t_xlim == 'auto':
        t_xlim = [np.nanpercentile(df[tx].values, 0.1), np.nanpercentile(df[tx].values, 99.9)]
        t_xlim[0] = t_xlim[0] - 0.1*(t_xlim[1] - t_xlim[0])
        t_xlim[1] = t_xlim[1] + 0.1*(t_xlim[1] - t_xlim[0])
        ax.set_xlim(t_xlim)
    elif t_xlim:
        ax.set_xlim(t_xlim)
    else:
        pass  

    if t_ylim == 'auto':
        t_ylim = [np.nanpercentile(df[ty].values, 0.1), np.nanpercentile(df[ty].values, 99.9)]
        t_ylim[0] = t_ylim[0] - 0.1*(t_ylim[1] - t_ylim[0])
        t_ylim[1] = t_ylim[1] + 0.1*(t_ylim[1] - t_ylim[0])
        ax.set_ylim(t_ylim)
    elif t_ylim:
        ax.set_ylim(t_ylim)
    else:
        pass

    return


def plot_tsne_labels(df, tx='tsne_x', ty='tsne_y', tc='cluster_ID', 
                    grey_label='unlabeled',
                    sample_frac=None,
                    sample_n=None,
                    legend_size=None,
                    legend_mode=0,
                    legend_kws=None,
                    s=1,
                    random_state=None,
                    output=None, show=True, close=False, 
                    t_xlim='auto', t_ylim='auto', title=None, figsize=(8,6),
                    legend_loc='lower right',
                    colors=['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C8', 'C9'], **kwargs):
    """
    tSNE plot

    xlim, ylim is set to facilitate displaying glial clusters only

    # avoid gray-like 'C7' in colors
    # color orders are arranged for exci-inhi-glia plot 11/1/2017
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    fig, ax = plt.subplots(figsize=figsize)

    myScatter(ax, df, tx, ty, tc,
             s=s,
             sample_frac=sample_frac,
             sample_n=sample_n,
             legend_size=legend_size,
             legend_kws=legend_kws,
             grey_label=grey_label,
             random_state=random_state, 
             legend_mode=legend_mode, 
             colors=colors, **kwargs)

    if title:
        ax.set_title(title)
    else:
        ax.set_title(tc)
    ax.set_xlabel(tx)
    ax.set_ylabel(ty)
    # ax.set_aspect('auto')

    if t_xlim == 'auto':
        t_xlim = [np.nanpercentile(df[tx].values, 0.1), np.nanpercentile(df[tx].values, 99.9)]
        t_xlim[0] = t_xlim[0] - 0.1*(t_xlim[1] - t_xlim[0])
        t_xlim[1] = t_xlim[1] + 0.1*(t_xlim[1] - t_xlim[0])
        ax.set_xlim(t_xlim)
    elif t_xlim:
        ax.set_xlim(t_xlim)
    else:
        pass  

    if t_ylim == 'auto':
        t_ylim = [np.nanpercentile(df[ty].values, 0.1), np.nanpercentile(df[ty].values, 99.9)]
        t_ylim[0] = t_ylim[0] - 0.1*(t_ylim[1] - t_ylim[0])
        t_ylim[1] = t_ylim[1] + 0.1*(t_ylim[1] - t_ylim[0])
        ax.set_ylim(t_ylim)
    elif t_ylim:
        ax.set_ylim(t_ylim)
    else:
        pass

    if output:
        fig.savefig(output)
        print('Saved to ' + output) 
    if show:
        plt.show()
    if close:
        plt.close(fig)

def plot_tsne_values_ax(df, ax, tx='tsne_x', ty='tsne_y', tc='mCH',
                    low_p=5, hi_p=95,
                    s=2,
                    cbar=True,
                    cbar_ax=None,
                    cbar_label=None,
                    t_xlim='auto', t_ylim='auto', title=None, **kwargs):
    """
    tSNE plot

    xlim, ylim is set to facilitate displaying glial clusters only

    """
    import matplotlib.pyplot as plt


    im = ax.scatter(df[tx], df[ty], s=s, 
        c=mcc_percentile_norm(df[tc].values, low_p=low_p, hi_p=hi_p), **kwargs)
    if title:
        ax.set_title(title)
    else:
        ax.set_title(tc)
    # ax.set_aspect('auto')
    if cbar:
        if cbar_ax:
            clb = plt.colorbar(im, cax=cbar_ax, shrink=0.4)
        else:
            clb = plt.colorbar(im, cax=ax, shrink=1)
        if cbar_label:
            clb.set_label(cbar_label, rotation=270, labelpad=10)

    if t_xlim == 'auto':
        t_xlim = [np.nanpercentile(df[tx].values, 0.1), np.nanpercentile(df[tx].values, 99.9)]
        t_xlim[0] = t_xlim[0] - 0.1*(t_xlim[1] - t_xlim[0])
        t_xlim[1] = t_xlim[1] + 0.1*(t_xlim[1] - t_xlim[0])
        ax.set_xlim(t_xlim)
    elif t_xlim:
        ax.set_xlim(t_xlim)
    else:
        pass  

    if t_ylim == 'auto':
        t_ylim = [np.nanpercentile(df[ty].values, 0.1), np.nanpercentile(df[ty].values, 99.9)]
        t_ylim[0] = t_ylim[0] - 0.1*(t_ylim[1] - t_ylim[0])
        t_ylim[1] = t_ylim[1] + 0.1*(t_ylim[1] - t_ylim[0])
        ax.set_ylim(t_ylim)
    elif t_ylim:
        ax.set_ylim(t_ylim)
    else:
        pass

    return im  


def get_mcc(df, base_call_cutoff=100, sufficient_coverage_fraction=1, suffix=True, fillna=True):
    """Get mcc matrix from mc_c matrix (filtering out low coverage gene or bins)
    """
    logging.info('Getting mcc matrix from mc and c') 
    logging.info('base_call_cutoff={}, sufficient_coverage_fraction={}'.format(
                base_call_cutoff, sufficient_coverage_fraction))
    
    df_c = df.filter(regex="_c$")
    df_c.columns = [col[:-len('_c')] for col in df_c.columns] 
    df_mc = df.filter(regex="_mc$")
    df_mc.columns = [col[:-len('_mc')] for col in df_mc.columns] 
    # a gene is sufficiently covered in % of cells 
    condition = (df_c > base_call_cutoff).sum(axis=1) >= sufficient_coverage_fraction*(df.shape[1])/2.0 

    logging.info("Matrix size before pruning... "+ str(df.shape))
    logging.info("Matrix size after pruning... "+ str(df.loc[condition].shape))
    
    # get mcc matrix with kept bins and nan values for low coverage sites
    df_c_nan = df_c.copy()
    df_c_nan[df_c < base_call_cutoff] = np.nan
    df_mcc = df_mc.loc[condition]/df_c_nan.loc[condition]
    logging.info(df_mcc.shape)

    # imputation (missing value -> mean value of all cells)
    if fillna:
        logging.info('Imputing data... (No effect if sufficient_coverage_fraction=1)')
        means = df_mcc.mean(axis=1)
        fill_value = pd.DataFrame({col: means for col in df_mcc.columns})
        df_mcc.fillna(fill_value, inplace=True)
    
    # add suffix
    if suffix:
        df_mcc.columns = df_mcc.columns.values + '_mcc'
    
    return df_mcc

def get_mcc_lite(mc_table, c_table, base_call_cutoff=100, sufficient_coverage_fraction=1, fillna=True):
    """Given 2 numpy array, return mcc table
    Gene/region by sample matrix
    """
    df_c = pd.DataFrame(c_table)
    df_mc = pd.DataFrame(mc_table)
    assert df_c.shape == df_mc.shape
    
    # a gene is sufficiently covered in % of cells 
    condition = (df_c > base_call_cutoff).sum(axis=1) >= sufficient_coverage_fraction*(df_c.shape[1])

    logging.info("Matrix size before pruning... "+ str(df_c.shape))
    logging.info("Matrix size after pruning... "+ str(df_c.loc[condition].shape))
    
    # get mcc matrix with kept bins and nan values for low coverage sites
    df_c_nan = df_c.copy()
    df_c_nan[df_c < base_call_cutoff] = np.nan
    df_mcc = df_mc.loc[condition]/df_c_nan.loc[condition]
    logging.info(df_mcc.shape)

    # imputation (missing value -> mean value of all cells)
    if fillna:
        logging.info('Imputing data... (No effect if sufficient_coverage_fraction=1)')
        means = df_mcc.mean(axis=1)
        fill_value = pd.DataFrame({col: means for col in df_mcc.columns})
        df_mcc.fillna(fill_value, inplace=True)
    
    # return matrix and index (regions)
    return df_mcc.values, df_mcc.index.values

def get_mcc_lite_v2(df_c, df_mc, base_call_cutoff):
    """
    """
    # get mcc matrix with kept bins and nan values for low coverage sites
    df_c_nan = df_c.copy()
    df_c_nan[df_c < base_call_cutoff] = np.nan
    df_mcc = df_mc/df_c_nan
    logging.info(df_mcc.shape)

    # imputation (missing value -> mean value of all cells)
    means = df_mcc.mean(axis=1)
    fill_value = pd.DataFrame({col: means for col in df_mcc.columns})
    df_mcc.fillna(fill_value, inplace=True)
    
    return df_mcc

def get_mcc_lite_v3(df_c, df_mc, base_call_cutoff):
    """
    """
    # get mcc matrix with kept bins and nan values for low coverage sites
    df_c_nan = df_c.copy()
    df_c_nan[df_c < base_call_cutoff] = np.nan
    df_mcc = df_mc/df_c_nan
    return df_mcc


def get_clusters_mc_c_worker(df_cells, df_input, cluster_col):
    """reduce gene*cell or bin*cell matrix to a gene*cluster or bin*cluster matrix
    Arguments:
        - df_cells: a dataframe indexed by 'cell_name', and have '$cluster_col' as column
        - df_input: a dataframe with 'sample_mc', 'sample_c' ... as columns
        sample names are cell names
    """
    # cluster mc_c
    df_c = df_input.filter(regex='_c$')
    df_mc = df_input.filter(regex='_mc$')

    df_mc_c = pd.DataFrame() 
    for label, df_sub in df_cells.groupby(cluster_col):
        samples = df_sub.index.values
        df_mc_c['{}_mc'.format(label)] = df_mc[samples+'_mc'].sum(axis=1)
        df_mc_c['{}_c'.format(label)] = df_c[samples+'_c'].sum(axis=1)

    logging.info("Output shape: {}".format(df_mc_c.shape))
    return df_mc_c

def rank_array(array):
    """Return ranking of each element of an array
    """
    array = np.array(array)
    temp = array.argsort()
    ranks = np.empty_like(temp)
    ranks[temp] = np.arange(len(array))
    return ranks

# added 4/5/2019
def rank_rows(matrix):
    """Return rankings of each rwo in a 2d array
    """
    matrix = np.array(matrix)
    return np.apply_along_axis(rank_array, 1, matrix) # row = 1

def spearman_corrcoef(X, Y):
    """return spearman correlation matrix for each pair of rows of X and Y
    """
    return np.corrcoef(rank_rows(X), rank_rows(Y))

def spearmanr_paired_rows(X, Y):
    from scipy import stats
    
    X = np.array(X)
    Y = np.array(Y)
    corrs = []
    ps = []
    for x, y in zip(X, Y):
        r, p = stats.spearmanr(x, y)
        corrs.append(r)
    return np.array(corrs), np.array(ps)

def get_index_from_array(arr, inqs, na_rep=-1):
    """Get index of array
    """
    arr = np.array(arr)
    arr = pd.Series(arr).reset_index().set_index(0)
    idxs = arr.reindex(inqs)['index'].fillna(na_rep).astype(int).values
    return idxs

def get_genomic_distance(sa, ea, sb, eb):
    """Get genomic distance
    """
    assert sa < ea and sb < eb
    if sa > sb:
        sa, sb = sb, sa
        ea, eb = eb, ea
        
    # sa <= sb
    distance = max(0, sb - ea)
    
    return distance

def get_reverse_comp(string):
    """Get reverse compliment of a string
    """
    comp_dict = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G',
        'N': 'N',
    }
    for char in set(string):
        if char not in ['A', 'C', 'G', 'T', 'N']:
            raise ValueError('Not allowed char in string')
            
    new_string = ''.join([comp_dict[char] for char in string[::-1]])
    return new_string
    
# added 4/11/2019
def save_gxc_matrix(gxc, f_mat, f_gene, f_cell):
    """
    """
    print("Deprecated: 6/29/2019 Use save_gc_matrix instead")
    sparse.save_npz(f_mat, gxc.data)
    with open(f_gene, 'w') as f:
        f.write('\n'.join(gxc.gene)+'\n')
    with open(f_cell, 'w') as f:
        f.write('\n'.join(gxc.cell)+'\n')
        
def save_gxc_matrix_methylation(gxc, f_mat_c, f_mat_mc, f_gene, f_cell):
    """
    """
    print("Deprecated: 6/29/2019 Use save_gc_matrix_methylation instead")
    sparse.save_npz(f_mat_c, gxc.data['c'])
    sparse.save_npz(f_mat_mc, gxc.data['mc'])
    with open(f_gene, 'w') as f:
        f.write('\n'.join(gxc.gene)+'\n')
    with open(f_cell, 'w') as f:
        f.write('\n'.join(gxc.cell)+'\n') 

def save_gc_matrix(gc_matrix, f_gene, f_cell, f_mat):
    """
    """
    sparse.save_npz(f_mat, gc_matrix.data)
    with open(f_gene, 'w') as f:
        f.write('\n'.join(gc_matrix.gene)+'\n')
    with open(f_cell, 'w') as f:
        f.write('\n'.join(gc_matrix.cell)+'\n')

def save_gc_matrix_methylation(gc_matrix, f_gene, f_cell, f_mat_mc, f_mat_c):
    """
    """
    sparse.save_npz(f_mat_mc, gc_matrix.data['mc'])
    sparse.save_npz(f_mat_c, gc_matrix.data['c'])
    with open(f_gene, 'w') as f:
        f.write('\n'.join(gc_matrix.gene)+'\n')
    with open(f_cell, 'w') as f:
        f.write('\n'.join(gc_matrix.cell)+'\n') 

def import_single_textcol(fname, header=None, col=0):
    return pd.read_csv(fname, header=header, sep='\t')[col].values

def export_single_textcol(fname, array):
    with open(fname, 'w') as f:
        f.write('\n'.join(array)+'\n')

def load_gc_matrix(f_gene, f_cell, f_mat):
    """
    """
    gene = import_single_textcol(f_gene)
    cell = import_single_textcol(f_cell)
    mat = sparse.load_npz(f_mat) 
    assert (len(gene), len(cell)) == mat.shape
    return GC_matrix(gene, cell, mat) 

def load_gc_matrix_methylation(f_gene, f_cell, f_mat_mc, f_mat_c):
    """
    """
    _gene = import_single_textcol(f_gene) 
    _cell = import_single_textcol(f_cell)
    _mat_mc = sparse.load_npz(f_mat_mc) 
    _mat_c = sparse.load_npz(f_mat_c) 
    gxc_raw = GC_matrix(_gene, _cell, 
                              {'c': _mat_c, 'mc': _mat_mc})
    return gxc_raw

def nondup_legends(ax='', **kwargs):
    """Assuming plt (matplotlib.pyplot) is imported
    """
    from collections import OrderedDict
    import matplotlib.pyplot as plt

    if ax == '':
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = OrderedDict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys(), **kwargs)
    else:
        handles, labels = ax.get_legend_handles_labels()
        by_label = OrderedDict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys(), **kwargs)
    return 

def dedup_array_elements(x, empty_string=''):
    """Replacing repeats with empty_string
    """
    newx = np.empty_like(x)
    newx[0] = x[0]
    for i in range(1, len(x)):
        if x[i-1] == x[i]:
            newx[i] = empty_string
        else:
            newx[i] = x[i]
    return newx

def vcorrcoef(X,Y):
    """Compute correlation coef for each rows of X and Y
    """
    assert X.shape == Y.shape
    Xm = np.mean(X,axis=1).reshape(-1,1)
    Ym = np.mean(Y,axis=1).reshape(-1,1)
    Xm = X-Xm
    Ym = Y-Ym
    
    r_num = np.sum(Xm*Ym,axis=1)
    r_den = np.sqrt(np.sum(Xm**2,axis=1)*np.sum(Ym**2, axis=1))
    r = r_num/r_den
    return r

def zscore(x, offset=1e-7, ddof=1):
    return (x - np.mean(x))/(np.std(x, ddof=ddof) + offset)


def clst_umap_pipe_lite(pcs, cells_all, 
                        resolution=1,
                        npc=50,
                        k=30,
                        verbose=False, seed=0, cluster_only=False, 
                       ):
    # clustering
    import CEMBA_clst_utils
    import CEMBA_run_tsne

    df_clst = CEMBA_clst_utils.clustering_routine(
                                    pcs, 
                                    cells_all, k, 
                                    verbose=verbose,
                                    resolution=resolution,
                                    seed=seed,
                                    metric='euclidean', option='plain', n_trees=10, search_k=-1)

    # umap
    if not cluster_only:
        df_tsne = CEMBA_run_tsne.run_umap_lite(
                    pcs, 
                    cells_all, 
                    verbose=verbose,
                    n_neighbors=30, min_dist=0.5, n_dim=2, 
                    random_state=1)

        df_summary = df_clst.join(df_tsne)
        return df_summary
    else:
        return df_clst

def gen_cdf(array, ax, x_range=[], n_points=1000, show=True, flip=False, **kwargs):
    """
    """
    x = np.sort(array)
    y = np.arange(len(array))/len(array)
    if flip:
        # x = x[::-1]
        y = 1 - y

    if not x_range:
        if show:
            ax.plot(x, y, **kwargs)
        return x, y 
    else:
        start, end = x_range
        xbins = np.linspace(start, end, n_points)
        ybins = np.interp(xbins, x, y)
        if show:
            ax.plot(xbins, ybins, **kwargs)
        return xbins, ybins 

def savefig(fig, path):
    """
    """
    fig.savefig(path, bbox_inches='tight', dpi=300)
    return 
