"""Utilities for SingleCellFusion
"""
from __init__ import *

import functools
import collections
import itertools
import re
from scipy import sparse
from scipy.stats import zscore
import fbpca
import sys
import logging

import basic_utils
import clst_utils

def sparse_adj_to_mat(adjs, row_size, col_size, dists=''):
    """Turn a knn adjacency matrix to a sparse matrix
    """
    n_obs, k = adjs.shape
    assert n_obs == row_size
    # row col 1/dist 
    row_inds = np.repeat(np.arange(row_size), k)
    col_inds = np.ravel(adjs)
    if isinstance(dists, np.ndarray):
        assert dists.shape == adjs.shape
        data = np.ravel(dists) 
    else:
        data = [1]*len(row_inds)
    knn_dist_mat = sparse.coo_matrix((data, (row_inds, col_inds)), shape=(row_size, col_size))
    return knn_dist_mat

# smooth-within modality
def smooth_in_modality(counts_matrix, norm_counts_matrix, k, ka, npc=100, sigma=1.0, p=0.1, drop_npc=0, 
    output_filename="",
    ):
    """Smooth a data matrix
    
    Arguments:
        - counts_matrix (pandas dataframe, feature by cell)
        - norm_counts_matrix (pandas dataframe, feature by cell) log10(CPM+1)
        - k (number of nearest neighbors)
        - ka normalized by the ka'th neighbor (including self)
    Return:
        - smoothed cells_matrix (pandas dataframe)
        - markov affinity matrix
    """
    # from sklearn.neighbors import NearestNeighbors
    import fbpca
    import clst_utils
    
    assert counts_matrix.shape[1] == norm_counts_matrix.shape[1] 

    c = norm_counts_matrix.columns.values
    N = len(c)

    # reduce dimension fast version
    U, s, Vt = fbpca.pca(norm_counts_matrix.T.values, k=npc)
    pcs = U.dot(np.diag(s))
    if drop_npc != 0:
        pcs = pcs[:, drop_npc:]

    # get k nearest neighbor distances fast version 
    inds, dists = clst_utils.gen_knn_annoy(pcs, k, form='list', 
                                                metric='euclidean', n_trees=10, search_k=-1, verbose=True, 
                                                include_distances=True)
    
    # remove itself
    dists = dists[:, 1:]
    inds = inds[:, 1:]

    # normalize by ka's distance 
    # (including itself; needs to -2; one to switch to 0-based index, one account for the remove itself step) 
    dists = (dists/(dists[:, ka-2].reshape(-1, 1)))

    # gaussian kernel
    adjs = np.exp(-((dists**2)/(sigma**2))) 

    # construct a sparse matrix 
    cols = np.ravel(inds)
    rows = np.repeat(np.arange(N), k-1) # remove itself
    vals = np.ravel(adjs)
    A = sparse.csr_matrix((vals, (rows, cols)), shape=(N, N))

    # Symmetrize A (union of connection)
    A = A + A.T

    # normalization fast (A is now a weight matrix excluding itself)
    degrees = A.sum(axis=1)
    A = sparse.diags(1.0/np.ravel(degrees)).dot(A)

    # include itself
    eye = sparse.identity(N)
    A = p*eye + (1-p)*A
    
    # smooth fast (future?)
    counts_matrix_smoothed = pd.DataFrame((A.dot(counts_matrix.T)).T, 
                                         columns=counts_matrix.columns, index=counts_matrix.index)
    if output_filename:
        # save A (name, axis?)
        sparse.save_npz(output_filename, A)
    return counts_matrix_smoothed, A

# impute across modality
def get_constrained_knn(mat_norm_j, mat_norm_i, knn, k_saturate, knn_speed_factor=10, metric='dot', verbose=False, 
    ):
    """Get constrained knn
    j <- i
    Look for kNN in i for each cell in j, cells in i are constrained to k_saturated
    
    get knn_speed_factor*knn number of nearest neighbors internally
    """
    ti = time.time()
    assert mat_norm_i.shape[1] == mat_norm_j.shape[1]
    knn = min(int(knn), len(mat_norm_i)) # min added 08/11/2020
    knn_speed_factor = int(knn_speed_factor)
    
    cells_i = np.arange(len(mat_norm_i))
    cells_j = np.arange(len(mat_norm_j))
    
    # record cells in j
    accepted_knn_ji = [] 
    accepted_cells = []
    rejected_cells = np.arange(len(cells_j))
    
    # record cell in i
    n_connects = np.zeros(len(cells_i)).astype(int) # record number of connection for each cell in i 
    unsaturated = (n_connects < k_saturate) # unsaturated bool 
    unsaturated_cells = np.arange(len(cells_i))[unsaturated]
    
    while rejected_cells.size != 0:
        if verbose:
            logging.info("**{}, {}, {}".format(len(rejected_cells), len(unsaturated_cells), time.time()-ti))
        
        np.random.shuffle(rejected_cells) # random order
        # do something to rejected cells and unsaturated cells
        # knn_ji # for each cell in j, its knn in i
        knn_ji = clst_utils.gen_knn_annoy_train_test(mat_norm_i.values[unsaturated_cells], # look for nearest neighbors in i 
                                                           mat_norm_j.values[rejected_cells], # for each row in j
                                                           min(knn*knn_speed_factor, len(unsaturated_cells)), #  
                                                           form='list', # adj matrix 
                                                           metric=metric, # correlation 
                                                           n_trees=10, search_k=-1, verbose=False, 
                                                           include_distances=False, # for now
                                                           ).astype(int)
        knn_ji = unsaturated_cells[knn_ji] # transform it to global index, need to check this like 
        
        rejected_local_idx = []
        # examine each cell in j
        for local_idx, cell in enumerate(rejected_cells):
            # get knn in i
            knn_in_i = knn_ji[local_idx]
            # filter out saturated ones
            knn_in_i = knn_in_i[unsaturated[knn_in_i]]
            
            if knn_in_i.size < knn:
                # reject
                rejected_local_idx.append(local_idx)
            else:
                # accept and update
                accepted_knn_ji.append(knn_in_i[:knn])
                accepted_cells.append(cell)
                n_connects[knn_in_i[:knn]] += 1 
                unsaturated = (n_connects < k_saturate) # unsaturated bool 
                
        unsaturated_cells = np.arange(len(cells_i))[unsaturated]
        rejected_cells = rejected_cells[rejected_local_idx]
    # break
                
    accepted_knn_ji = pd.DataFrame(np.vstack(accepted_knn_ji), index=accepted_cells)
    accepted_knn_ji = accepted_knn_ji.sort_index().values

    
    return accepted_knn_ji
# 
def impute_1pair_cca(mod_i, mod_j, 
                     smoothed_features_i, smoothed_features_j,
                     settings,
                     knn,
                     relaxation,
                     n_cca,
                     output_knn_mat_ij='',
                     output_knn_mat_ji='',
                     impute_j=True,
                    ):
    """
    """
    # set up
    direct_i, direct_j = settings[mod_i].mod_direction, settings[mod_j].mod_direction
    
    mat_ii = smoothed_features_i.T # cell in mod i; gene in mod i
    mat_jj = smoothed_features_j.T # cell in mod j; gene in mod j
    
    genes_i = mat_ii.columns.values
    genes_j = mat_jj.columns.values
    genes_common = np.intersect1d(genes_i, genes_j)
    
    cells_i = mat_ii.index.values
    cells_j = mat_jj.index.values
    
    ## CCA euclidean distance 
    # normalize the feature matrix
    X = mat_ii[genes_common].T.apply(basic_utils.zscore, axis=0)*direct_i # gene by cell, zscore across genes
    Y = mat_jj[genes_common].T.apply(basic_utils.zscore, axis=0)*direct_j
    U, s, Vt = fbpca.pca(X.T.values.dot(Y.values), k=n_cca)
    del X, Y

    mat_norm_i = pd.DataFrame(U, index=mat_ii.index)
    maxk_i = int((len(cells_j)/len(cells_i))*knn*relaxation)+1 # max number of NN a cell in i can get 
    mat_norm_j = pd.DataFrame(Vt.T, index=mat_jj.index)
    maxk_j = int((len(cells_i)/len(cells_j))*knn*relaxation)+1 # max number of NN a cell in j can get 
    
    if impute_j:
        # knn_i and knn_j
        # j <- i for each j, get kNN in i
        knn_ji = get_constrained_knn(mat_norm_j, mat_norm_i, knn=knn, k_saturate=maxk_i, metric='euclidean')
        mat_knn_ji = sparse_adj_to_mat(knn_ji, len(cells_j), len(cells_i))
        
        if output_knn_mat_ji:
            sparse.save_npz(output_knn_mat_ji, mat_knn_ji)
        
        # normalize 
        degrees_j = np.ravel(mat_knn_ji.sum(axis=1)) # for each cell in j, how many cells in i it connects to 
        mat_knn_ji = sparse.diags(1.0/(degrees_j+1e-7)).dot(mat_knn_ji) 
        
        # imputation both across and within modality
        mat_ji = mat_knn_ji.dot(mat_ii) # cell in mod j, gene in mod i
    
    
    # i <- j
    knn_ij = get_constrained_knn(mat_norm_i, mat_norm_j, knn=knn, k_saturate=maxk_j, metric='euclidean')
    mat_knn_ij = sparse_adj_to_mat(knn_ij, len(cells_i), len(cells_j))
    
    if output_knn_mat_ij:
        sparse.save_npz(output_knn_mat_ij, mat_knn_ij)
    
    degrees_i = np.ravel(mat_knn_ij.sum(axis=1)) # for each cell in i, how many cells in j it connects to 
    mat_knn_ij = sparse.diags(1.0/(degrees_i+1e-7)).dot(mat_knn_ij) 
    
    mat_ij = mat_knn_ij.dot(mat_jj) # cell in mod i, gene in mod j
    
    if impute_j:
        return mat_ij, mat_ji
    else:
        return mat_ij

def impute_1pair(mod_i, mod_j, 
                 smoothed_features_i, smoothed_features_j,
                 settings,
                 knn, # 20
                 relaxation, # 3
                 output_knn_mat_ij='',
                 output_knn_mat_ji='',
                 impute_j=True,
                 ):
    """
    """
    # set up
    direct_i, direct_j = settings[mod_i].mod_direction, settings[mod_j].mod_direction
    
    mat_ii = smoothed_features_i.T # cell in mod i; gene in mod i
    mat_jj = smoothed_features_j.T # cell in mod j; gene in mod j
    
    genes_i = mat_ii.columns.values
    genes_j = mat_jj.columns.values
    genes_common = np.intersect1d(genes_i, genes_j)
    
    cells_i = mat_ii.index.values
    cells_j = mat_jj.index.values
    
    ## spearman correlation as distance  (rank -> zscore -> (flip sign?) -> "dot" distance) 
    # normalize the feature matrix
    mat_norm_i = (mat_ii[genes_common].rank(pct=True, axis=1)
                                      .apply(basic_utils.zscore, axis=1)
                                      *direct_i
                 )
    mat_norm_j = (mat_jj[genes_common].rank(pct=True, axis=1)
                                      .apply(basic_utils.zscore, axis=1)
                                      *direct_j
                 )
    maxk_i = int((len(cells_j)/len(cells_i))*knn*relaxation)+1 # max number of NN a cell in i can get 
    maxk_j = int((len(cells_i)/len(cells_j))*knn*relaxation)+1 # max number of NN a cell in j can get 
    
    if impute_j:
        # knn_i and knn_j
        # j <- i for each j, get kNN in i
        knn_ji = get_constrained_knn(mat_norm_j, mat_norm_i, knn=knn, k_saturate=maxk_i, metric='dot')
        mat_knn_ji = sparse_adj_to_mat(knn_ji, len(cells_j), len(cells_i))
        
        if output_knn_mat_ji:
            sparse.save_npz(output_knn_mat_ji, mat_knn_ji)
        
        # normalize 
        degrees_j = np.ravel(mat_knn_ji.sum(axis=1)) # for each cell in j, how many cells in i it connects to 
        mat_knn_ji = sparse.diags(1.0/(degrees_j+1e-7)).dot(mat_knn_ji) 
        
        # imputation both across and within modality
        mat_ji = mat_knn_ji.dot(mat_ii) # cell in mod j, gene in mod i
    
    
    # i <- j
    knn_ij = get_constrained_knn(mat_norm_i, mat_norm_j, knn=knn, k_saturate=maxk_j, metric='dot')
    mat_knn_ij = sparse_adj_to_mat(knn_ij, len(cells_i), len(cells_j))

    if output_knn_mat_ij:
        sparse.save_npz(output_knn_mat_ij, mat_knn_ij)
    
    degrees_i = np.ravel(mat_knn_ij.sum(axis=1)) # for each cell in i, how many cells in j it connects to 
    mat_knn_ij = sparse.diags(1.0/(degrees_i+1e-7)).dot(mat_knn_ij) 
    
    mat_ij = mat_knn_ij.dot(mat_jj) # cell in mod i, gene in mod j
    
    if impute_j:
        return mat_ij, mat_ji
    else:
        return mat_ij

def core_scf_routine(mods_selected, features_selected, settings, 
                    metas, gxc_hvftrs, 
                    ps, drop_npcs,
                    cross_mod_distance_measure, knn, relaxation, n_cca,
                    npc,
                    output_pcX_all, output_cells_all,
                    output_imputed_data_format,
                    ka_smooth=5,
                    save_knn=False,
                    output_knn_within="",
                    output_knn_across="",
                    ):
    """smooth within modality, impute across modalities, and construct a joint PC matrix
    """
    # GENE * CELL !!!!
    smoothed_features = collections.OrderedDict()
    logging.info("Smoothing within modalities...")
    for mod in mods_selected:
        ti = time.time()
        if settings[mod].mod_category == 'mc':
            _df = gxc_hvftrs[mod]
        else:
            _mat = gxc_hvftrs[mod].data.todense()
            _df = pd.DataFrame(_mat, 
                              index=gxc_hvftrs[mod].gene, 
                              columns=gxc_hvftrs[mod].cell, 
                              ) 
        if len(metas[mod]) < 2:
            raise ValueError("Dataset too small (contains <2 cells)")

        npc = min(len(metas[mod]), npc)
        # added 08/03/2020
        ka = min(len(metas[mod]), ka_smooth) 
        k_smooth = min(len(metas[mod]), 3*ka)

        if save_knn: # added 07/27/2020
            output_knn_mat = output_knn_within.format(mod)
        else:
            output_knn_mat = ""

        # if k_smooth >= 2*ka:
        mat_smoothed, mat_knn = smooth_in_modality(_df, _df, k=k_smooth, ka=ka, npc=npc, 
                                                     p=ps[settings[mod].mod_category], 
                                                     drop_npc=drop_npcs[settings[mod].mod_category], 
                                                     output_filename=output_knn_mat,
                                                     )
        smoothed_features[mod] = mat_smoothed
        # else:
        #     smoothed_features[mod] = _df

        logging.info("{}: {}".format(mod, time.time()-ti))
    # delete
    del gxc_hvftrs[mod]

    # construct a joint matrix (PCA)
    logging.info("Constructing a joint matrix...")
    cells_all = np.hstack([metas[mod].index.values for mod in mods_selected]) # cell (all mods)  
    pcX_all = []
    for mod_y in features_selected: ## to 
        logging.info("Imputing into {} space...".format(mod_y))
        # get all_features
        X = []
        for mod_x in mods_selected:
            logging.info("for {} cells...".format(mod_x))
            if mod_x == mod_y:
                smoothed_yy = smoothed_features[mod_y].T # gene by cell !!! VERY IMPORTANT
                X.append(smoothed_yy)
            else:
                if save_knn: # added 07/27/2020
                    output_knn_mat_ij = output_knn_across.format(mod_x, mod_y)
                    output_knn_mat_ji = output_knn_across.format(mod_y, mod_x)
                else:
                    output_knn_mat_ij = "" 
                    output_knn_mat_ji = "" 
                # impute x cells y space
                smoothed_features_x = smoothed_features[mod_x]
                smoothed_features_y = smoothed_features[mod_y]
                if cross_mod_distance_measure == 'correlation':
                    imputed_xy = impute_1pair(mod_x, mod_y, 
                                              smoothed_features_x, smoothed_features_y,
                                              settings,
                                              knn=knn,
                                              relaxation=relaxation,
                                              impute_j=False,
                                              output_knn_mat_ij=output_knn_mat_ij,
                                              output_knn_mat_ji=output_knn_mat_ji,
                                              )
                elif cross_mod_distance_measure == 'cca':
                    imputed_xy = impute_1pair_cca(mod_x, mod_y, 
                                                 smoothed_features_x, smoothed_features_y,
                                                 settings,
                                                 knn=knn,
                                                 relaxation=relaxation,
                                                 n_cca=n_cca,
                                                 impute_j=False,
                                                 output_knn_mat_ij=output_knn_mat_ij,
                                                 output_knn_mat_ji=output_knn_mat_ji,
                                                )
                else:
                    raise ValueError("Choose from correlation and cca")
                X.append(imputed_xy)
        X = np.vstack(X) # cell (all mods) by gene (mod_y) 
        # save X (imputed counts)
        np.save(output_imputed_data_format.format(mod_y), X)
        # PCA
        U, s, V = fbpca.pca(X, npc)
        del X
        pcX = U.dot(np.diag(s))
        # normalize PCs
        sigma = np.sqrt(np.sum(s*s)/(pcX.shape[0]*pcX.shape[1]))
        pcX = pcX/sigma
        pcX_all.append(pcX)
        
    pcX_all = np.hstack(pcX_all)
    # save pcX_all
    np.save(output_pcX_all, pcX_all)
    np.save(output_cells_all, cells_all)
    logging.info("Saved output to: {}".format(output_pcX_all))
    logging.info("Saved output to: {}".format(output_cells_all))
    return pcX_all, cells_all

def clustering_umap_routine(pcX_all, cells_all, mods_selected, metas, 
                            resolutions, k, 
                            umap_neighbors, min_dist, 
                            output_clst_and_umap,
                            cluster_only=False,
                            ):
    """
    """

    # clustering
    df_clsts = []
    for resolution in resolutions:
        logging.info('resolution r: {}'.format(resolution))
        df_clst = clst_utils.clustering_routine(
                                        pcX_all, 
                                        cells_all, k, 
                                        resolution=resolution,
                                        metric='euclidean', option='plain', n_trees=10, search_k=-1, verbose=False)
        df_clsts.append(df_clst.rename(columns={'cluster': 
                                                'cluster_joint_r{}'.format(resolution)
                                               }))
    df_clst = pd.concat(df_clsts, axis=1) 

    if cluster_only:
        # save results
        df_summary = df_clst
        df_summary['modality'] = ''
        for mod in mods_selected:
            _cells = metas[mod].index.values
            df_summary.loc[_cells, 'modality'] = mod
        df_summary.to_csv(output_clst_and_umap, sep='\t', header=True, index=True)    
        return df_summary

    else:
        # umap
        df_tsne = clst_utils.run_umap_lite(
                    pcX_all, 
                    cells_all, n_neighbors=umap_neighbors, min_dist=min_dist, n_dim=2, 
                    random_state=1)

        # summary
        df_summary = df_clst.join(df_tsne)
        df_summary['modality'] = ''
        for mod in mods_selected:
            _cells = metas[mod].index.values
            df_summary.loc[_cells, 'modality'] = mod
        # save results
        df_summary.to_csv(output_clst_and_umap, sep='\t', header=True, index=True)    
        return df_summary
