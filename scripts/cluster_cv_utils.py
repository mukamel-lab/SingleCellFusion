# import sys
# sys.path.insert(0, '/cndd/fangming/CEMBA/snmcseq_dev')

from __init__ import *
# from __init__jupyterlab import *

from scipy import sparse
import collections
import itertools
import re
# import fbpca
# from sklearn.model_selection import KFold

import basic_utils
# import CEMBA_run_tsne
# import CEMBA_clst_utils
# import SCF_utils

# general utils
def split_genes(gene_chrom_lookup, split_frac=0.5, report_chroms=False):
    """Gene_chrom_lookup: pd.Series gene->chrom (int, 1-20)
    """
    gene_chrom_nums = gene_chrom_lookup.value_counts()
    chrom_order = np.random.permutation(gene_chrom_nums.index.values)
    gene_chrom_nums_cumsum = np.cumsum(gene_chrom_nums.loc[chrom_order].values)
    cond = gene_chrom_nums_cumsum > gene_chrom_nums_cumsum[-1]*split_frac # if this condition is not satisfied, go to dataset 0, which account for split_frac of features
    chrom_set_lookup = {chrom: bl.astype(int) for chrom, bl in zip(chrom_order, cond)}
    gene_set_lookup = gene_chrom_lookup.apply(lambda x: chrom_set_lookup[x])

    chrom_order_0 = chrom_order[~cond]
    chrom_order_1 = chrom_order[cond]
    logging.info('Cluster chrom: {}, feature chrom {}'.format(chrom_order_0, chrom_order_1))
    
    if not report_chroms:
        return gene_set_lookup
    else:
        return gene_set_lookup, chrom_order_0, chrom_order_1


def split_features_routine_lite(data_matrix, gene_set_lookup):
    """
    - data_matrix: DataFrame, gene by cell
    - gene_set_lookup: Series, indexed by gene, binary lookup (0, 1)
    
    Returns:
    - sub_g0, sub_g1: DataFrame, gene by cell
    """
    _genes = data_matrix.index.values
    _lookup = gene_set_lookup.reindex(_genes).fillna(-1).values
    _genes_set0 = _genes[_lookup == 0]
    _genes_set1 = _genes[_lookup == 1]
    data_matrix_sub_g0 = data_matrix.loc[_genes_set0]
    data_matrix_sub_g1 = data_matrix.loc[_genes_set1]

    logging.info('Finish split: {} -> {} {}'.format(data_matrix.shape, 
                                                      data_matrix_sub_g0.shape, data_matrix_sub_g1.shape))
    return data_matrix_sub_g0, data_matrix_sub_g1

def split_features_routine(mods_selected, settings, gxc_hvftrs_sub, gene_set_lookup):
    """Do many datasets at the same time
    """
    gxc_hvftrs_sub_g0 = collections.OrderedDict()
    gxc_hvftrs_sub_g1 = collections.OrderedDict()
    for mod in mods_selected: 
        # split gxc_hvftrs
        if settings[mod].mod_category == 'mc':
            _genes = gxc_hvftrs_sub[mod].index.values
            _lookup = gene_set_lookup.reindex(_genes).fillna(-1).values
            _genes_set0 = _genes[_lookup == 0]
            _genes_set1 = _genes[_lookup == 1]
            print(len(_genes_set0), len(_genes_set1))
            gxc_hvftrs_sub_g0[mod] = gxc_hvftrs_sub[mod].loc[_genes_set0]
            gxc_hvftrs_sub_g1[mod] = gxc_hvftrs_sub[mod].loc[_genes_set1]

            logging.info('Finish split {}: {} -> {} {}'.format(mod, gxc_hvftrs_sub[mod].shape, 
                                                               gxc_hvftrs_sub_g0[mod].shape, gxc_hvftrs_sub_g1[mod].shape))
            continue

        _genes = gxc_hvftrs_sub[mod].gene
        _lookup = gene_set_lookup.reindex(_genes).fillna(-1).values
        _genes_set0 = _genes[_lookup == 0]
        _genes_set0_index = basic_utils.get_index_from_array(_genes, _genes_set0)
        _genes_set1 = _genes[_lookup == 1]
        _genes_set1_index = basic_utils.get_index_from_array(_genes, _genes_set1)
        print(len(_genes_set0), len(_genes_set1))
        gxc_hvftrs_sub_g0[mod] = GC_matrix(
                                    _genes_set0,
                                    gxc_hvftrs_sub[mod].cell,
                                    gxc_hvftrs_sub[mod].data.tocsr()[_genes_set0_index,:],
                                    )
        gxc_hvftrs_sub_g1[mod] = GC_matrix(
                                    _genes_set1,
                                    gxc_hvftrs_sub[mod].cell,
                                    gxc_hvftrs_sub[mod].data.tocsr()[_genes_set1_index,:],
                                    )

        logging.info('Finish split {}: {} -> {} {}'.format(mod, gxc_hvftrs_sub[mod].data.shape, 
                                                           gxc_hvftrs_sub_g0[mod].data.shape, gxc_hvftrs_sub_g1[mod].data.shape))
        
    return gxc_hvftrs_sub_g0, gxc_hvftrs_sub_g1

def subsampling_lite(metadata, data_matrix, p=0, n=0):
    """Do many datasets at the same time
    p - fraction of cells from each dataset to be included
    """
    # subsample meta
    if p != 0: 
        cells_included = metadata.index.values[np.random.rand(len(metadata))<p]
    elif n != 0:
        if n < len(metadata):
            cells_included = metadata.sample(n=n).index.values
        else:
            cells_included = metadata.index.values
    metadata_sub = metadata.loc[cells_included]
    
    # subsample data_matrix
    data_matrix_sub = data_matrix[cells_included]
    print(metadata_sub.shape, data_matrix_sub.shape)
    return metadata_sub, data_matrix_sub

def subsampling(mods_selected, settings, metas, gxc_hvftrs, p=0, n=0):
    """Do many datasets at the same time
    p - fraction of cells from each dataset to be included
    """
    metas_sub = collections.OrderedDict()
    gxc_hvftrs_sub = collections.OrderedDict()
    for mod in mods_selected: 
        # subsample meta
        if p != 0: 
            cells_included = metas[mod].index.values[np.random.rand(len(metas[mod]))<p]
        elif n != 0:
            if n > len(metas[mod]):
                n = len(metas[mod])
                
            cells_included = metas[mod].sample(n=n).index.values
        metas_sub[mod] = metas[mod].loc[cells_included]

        # subsample gxc_hvftrs
        if settings[mod].mod_category == 'mc':
            gxc_hvftrs_sub[mod] = gxc_hvftrs[mod][cells_included]
            logging.info("{} {} {}".format(mod, metas_sub[mod].shape, gxc_hvftrs_sub[mod].shape))
            continue

        cells_included_idx = basic_utils.get_index_from_array(gxc_hvftrs[mod].cell, cells_included)
        gxc_hvftrs_sub[mod] = GC_matrix(
                                        gxc_hvftrs[mod].gene,
                                        cells_included,
                                        gxc_hvftrs[mod].data.tocsc()[:, cells_included_idx],
                                        )
        logging.info("{} {} {}".format(mod, metas_sub[mod].shape, gxc_hvftrs_sub[mod].data.shape))
    return metas_sub, gxc_hvftrs_sub


# # cross validation utils
# def centroid_prediction_model_cv(metadata, features_y, n_repeats, nfolds):
#     """Prediction of y with centroids (x is the cluster label), return training and test MSEs (cross-validated) for n_repeats 
#     - metadata: DataFrame, indexed by cell_name, "cluster_cv" as a column
#     - features_y: Numpy array, Feature matrix cell by feature
#     - n_repeats: number of shuffles on cells
#     - nfolds: k fold cross-validation
#     """
#     # Assuming metadata and features_y has the same order of cells
#     assert isinstance(features_y, np.ndarray)
    
#     kl = KFold(n_splits=nfolds)
#     ncells = len(metadata)
#     res = []
#     cells_clst = metadata['cluster_cv']
    
#     for i_repeat in range(n_repeats):
#         print('.', end='')
#         # shuffle data
#         cells_shuffled_idx = np.random.permutation(np.arange(ncells))
#         metadata_shuffled = metadata.iloc[cells_shuffled_idx, :].copy()
#         metadata_shuffled['cell_idx'] = np.arange(ncells)
#         features_y_shuffled = features_y[cells_shuffled_idx, :]

#         # split training and test 
#         for i_fold, (train_idx, test_idx) in enumerate(kl.split(np.arange(ncells))):
#             ti = time.time()
#             # compute cluster centroids for training cells 
#             metadata_train = metadata_shuffled.iloc[train_idx]
#             clsts_in_train = np.unique(metadata_train['cluster_cv'].values)
#             clsts_not_in_train = np.unique(metadata_shuffled['cluster_cv'].values).tolist()
#             y_centroids = np.zeros((len(clsts_in_train), features_y_shuffled.shape[1]))
#             cluster_to_idx_lookup = {}
#             for count_idx, (clst, df_sub) in enumerate(metadata_train.groupby('cluster_cv')):
#                 cells_sub_idx = df_sub['cell_idx'].values
#                 y_centroids[count_idx, :] = features_y_shuffled[cells_sub_idx, :].mean(axis=0)
#                 cluster_to_idx_lookup[clst] = count_idx
#                 clsts_not_in_train.remove(clst)
#             for clst in clsts_not_in_train:
#                 cluster_to_idx_lookup[clst] = -1

#             # compute MSE for test cells
#             cells_j = metadata_shuffled.index.values[test_idx]
#             clsts_i = cells_clst[cells_j]
#             clsts_i_idx = np.array([cluster_to_idx_lookup[clst] for clst in clsts_i])
#             cond = (clsts_i_idx != -1)  # test if clsts_i in clsts_in_train
#             test_idx, cells_j, clsts_i, clsts_i_idx = test_idx[cond], cells_j[cond], clsts_i[cond], clsts_i_idx[cond]
#             diff = features_y_shuffled[test_idx, :] - y_centroids[clsts_i_idx, :]
#             mse = (diff**2).sum(axis=1).mean()

#             # compute MSE for training cells 
#             cells_j = metadata_train.index.values
#             clsts_i = cells_clst[cells_j]
#             clsts_i_idx = np.array([cluster_to_idx_lookup[clst] for clst in clsts_i])
#             diff = features_y_shuffled[train_idx, :] - y_centroids[clsts_i_idx, :]
#             mse_t = (diff*diff).sum(axis=1).mean()
            
#             # collect 1 data point
#             res.append({'mse': mse, 
#                         'mse_t': mse_t, 
#                         'i_fold': i_fold, 
#                         'i_repeat': i_repeat,
#                        })
#         # end of n-fold training test for 
#     # end of n-repeats for
#     res = pd.DataFrame(res)
#     return res 

# def nfoldcv_random_features_split(data_matrix, resolutions, gene_chrom_lookup,
#                                   output_prefix,
#                                   k=30, 
#                                   reduce_dim=0,
#                                   nfolds=5, n_repeats=5, n_splits=5, split_frac=0.8):
#     """cv for random features
#     - data_matrix: DataFrame gene by cell
#     """
#     cell_list = data_matrix.columns.values 
#     metadata = pd.DataFrame(index=cell_list)
    
#     res = []
#     res_nclsts = []
    
#     # n_splits (split and cluster)
#     for n_split in np.arange(n_splits):
#         gene_set_lookup = split_genes(gene_chrom_lookup, split_frac=split_frac)
#         data_matrix_sub_g0, data_matrix_sub_g1 = split_features_routine_lite(data_matrix, gene_set_lookup)
        
#         # PCA on g0
#         U, s, Vt = fbpca.pca(data_matrix_sub_g0.T.values, k=50) # transpose very important!
#         pcX_0 = U.dot(np.diag(s))
#         # cluster on PCA_g0 with different resolutions
#         df_clsts = CEMBA_clst_utils.clustering_routine_multiple_resolutions(pcX_0, cell_list, k, 
#                                                     seed=1, verbose=False,
#                                                     resolutions=resolutions, metric='euclidean', option='plain', n_trees=10, search_k=-1, num_starts=None)
#         output = "{}_{}.tsv".format(output_prefix, n_split) 
#         df_clsts.to_csv(output, sep="\t", header=True, index=True)
#         logging.info("Saved to {}".format(output))
            
#         # test different resolution
#         for resolution in resolutions:
#             print(resolution, end='')
        
#             df_clst = df_clsts[['cluster_r{}'.format(resolution)]].rename(columns={'cluster_r{}'.format(resolution): 'cluster'})
#             nclsts = len(df_clst['cluster'].unique())
#             res_nclsts.append({
#                         'i_split': n_split,
#                         'resolution': resolution,
#                         'nclsts': nclsts,
#                         }) # record number of clusters

#             # set up cross-validation on g1 (features_y)
#             metadata_copy = metadata.copy()
#             metadata_copy['cluster_cv'] = df_clst.loc[metadata_copy.index, 'cluster'] 
#             assert np.all(metadata_copy.index.values == data_matrix_sub_g1.columns.values)
#             features_y = data_matrix_sub_g1.T.values # transpose very important!
#             if reduce_dim:
#                 U, s, Vt = fbpca.pca(features_y, k=reduce_dim)
#                 features_y = U.dot(np.diag(s))
#             # one set CV
#             res_cv = centroid_prediction_model_cv(metadata_copy, features_y, n_repeats, nfolds)
#             res_cv['i_split'] = n_split 
#             res_cv['resolution'] =  resolution
#             res.append(res_cv)
            
#         # end of resolution for
#     # end of n-split for 
#     res = pd.concat(res, axis=0)
#     res_nclsts = pd.DataFrame(res_nclsts)
#     return res_nclsts, res 

# def nfoldcv_fixed_features_split(data_matrix_sub_g0, data_matrix_sub_g1, resolutions, 
#                                 k=30, 
#                                 reduce_dim=0,
#                                 n_pc=50,
#                                 nfolds=5, n_repeats=5, n_clstseeds=5):
#     """Cross validation for fixed features
#     - data_matrix: DataFrame gene by cell
#     """
#     cell_list = data_matrix_sub_g0.columns.values 
#     assert np.all(cell_list==data_matrix_sub_g1.columns.values)
#     metadata = pd.DataFrame(index=cell_list)
#     res = []
#     res_nclsts = []

#     # PCA on g0
#     U, s, Vt = fbpca.pca(data_matrix_sub_g0.T.values, k=n_pc) # transpose very important!
#     pcX_0 = U.dot(np.diag(s))
#     # cluster on PCA_g0 with different resolutions and seeds
#     df_clsts_all = []
#     for i_seed in range(n_clstseeds):
#         seed = int(np.random.rand()*1e5)
#         df_clsts = CEMBA_clst_utils.clustering_routine_multiple_resolutions(pcX_0, cell_list, k, 
#                                                     seed=seed, verbose=False,
#                                                     resolutions=resolutions, 
#                                                     metric='euclidean', option='plain', n_trees=10, search_k=-1, num_starts=None)
#         df_clsts_all.append(df_clsts)
        
#     # test different resolution
#     for resolution in resolutions:
#         print(resolution, end='')
    
#         for i_seed in range(n_clstseeds):
#             df_clst = (df_clsts_all[i_seed][['cluster_r{}'.format(resolution)]]
#                                             .rename(columns={'cluster_r{}'.format(resolution): 'cluster'})
#                       )
#             nclsts = len(df_clst['cluster'].unique())
#             res_nclsts.append({
#                         'resolution': resolution,
#                         'nclsts': nclsts,
#                         'i_seed': i_seed,
#                         }) # record number of clusters

#             # set up cross-validation on g1 (features_y)
#             metadata_copy = metadata.copy()
#             metadata_copy['cluster_cv'] = df_clst.loc[metadata_copy.index, 'cluster'] 
#             assert np.all(metadata_copy.index.values == data_matrix_sub_g1.columns.values)
#             features_y = data_matrix_sub_g1.T.values # transpose very important!
#             if reduce_dim:
#                 U, s, Vt = fbpca.pca(features_y, k=reduce_dim)
#                 features_y = U.dot(np.diag(s))
#             # one set CV
#             res_cv = centroid_prediction_model_cv(metadata_copy, features_y, n_repeats, nfolds)
#             res_cv['resolution'] =  resolution
#             res_cv['i_seed'] = i_seed 
#             res.append(res_cv)
#         # end of n_clstseeds for
        
#     # end of resolution for
        
#     res = pd.concat(res, axis=0)
#     res_nclsts = pd.DataFrame(res_nclsts)
#     return res_nclsts, res 


# def nfoldcv_scf_random_features_split(gxc_hvftrs, resolutions, gene_chrom_lookup,
#                                       mods_selected, metas, 
#                                       features_selected, settings, 
#                                       ps, drop_npcs,
#                                       cross_mod_distance_measure, knn, relaxation, n_cca,
#                                       npc,
#                                       output_pcX_all, output_cells_all, output_clst_and_umap_prefix,
#                                       output_imputed_data_format,
#                                       k=30, 
#                                       reduce_dim=0,
#                                       nfolds=5, n_repeats=5, n_splits=5, split_frac=0.8):
#     """Combining cv with SCF
#     - data_matrices: dict of DataFrame gene by cell
#     """
#     # check cell order
#     for mod in mods_selected:
#         metadata = metas[mod]
#         if settings[mod].mod_category == 'mc':
#             cell_array = gxc_hvftrs[mod].columns.values 
#         else:
#             cell_array = gxc_hvftrs[mod].cell 
#         assert np.all(metadata.index.values == cell_array)
    
#     res = []
#     res_nclsts = []
#     # n_splits (split and cluster)
#     for n_split in np.arange(n_splits):
#         # split features
#         gene_set_lookup = split_genes(gene_chrom_lookup, split_frac=split_frac)
#         gxc_hvftrs_sub_g0, gxc_hvftrs_sub_g1 = split_features_routine(mods_selected, settings, gxc_hvftrs, gene_set_lookup)
            
#         # run scf 
#         # cluster on g0 with different resolutions
#         pcX_all, cells_all = SCF_utils.core_scf_routine(mods_selected, features_selected, settings, 
#                                                         metas, gxc_hvftrs_sub_g0, 
#                                                         ps, drop_npcs,
#                                                         cross_mod_distance_measure, knn, relaxation, n_cca,
#                                                         npc,
#                                                         output_pcX_all, output_cells_all,
#                                                         output_imputed_data_format,
#                                                         )
#         umap_neighbors = 60
#         min_dist = 0.5
#         output_clst_and_umap = "{}_{}.tsv".format(output_clst_and_umap_prefix, n_split)
#         df_clsts = SCF_utils.clustering_umap_routine(pcX_all, cells_all, mods_selected, metas,
#                                                      resolutions, k, 
#                                                      umap_neighbors, min_dist,
#                                                      output_clst_and_umap,
#                                                      cluster_only=True,
#                                                      )
    
#         # test different resolution
#         for resolution in resolutions:
#             print(resolution, end='')
#             df_clst = (df_clsts[['cluster_joint_r{}'.format(resolution)]]
#                                 .copy()
#                                 .rename(columns={'cluster_joint_r{}'.format(resolution): 'cluster'})
#                       )
#             nclsts = len(df_clst['cluster'].unique())
#             res_nclsts.append({
#                         'i_split': n_split,
#                         'resolution': resolution,
#                         'nclsts': nclsts,
#                         }) # record number of clusters

#             # set up cross-validation on g1 (features_y) with every mod
#             for mod in mods_selected:
#                 metadata_copy = metas[mod].copy()
#                 metadata_copy['cluster_cv'] = df_clst.loc[metadata_copy.index, 'cluster'] 
#                 gxc_hvftr_sub_g1 = gxc_hvftrs_sub_g1[mod]
                
#                 if settings[mod].mod_category == 'mc':
#                     assert np.all(metadata_copy.index.values == gxc_hvftr_sub_g1.columns.values)
#                     features_y = gxc_hvftr_sub_g1.T.values # transpose very important!
#                 else:
#                     assert np.all(metadata_copy.index.values == gxc_hvftr_sub_g1.cell)
#                     features_y = np.asarray(gxc_hvftr_sub_g1.data.todense().T) # transpose very important, turn into a ndarray is super important!
#                 if reduce_dim:
#                     U, s, Vt = fbpca.pca(features_y, k=reduce_dim)
#                     features_y = U.dot(np.diag(s))
                    
#                 # one set CV
#                 res_cv = centroid_prediction_model_cv(metadata_copy, features_y, n_repeats, nfolds)
#                 res_cv['mod'] = mod
#                 res_cv['i_split'] = n_split 
#                 res_cv['resolution'] =  resolution
#                 res.append(res_cv)
#             # end of modality for 
#         # end of resolution for
#     # end of n-split for 
#     res = pd.concat(res, axis=0)
#     res_nclsts = pd.DataFrame(res_nclsts)
#     return res_nclsts, res 

# # plotting utils
# def plot_errorbar_ax(ax, x, y, yerr, color='C0', label=''):
#     """Plot a line with errorbar 
#     """
#     x = np.array(x)
#     y = np.array(y)
#     yerr = np.array(yerr)
    
#     ax.plot(x, y, '-o', 
#            markersize=5,
#            color=color,
#            label=label,
#            )
#     ax.fill_between(x, y-yerr, y+yerr, 
#                     color=color,
#                     alpha=0.3,
#                     zorder=0,
#                    )
#     return

# def plot_errorbar_fancymin_ax(ax, x, y, yerr, color='C0', label=''):
#     """Plot a line with errorbar + min position and min-se position
#     """
#     from scipy import optimize
#     x = np.array(x)
#     y = np.array(y)
#     yerr = np.array(yerr)
#     plot_errorbar_ax(ax, x, y, yerr, color=color, label=label)
    
#     # get minimum and plot
#     min_arg = np.nanargmin(y) # update 11/20/2019
#     min_x = x[min_arg]
#     min_y = y[min_arg]
#     ax.plot(min_x, min_y, '^',
#                markersize=12,
#                color=color,
#                )
    
#     # get minimum + se and plot
#     epsilon = 0
#     if min_arg == 0:
#         min_x_se = min_x

#     else:
#         f = lambda _x: np.interp(_x, x[:min_arg], (y-yerr)[:min_arg]) - (min_y + epsilon)
#         try:
#             res_root = optimize.root_scalar(f, bracket=(x[0], min_x))
#             min_x_se = int(res_root.root+0.5)
#         except: 
#             if np.all(f(x[:min_arg])<0):
#                 min_x_se = x[0]
#             elif np.all(f(x[:min_arg])>0):
#                 min_x_se = x[min_arg]
#             else:
#                 raise ValueError("Dont understand f: {}".format(f(x[:min_arg])))
#         ax.plot(min_x_se, min_y, 's', 
#                    markersize=10,
#                    color=color,
#                )
        
#     return int(min_x_se), int(min_x), min_y
    
# def plot_bi_cv_ax(ax, x, y, yerr, color='C0', mod="", ylabel="MSE +/- SEM Normalized"):
#     """
#     """
#     min_x_se, min_x, min_y = plot_errorbar_fancymin_ax(ax, x, y, yerr, color=color,)
#     ax.set_title("{}: {} - {}".format(mod, min_x_se, min_x))
#     ax.set_ylabel(ylabel)
#     return min_x_se, min_x, min_y

# def plot_bi_cv_subfig(ax, x1, y1, yerr1, y1_tr, yerr1_tr, color1, mod1,
#                     xlabel='Number of clusters',
#                     ylabel='MSE +/- SEM Normalized',
#                    ):
#     from matplotlib.ticker import ScalarFormatter
    
#     plot_errorbar_ax(ax, x1, y1_tr, yerr1_tr, color='black', label='Training error')
#     min_x_se, min_x, min_y = plot_bi_cv_ax(ax, x1, y1, yerr1, color=color1, mod=mod1, ylabel=ylabel)
        
#     ax.set_xscale('log', basex=2)
#     ax.xaxis.set_major_formatter(ScalarFormatter())
#     ax.set_xlabel(xlabel)
    
#     return min_x_se, min_x, min_y

# # tmp
# def plot_bi_cv_fig_mCT(x1, y1, yerr1, y1_tr, yerr1_tr, color1, mod1,
#                     x2, y2, yerr2, y2_tr, yerr2_tr, color2, mod2,
#                     output, 
#                     xlabel='Number of clusters',
#                    ):
#     from matplotlib.ticker import ScalarFormatter
    
#     fig, axs = plt.subplots(2, 1, figsize=(6, 4*2), sharex=True)
#     ax = axs[0]
#     plot_errorbar_ax(ax, x1, y1_tr, yerr1_tr, color='black', label='Training error')
#     plot_bi_cv_ax(ax, x1, y1, yerr1, color=color1, mod=mod1)
        
#     ax = axs[1]
#     plot_errorbar_ax(ax, x2, y2_tr, yerr2_tr, color='black', label='Training error')
#     plot_bi_cv_ax(ax, x2, y2, yerr2, color=color2, mod=mod2)
#     ax.set_xscale('log', basex=2)
#     ax.xaxis.set_major_formatter(ScalarFormatter())
    
#     ax.set_xlabel(xlabel)
    
#     fig.tight_layout()
#     fig.savefig(output, bbox_inches='tight')
#     return
