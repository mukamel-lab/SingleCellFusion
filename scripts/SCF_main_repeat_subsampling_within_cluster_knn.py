#!/usr/bin/env python3
"""SingleCellFusion main rontine"""

from __init__ import *

from scipy import sparse
import collections
import itertools
import sys
import pickle
import argparse
import glob

import basic_utils
import SCF_utils


def create_parser():
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config_py", help="Configuration file", required=True)
    parser.add_argument("-s", "--subsample_fraction", help="0~1", type=float, required=True)
    parser.add_argument("-sn", "--subsample_times", help=">1", type=int, required=True)
    parser.add_argument("-f", "--fname_cluster", help="Cluster file", type=str, required=True)
    parser.add_argument("-col", "--column_cluster", help="the cluster column", type=str, required=True)
    return parser

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

def get_oneclst(mods_selected, settings, metas, gxc_hvftrs, df_clst_oneclst, num_threshold=0):
    """Do many datasets at the same time
    p - fraction of cells from each dataset to be included
    """
    metas_sub = collections.OrderedDict()
    gxc_hvftrs_sub = collections.OrderedDict()
    for mod in mods_selected: 
        # subsample meta
        cells_included = np.intersect1d(df_clst_oneclst.index.values, metas[mod].index.values)
        if len(cells_included) < num_threshold:
            # one modality doesn't meet the requirement
            return None, None, False 
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

    # all modalities meet the requirement
    return metas_sub, gxc_hvftrs_sub, True


log = basic_utils.create_logger()


parser = create_parser()
args = parser.parse_args()

config_dirc, config_py = os.path.split(args.config_py)
logging.info("{} {}".format(config_dirc, config_py))
if config_py.endswith('.py'):
    config_py = config_py[:-3]
if os.path.isdir(config_dirc):
    logging.info('Adding {} to python path'.format(config_dirc))
    sys.path.insert(0, config_dirc)
exec("from {} import *".format(config_py))

if not os.path.isdir(outdir):
    os.makedirs(outdir)
# end of configurations

### ---- fixed after ----
# ## Read in data 
logging.info('* Begin integration')

metas = collections.OrderedDict()
for mod in mods_selected:
    metas[mod] = pd.read_csv(meta_f.format(mod), sep="\t").reset_index().set_index(settings[mod].cell_col)
    logging.info("Metadata {} {}".format(mod, metas[mod].shape))

gxc_hvftrs = collections.OrderedDict()
for mod in mods_selected:
    if settings[mod].mod_category == 'mc':
        f_mat = hvftrs_f.format(mod, 'tsv')
        gxc_hvftrs[mod] = pd.read_csv(f_mat, sep='\t', header=0, index_col=0) 
        logging.info("Feature matrix {} {}".format(mod, gxc_hvftrs[mod].shape))
        assert np.all(gxc_hvftrs[mod].columns.values == metas[mod].index.values) # make sure cell name is in the sanme order as metas (important if save knn mat)
        continue
        
    f_mat = hvftrs_f.format(mod, 'npz')
    f_gene = hvftrs_gene.format(mod)
    f_cell = hvftrs_cell.format(mod)
    _gxc_tmp = basic_utils.load_gc_matrix(f_gene, f_cell, f_mat)
    _gene = _gxc_tmp.gene
    _cell = _gxc_tmp.cell
    _mat = _gxc_tmp.data

    gxc_hvftrs[mod] = GC_matrix(_gene, _cell, _mat)
    assert np.all(gxc_hvftrs[mod].cell == metas[mod].index.values) # make sure cell name is in the sanme order as metas (important if save knn mat)
    logging.info("Feature matrix {} {}".format(mod, gxc_hvftrs[mod].data.shape))
logging.info('Done reading data')


subsample_fraction = args.subsample_fraction
subsample_times = args.subsample_times
fname_cluster = args.fname_cluster
column_cluster = args.column_cluster

df_clst = pd.read_csv(fname_cluster, sep='\t', index_col=0)[[column_cluster]] # cluster table
uniq_clusters = np.sort(df_clst[column_cluster].unique())
# subsampling 80% of cells from each modality
for i in range(subsample_times):
    logging.info("Subsampling {}/{}".format(i+1, subsample_times))
    output_tag = ".{}".format(i)
    metas_sub, gxc_hvftrs_sub = subsampling(
                mods_selected, settings, metas, gxc_hvftrs, p=subsample_fraction,
                )

    # run SCF for each cluster then patch the results together
    knn_xx_all = []
    knn_xy_all = []
    cells_ordered = {mod: [] for mod in mods_selected}
    for clst in uniq_clusters:
        logging.info(clst)
        # get cells and further split matrices
        df_clst_oneclst = df_clst[df_clst[column_cluster]==clst]

        # prep new metas_sub, gxc_hvftrs_sub
        metas_sub_clst, gxc_hvftrs_sub_clst, enough_cells = get_oneclst(
                mods_selected, settings, metas_sub, gxc_hvftrs_sub, df_clst_oneclst,
                num_threshold=2, # should have at least 2 cells to proceed
                )
        if not enough_cells:
            logging.info("! not enough cell for cluster {}, skipped".format(clst))
            continue # skip this cluster 

        # record the new cell orders
        for mod in mods_selected:
            cells_ordered[mod].append(metas_sub_clst[mod].index.values)

        tmp_output_tag = output_tag+"__TMP__"+clst+"__"

        # ## run SCF
        pcX_all, cells_all = SCF_utils.core_scf_routine(mods_selected, features_selected, settings, 
                                                        metas_sub_clst, gxc_hvftrs_sub_clst, # sub 
                                                        ps, drop_npcs,
                                                        cross_mod_distance_measure, knn, relaxation, n_cca,
                                                        npc,
                                                        output_pcX_all+tmp_output_tag, output_cells_all+tmp_output_tag,
                                                        output_imputed_data_format+tmp_output_tag,
                                                        ka_smooth=ka_smooth,
                                                        save_knn=save_knn,
                                                        output_knn_within=output_knn_within+tmp_output_tag,
                                                        output_knn_across=output_knn_across+tmp_output_tag,
                                                        )
        # collect all 
        # according to the cell order already saved (cells_mod)
        mod_y = features_selected[0]
        mod_xs = [mod for mod in mods_selected if mod != mod_y]
        mod_x = mod_xs[0]
        # mod_y = 'snmcseq_gene'
        # mod_x = 'smarter_cells'
        # 20 cells in mod_y for each cell in mod_x

        fname = (output_knn_across+tmp_output_tag+'.npz').format(mod_x, mod_y)
        knn_xy_all.append(sparse.load_npz(fname))  

        fname = (output_knn_within+tmp_output_tag+'.npz').format(mod_x)
        knn_xx_all.append(sparse.load_npz(fname)) 

        # remove __TMP__ tags
        to_removes = glob.glob(outdir+"/*{}*".format(tmp_output_tag))
        logging.info("removing {} files: {}".format(len(to_removes), to_removes))
        for to_remove in to_removes:
            os.remove(to_remove) 

    # glue together clusters for each set of matrices
    # output this
    knn_xx_all = sparse.block_diag(knn_xx_all) # patch along the diagnals
    sparse.save_npz((output_knn_within+output_tag).format(mod_x), knn_xx_all)

    knn_xy_all = sparse.block_diag(knn_xy_all) # patch along the diagnals
    sparse.save_npz((output_knn_across+output_tag).format(mod_x, mod_y), knn_xy_all)

    # output cells in order
    for mod in mods_selected:
        cells_ordered[mod] = np.hstack(cells_ordered[mod])
        fout = (output_cells+output_tag).format(mod)
        np.save(fout, cells_ordered[mod])






