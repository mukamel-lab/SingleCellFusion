#!/usr/bin/env python3
"""SingleCellFusion main routine"""

from __init__ import *

from scipy import sparse
import collections
import itertools
import sys
import pickle
import argparse
import anndata

import basic_utils
import SCF_utils

log = basic_utils.create_logger()

def create_parser():
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config_py", help="Configuration file (full path)", required=True)
    return parser
parser = create_parser()
args = parser.parse_args()

# replace the config.py with argparse command-line api
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
# 


# ## Read in data 
logging.info('* Begin integration')

### read in data (h5ad)
metas = collections.OrderedDict()
gxc_hvftrs = collections.OrderedDict()
for mod in mods_selected:
    print(mod)
    if settings[mod].mod_category == 'mc':
        # path
        _file = data_f.format(mod)
        # read 
        print(_file)
        h5ad_mat = anndata.read_h5ad(_file) 
        # convert
        meta, mat = basic_utils.h5ad_to_scf_mc_format(h5ad_mat)
        
        metas[mod] = meta
        gxc_hvftrs[mod] = mat

        assert np.all(gxc_hvftrs[mod].columns.values == metas[mod].index.values) # make sure cell name is in the sanme order as metas (important if save knn mat)

        logging.info("Metadata {} {}".format(mod, metas[mod].shape))
        logging.info("Feature matrix {} {}".format(mod, gxc_hvftrs[mod].shape))
        continue
        
    # path 
    _file = data_f.format(mod)
    # read 
    print(_file)
    h5ad_mat = anndata.read_h5ad(_file) 
    # convert
    meta, gc_mat = basic_utils.h5ad_to_scf_rna_format(h5ad_mat)
    
    metas[mod] = meta
    gxc_hvftrs[mod] = gc_mat

    assert np.all(gxc_hvftrs[mod].cell == metas[mod].index.values) # make sure cell name is in the sanme order as metas (important if save knn mat)
    logging.info("Feature matrix {} {}".format(mod, gxc_hvftrs[mod].data.shape))
logging.info('Done reading data')


# ## run SCF
pcX_all, cells_all = SCF_utils.core_scf_routine(mods_selected, features_selected, settings, 
                                                metas, gxc_hvftrs, 
                                                ps, drop_npcs,
                                                cross_mod_distance_measure, knn, relaxation, n_cca,
                                                npc,
                                                output_pcX_all, output_cells_all,
                                                output_imputed_data_format,
                                                )
logging.info('Done integration into a common PC space')


df_summary = SCF_utils.clustering_umap_routine(pcX_all, cells_all, mods_selected, metas,
                                               resolutions, k, 
                                               umap_neighbors, min_dist, 
                                               output_clst_and_umap,
                                              )
logging.info('Done clustering and UMAP')
