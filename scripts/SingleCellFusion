#!/usr/bin/env python3
"""SingleCellFusion main routine"""

from __init__ import *
# public packages 
from scipy import sparse
import collections
import itertools
import sys
import os
import pickle
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)
    import anndata
import pandas
# scripts from this package
import cli_parser
import basic_utils
import SCF_utils

# TODO:
# -- Proper error messages (replace assert() with helpful messages)

parser = cli_parser.create_parser()
args = parser.parse_args()

log = basic_utils.create_logger()
logging.info('* Parsing Command Line Arguments')

# specify output filenames
outdir = args.output_dir
if not os.path.isdir(outdir):
    os.makedirs(outdir)

name = args.output_prefix
output_pcX_all = outdir + '/{}_pcX_all.npy'.format(name)
output_cells_all = outdir + '/{}_cells_all.npy'.format(name)
output_imputed_data_format = outdir + '/{}_imputed_data_{{}}.npy'.format(name)
output_clst_and_umap = outdir + '/{}_intg_summary.tsv'.format(name)
output_cluster_centroids = outdir + '/{}_centroids.pkl'.format(name)
output_figures = outdir + '/{}_{{}}.{{}}'.format(name)

# get input files, modaltiies (internal rep of input files), and feature datasets
data_files = args.input_datasets
feature_files = args.feature_datasets
mods_selected = [cli_parser.parse_filename(data_file) for data_file in data_files]
features_selected = [cli_parser.parse_filename(data_file) for data_file in feature_files]
for features_modality in features_selected:
    assert (features_modality in mods_selected)

# get dataset metadata
mod_catgories = args.input_modalities
for mod_category in mod_catgories:
    assert (mod_category in ['mc', 'atac', 'rna'])
settings = collections.OrderedDict()
Mod_info = collections.namedtuple('Mod_info', ['mod', 'mod_category', 'mod_direction',])
for mod, mod_category in zip(mods_selected, mod_catgories):
    mod_direction = cli_parser.modality_default_options(mod_category)
    settings[mod] = Mod_info(mod, mod_category, mod_direction,)

# parameters
# Within modality
ps = {
    'rna': args.smoothing_fractions[0],
    'atac': args.smoothing_fractions[1],
    'mc': args.smoothing_fractions[2],
}
# across modality
knn = args.nearest_neighbors
relaxation = args.relaxation
# PCA
npc = args.num_pcs
# clustering
k = args.leiden_n_neighbors
resolutions = args.leiden_resolutions
# umap
umap_neighbors = args.umap_n_neighbors
min_dist = args.umap_min_dist

### --- deprecated arguments (for testing; not open to general users)
n_cca = 0 # deprecated args.n_cca
drop_npcs = {
      'mc': 0, 
      'rna': 0, 
      'atac': 0, 
     } 
cross_mod_distance_measure = 'correlation' # or 'cca' 
### --- end of deprecation



# ## Read in data 
logging.info('* Begin integration')
### read in data (h5ad)
metas = collections.OrderedDict()
gxc_hvftrs = collections.OrderedDict()
for mod, _file in zip(mods_selected, data_files):
    print(mod)
    if settings[mod].mod_category == 'mc':
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
