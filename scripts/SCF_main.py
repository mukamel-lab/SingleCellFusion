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
import anndata
import pandas
# scripts from this package
import cli_parser
import basic_utils
import SCF_utils

# TODO:
# -- Clean up command line argument descriptions
# -- Figure out where to put default values
# -- Fix __init__dataset.py rept
# -- Proper error messages (replace assert() with helpful messages)
# -- Simplify API to take care of redundant information
# -- Make the assumed input structures explicit (./datasets/dataset_metadata.csv)

log = basic_utils.create_logger()

parser = cli_parser.create_parser()
args = parser.parse_args()

logging.info('* Parsing Command Line Arguments')

# Get input and output directories
DATA_DIR = args.dataset
outdir = args.output

if not os.path.isdir(outdir):
    os.makedirs(outdir)

# output file names
name = args.name

output_pcX_all = outdir + '/pcX_all_{}.npy'.format(name)
output_cells_all = outdir + '/cells_all_{}.npy'.format(name)
output_imputed_data_format = outdir + '/imputed_data_{}_{{}}.npy'.format(name)
output_clst_and_umap = outdir + '/intg_summary_{}.tsv'.format(name)
output_cluster_centroids = outdir + '/centroids_{}.pkl'.format(name)
output_figures = outdir + '/{}_{{}}.{{}}'.format(name)

# Get dataset configuration
settings = collections.OrderedDict()

# Load in dataset metadata
dataset_metadata_filename = os.path.join(DATA_DIR, "dataset_metadata.csv")
assert(os.path.isfile(dataset_metadata_filename))    
metadata = pd.read_csv(dataset_metadata_filename)

# Load metadata file into settings OrderedDict
Mod_info = collections.namedtuple('Mod_info', metadata.columns)
for _, row in metadata.iterrows():
    sample_name = row['mod']
    settings[sample_name] = Mod_info(*row)

# Select modalities and features
mods_selected = args.mods
features_selected = args.features

for features_modality in features_selected:
    assert (features_modality in mods_selected)

data_f = os.path.join(DATA_DIR, "{0}.h5ad")

# Within modality
assert(len(args.ps) == 3)
assert(all(isinstance(ps, float) for ps in args.ps))
assert(len(args.n_dropped_pcs) == 3)
assert(all(isinstance(n_dropped_pcs, int) for n_dropped_pcs in args.n_dropped_pcs))

ps = {'mc': args.ps[0],
      'rna': args.ps[1],
      'atac': args.ps[2],
     }
drop_npcs = {
      'mc': args.n_dropped_pcs[0],
      'rna': args.n_dropped_pcs[1],
      'atac': args.n_dropped_pcs[2],
     }

# across modality
cross_mod_distance_measure = args.cross_mod_distance
knn = args.nearest_neighbors
relaxation = args.relaxation
n_cca = args.n_cca

# PCA
npc = args.num_pcs

# clustering
k = args.k_clusters
resolutions = args.resolutions
# umap
umap_neighbors = args.umap_neighbors
min_dist = args.min_distance

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
