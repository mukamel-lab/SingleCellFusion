#!/usr/bin/env python3
"""SingleCellFusion main routine"""

from __init__ import *

from scipy import sparse
import collections
import itertools
import sys
import os
import pickle
import argparse
import anndata

import pandas # To import dataset_metadata.csv

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

DESCRIPTION="""
SCF is a computational tool to integrate single-cell transcriptome and epigenome datasets. 

Contributors: Fangming Xie, Aditya Chandrasekar, Wayne I. Doyle, Ethan J. Armand, Eran Mukamel
Contact: Eran Mukamel (emukamel@ucsd.edu)
"""

def create_parser():
    """
    """
    parser = argparse.ArgumentParser(prog="SingleCellFusion",
                                     description=DESCRIPTION,
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    ## ARGUMENTS DIRECTLY FED INTO SingleCellFusion CLI 

    # Name of settings
    # where does it show up in the results?
    parser.add_argument(
        "-n", "--name", 
        metavar="NAME", 
        type=str,
        default="test_scf_biccn2",
        help="Name of current configuration",
    )

    # Input/Output Dataset Settings
    # what's required false/true?
    parser.add_argument(
        "-d", "--dataset", 
        metavar="DATASET_DIR", 
        type=str,
        default="./datasets",
        required=False,
        help="Directory containing input gene x cell tables as h5ad",
    )

    parser.add_argument(
        "-m", "--mods", 
        metavar="MODALITY", 
        type=str,
        nargs="+",
        required=True,
        help="Modalities inside the data directory to be integrated, \
            entered as a space seperated list of filenames.",
    )
    
    parser.add_argument(
        "-f", "--features", 
        metavar="FEATURE", 
        type=str,
        nargs="+",
        required=True,
        help="Modalities to impute into, entered as a space-separated \
            list of filenames. The features of these modalities will \
            be the features kept in the output imputed data table.",
    )
    
    parser.add_argument(
        "-o", "--output", 
        metavar="OUT_DIR", 
        type=str,
        default="./results",
        help="Directory to store output of SingleCellFusion",
    )

    # Measures Within Each Modality
    parser.add_argument(
        "--num_pcs", 
        metavar="PCs", 
        type=int,
        default=50,
        help="Number of Principal Components to use for main function",
    )

    parser.add_argument(
        "--ps", 
        metavar="PS",
        default=[0.9, 0.7, 0.1],
        help="TODO: Specified list of PS for [mc, rna, atac] \
              data respectively",
    )


    parser.add_argument(
        "--cross_mod_distance", 
        metavar="DISTANCE_METRIC", 
        type=str,
        default="correlation",
        help="Distance for comparisons across modalities",
    ) #correlation or cca

    parser.add_argument(
        "--nearest_neighbors", 
        metavar="kNN", 
        type=int,
        default=20,
        help="Number of nearest neighbors used to impute data",
    )

    parser.add_argument(
        "--relaxation", 
        metavar="RELAXATION", 
        type=int,
        default=3,
        help="A integer between 1 to infinity. It specifies how much the datasets should be \
            assumed to be structurally similar. \
            relaxation=1 enforces a hard limit that every cell receives equal number of nearest neighbors\
            relaaxation=infinity approaches traditional kNN",
    )

    parser.add_argument(
        "--n_cca", 
        metavar="NUMBER_CCA", 
        type=int,
        default=0,
        help="TODO",
    )

    # Number of PCs
    parser.add_argument(
        "--n_dropped_pcs", 
        metavar="NUM_DROPPED_PCS", 
        type=list,
        default=[0, 0, 0],
        help="Specified list of Principal Components to drop for \
                [mc, rna, atac] respectively",
    )
   
    # Arguments for Clustering
    parser.add_argument(
        "-k", "--k_clusters", 
        metavar="K_CLUSTERS", 
        type=int,
        default=30,
        help="Number of Clusters for late clustering stage of SingleCellFusion",
    )
   
    parser.add_argument(
        "--resolutions", 
        metavar="RESOLUTIONS", 
        type=list,
        default=[0.1, 0.2, 0.4, 0.8],
        help="Series of Resolutions to be used for Leiden Clustering Algorithm. \
                Should be given as list of resolutions, to be iterated through and used \
                from left to right.",
        )

    # Arguments for UMAP
    # TODO: get documentation for --n_neighbors & --min_dist respectively
    parser.add_argument(
        "-u", "--umap_neighbors", 
        metavar="UMAP_NEIGHBORS", 
        type=int,
		default=60,
        help="TODO: Number of Neighbors used for UMAP stage of SingleCellFusion",
    )
   
    parser.add_argument(
        "--min_distance", 
        metavar="MIN_UMAP_DISTANCE", 
        type=float,
        default=0.5,
        help="TODO: Minimum Distance bused be used for UMAP stage of SingleCellFusion",
        )

    return parser

parser = create_parser()
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
