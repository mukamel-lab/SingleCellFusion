#!/usr/bin/env python3
"""An example configuration file
"""
import sys
import os

# EAM: Figure out the data path on Google workspace:
# TEMPORARY HACK
sample = 'example_l5pt'

import pathlib 
currdir=str(pathlib.Path().absolute())
DATA_DIR = currdir+'/SingleCellFusion/'+sample+'/datasets'

# # Configs  
name = 'test_scf_mouse_mop'
outdir = currdir+'/SingleCellFusion/'+sample+'/results'
output_pcX_all = outdir + '/pcX_all_{}.npy'.format(name)
output_cells_all = outdir + '/cells_all_{}.npy'.format(name)
output_imputed_data_format = outdir + '/imputed_data_{}_{{}}.npy'.format(name)
output_clst_and_umap = outdir + '/intg_summary_{}.tsv'.format(name)
output_cluster_centroids = outdir + '/centroids_{}.pkl'.format(name)
output_figures = outdir + '/{}_{{}}.{{}}'.format(name)

# fixed dataset configs
sys.path.insert(0, DATA_DIR)
from __init__datasets import *

meta_f = os.path.join(DATA_DIR, '{0}_metadata.tsv')
hvftrs_f = os.path.join(DATA_DIR, '{0}_hvfeatures.{1}')
hvftrs_gene = os.path.join(DATA_DIR, '{0}_hvfeatures.gene')
hvftrs_cell = os.path.join(DATA_DIR, '{0}_hvfeatures.cell')

mods_selected = [
    'snmcseq_gene',
    'smarter_cells',
    'smarter_nuclei',
    '10x_cells_v2', 
    ]
features_selected = ['10x_cells_v2']
# check features
for features_modality in features_selected:
    assert (features_modality in mods_selected)

# within modality
ps = {'mc': 0.9,
      'rna': 0.7,
     }
drop_npcs = {
      'mc': 0,
      'rna': 0,
     }

# across modality
cross_mod_distance_measure = 'correlation' # cca
knn = 20 
relaxation = 3
n_cca = 0

# PCA
npc = 50

# clustering
k = 30
resolutions = [0.1, 0.2, 0.4, 0.8]
# umap
umap_neighbors = 60
min_dist = 0.5
