"""Set up dataset-specific settings
"""

import collections

# meta settings
mods = (
    'snmcseq_gene',
    'smarter_cells',
    'smarter_nuclei',
    '10x_cells_v2', 
    'snatac',
)

Mod_info = collections.namedtuple('Mod_info', [
    'mod', # name of the data modality (used internally as prefix in every related file name) 
    'name', # name of the data modality (used to make plots)
    'mod_category', # choose from rna/atac/mc
    'norm_option', # choose from mc/cpm/tpm
    'mod_direction', # +1 or -1; +1 for RNA/ATAC/MERFISH; -1 for DNA methylation
    'cell_col', # cell id
    'cluster_col', # cluster label 
    'annot_col', # cluster annotation 
    'global_mean', # for DNA methylation, this is global mean mCH level. Leave blank for other modalities
    'global_mean_mcg', # for mcg only; leave it blank if it's not DNA methylation data
    'color',
    'species',
],
defaults=('ID','NAME','rna','tpm',1,'','','','','','black','mouse')
)

# settngs
settings_mc = Mod_info(
    mods[0],
    'DNA methylation',
    'mc',
    'mc', 
    -1, # negative direction 
    'cell',
    'SubCluster',
    'SubCluster', #'major_clusters' 'sub_cluster' 
    'CH_Rate',
    'CG_Rate',
    '#6DC8BF',
    'mouse',
)


settings_sc = Mod_info(
    mods[1],
    'SmartSeq_cells_AIBS',
    'rna',
    'tpm', 
    +1, # direction 
    'cell',
    'cluster_id',
    'cluster_label',
    '',
    '',
    '#FF3333', 
    'mouse',
)

settings_sn = Mod_info(
    mods[2],
    'SmartSeq_nuclei_AIBS',
    'rna',
    'tpm', 
    +1, # direction 
    'cell',
    'cluster_id',
    'cluster_label',
    '',
    '',
    '#7D9B3D',
    'mouse',
)

settings_10xc = Mod_info(
    mods[3],
    '10X_cells_v2_AIBS',
    'rna',
    'cpm', 
    +1, # direction 
    'cell',
    'cluster_id',
    'cluster_label',
    '',
    '',
    '#F68B1F',
    'mouse',
)

settings_snatac = Mod_info(
    mods[4],
    'snATAC-seq',
    'atac',
    'tpm', 
    +1, # direction 
    'cell',
    '',
    '',
    '',
    '',
    'C1',
    'mouse',
)


settings = collections.OrderedDict({
    mods[0]: settings_mc,
    mods[1]: settings_sc,
    mods[2]: settings_sn,
    mods[3]: settings_10xc,
    mods[4]: settings_snatac,
})

