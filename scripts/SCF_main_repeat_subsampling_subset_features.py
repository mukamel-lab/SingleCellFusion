#!/usr/bin/env python3
"""SingleCellFusion main rontine"""

from __init__ import *

from scipy import sparse
import collections
import itertools
import sys
import pickle
import argparse

import basic_utils
import SCF_utils
import cluster_cv_utils


def create_parser():
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config_py", help="Configuration file", required=True)
    parser.add_argument("-s", "--subsample_fraction", help="0~1", type=float, required=True)
    parser.add_argument("-sn", "--subsample_times", help=">1", type=int, required=True)
    parser.add_argument("-fs", "--feature_subsample_fraction", help="0~1", type=float, required=True)
    parser.add_argument("-fsn", "--feature_subsample_times", help=">1", type=int, required=True)
    parser.add_argument("-ga", "--gene_annotation", help="gene annotation file (.bed)", type=str, required=True)
    return parser

def subsampling(mods_selected, settings, metas, gxc_hvftrs, p=0, n=0):
    """Do many datasets at the same time
    p - fraction of cells from each dataset to be included
    """
    if p < 1:
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
    else:
        return metas, gxc_hvftrs



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

# get gene chrom lookup table (Dec 29, 2020)
gene_chrom_lookup = pd.read_csv(args.gene_annotation, sep='\t', header=None)[[0, 3]]
gene_chrom_lookup['gid'] = gene_chrom_lookup[3].apply(lambda x: x.split('.')[0])
gene_chrom_lookup = gene_chrom_lookup.set_index('gid')[0]
logging.info("gene chrom lookup table: {}".format(gene_chrom_lookup.shape))

# output chroms used (Dec 29, 2020)
output_chroms = outdir + "/chroms_{}.txt".format(name)

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

feature_subsample_fraction = args.feature_subsample_fraction
feature_subsample_times = args.feature_subsample_times

# subsampling 80% of cells from each modality
for i in range(subsample_times):
    output_tag = ".{}".format(i)

    logging.info("Subsampling {}/{}".format(i+1, subsample_times))
    # subset cells
    metas_sub, gxc_hvftrs_sub = subsampling(
                mods_selected, settings, metas, gxc_hvftrs, p=subsample_fraction,
                )
    for mod in mods_selected:
        fout = (output_cells+output_tag).format(mod)
        cells_mod = metas_sub[mod].index.values
        np.save(fout, cells_mod)

    # subset features
    for n_split in np.arange(feature_subsample_times):
        # split features
        gene_set_lookup, chroms_0, chroms_1 = cluster_cv_utils.split_genes(
                                                        gene_chrom_lookup, 
                                                        split_frac=feature_subsample_fraction, 
                                                        report_chroms=True,
                                                        )
        gxc_hvftrs_sub_g0, gxc_hvftrs_sub_g1 = cluster_cv_utils.split_features_routine(mods_selected, settings, gxc_hvftrs_sub, gene_set_lookup)
        # use g0; record using g0 
        output_tag2 = ".f{}".format(n_split)
        fout = (output_chroms+output_tag+output_tag2)
        basic_utils.export_single_textcol(fout, chroms_0)
        logging.info('Saved {}'.format(fout))

        # ## run SCF
        pcX_all, cells_all = SCF_utils.core_scf_routine(mods_selected, features_selected, settings, 
                                                        metas_sub, gxc_hvftrs_sub_g0, # sub 
                                                        ps, drop_npcs,
                                                        cross_mod_distance_measure, knn, relaxation, n_cca,
                                                        npc,
                                                        output_pcX_all+output_tag+output_tag2, 
                                                        output_cells_all+output_tag+output_tag2,
                                                        output_imputed_data_format+output_tag+output_tag2,
                                                        ka_smooth=ka_smooth,
                                                        save_knn=save_knn,
                                                        output_knn_within=output_knn_within+output_tag+output_tag2,
                                                        output_knn_across=output_knn_across+output_tag+output_tag2,
                                                        )
        logging.info('Done integration into a common PC space')
