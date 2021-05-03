#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
from scipy import sparse
import time
import re
import anndata
import logging
import os

import basic_utils
import preproc_utils
import cli_parser

def get_gene_annotation(gene_annotation_file):
    """
    """
    genes = pd.read_csv(
        gene_annotation_file,
        sep='\t',
        header=None,
        usecols=[0,1,2,3,4],
        ).rename(columns={0: 'chr', 
                            1: 'start',
                            2: 'end',
                            3: 'ensid',
                            4: 'gene',
                        })
    return genes

def preproc_sparse(
    f_data, 
    f_hvftr_data, 
    normalization_option, 
    gene_lengths_base='', # required if normalization option == "tpm"
    gid_col='', 
    cid_col='',
    ):
    """Generate normalized HVG matrices from raw count matrices
    """
    # # highly variable features
    ti = time.time()
    logging.info("Preprocessing {}".format(f_data))

    # read data matrix
    if normalization_option == 'mc':
        f_data = f_data_format.format(SRC_DIR, mod)
        
        # read in files
        print(mod, "Reading in files {}".format(time.time()-ti))
        gxc_raw = snmcseq_utils.load_gc_matrix_methylation(f_data_gene, f_data_cell, f_data_mc, f_data_c)
        print(gxc_raw.data['mc'].shape, gxc_raw.data['c'].shape)
        print(time.time()-ti)
        
        # output file
        f_hvftr_data_methylation = f_hvftr_format.format(DST_DIR, mod, 'tsv') 
        print(time.time()-ti)
        
        # check meta cells agree with gxc cells
        assert np.all(meta.index.values == gxc_raw.cell)
        # check genes are uniq 
        assert len(gxc_raw.gene) == len(np.unique(gxc_raw.gene)) 
        # do
        gxc_hvftr = preproc_utils.preproc_methylation(
                                                    gxc_raw,
                                                    meta,
                                                    global_value_col=settings[mod].global_mean, 
                                                    base_call_cutoff=20, 
                                                    sufficient_coverage_fraction=0.95,
                                                    hv_percentile=30,
                                                    n_qcut=10,
                                                    )
        # save
        print(mod, "Saving to files {}".format(time.time()-ti))
#         gxc_hvftr.to_csv(f_hvftr_data_methylation, sep="\t", header=True, index=True, na_rep='NA')
        h5ad_mat_hvftr.write(f_hvftr_data, compression='gzip')
        
    else:
        # read in files
        logging.info("Reading in file {}".format(f_data))
        h5ad_mat = anndata.read_h5ad(f_data)
        meta, gxc_raw = basic_utils.h5ad_to_scf_rna_format(h5ad_mat, gid_col, cid_col)
        
        # check meta cells agree with gxc cells
        assert np.all(meta.index.values == gxc_raw.cell)
        # check genes are uniq 
        assert len(gxc_raw.gene) == len(np.unique(gxc_raw.gene)) 
    
        # get hvftrs
        logging.info("Preproc and get highly variable genes {}".format(f_data))
        if normalization_option == 'cpm':
            gxc_hvftr = preproc_utils.preproc_rna_cpm_based(
                                            gxc_raw, 
                                            sufficient_cell_coverage=0.01, 
                                            hv_percentile=30, hv_ncut=10)
        elif normalization_option == 'tpm':
            gene_lengths = gene_lengths_base.reindex(gxc_raw.gene)
            gxc_hvftr = preproc_utils.preproc_rna_tpm_based(
                                            gxc_raw, gene_lengths, impute_gene_lengths=True, 
                                            sufficient_cell_coverage=0.01, 
                                            hv_percentile=30, hv_ncut=10)
    
        # save
        logging.info("Saving to file {}".format(f_hvftr_data))
        h5ad_mat_hvftr = basic_utils.scf_rna_format_to_h5ad(meta, gxc_hvftr)
        h5ad_mat_hvftr.write(f_hvftr_data, compression='gzip')
    
    print("Done: {} -> {}".format(f_data, f_hvftr_data))
    return 

def preproc_dense(
    f_data, 
    f_hvftr_data, 
    normalization_option, 
    gene_lengths_base='', # required if normalization option == "tpm"
    gid_col='', 
    cid_col='',
    ):
    """Generate normalized HVG matrices from raw count matrices
    """
    # # highly variable features
    ti = time.time()
    logging.info("Preprocessing {}".format(f_data))

    # read data matrix
    if normalization_option == 'mc':
        f_data = f_data_format.format(SRC_DIR, mod)
        
        # read in files
        print(mod, "Reading in files {}".format(time.time()-ti))
        gxc_raw = snmcseq_utils.load_gc_matrix_methylation(f_data_gene, f_data_cell, f_data_mc, f_data_c)
        print(gxc_raw.data['mc'].shape, gxc_raw.data['c'].shape)
        print(time.time()-ti)
        
        # output file
        f_hvftr_data_methylation = f_hvftr_format.format(DST_DIR, mod, 'tsv') 
        print(time.time()-ti)
        
        # check meta cells agree with gxc cells
        assert np.all(meta.index.values == gxc_raw.cell)
        # check genes are uniq 
        assert len(gxc_raw.gene) == len(np.unique(gxc_raw.gene)) 
        # do
        gxc_hvftr = preproc_utils.preproc_methylation(
                                                    gxc_raw,
                                                    meta,
                                                    global_value_col=settings[mod].global_mean, 
                                                    base_call_cutoff=20, 
                                                    sufficient_coverage_fraction=0.95,
                                                    hv_percentile=30,
                                                    n_qcut=10,
                                                    )
        # save
        print(mod, "Saving to files {}".format(time.time()-ti))
#         gxc_hvftr.to_csv(f_hvftr_data_methylation, sep="\t", header=True, index=True, na_rep='NA')
        h5ad_mat_hvftr.write(f_hvftr_data, compression='gzip')
        
    else:
        # read in files
        logging.info("Reading in file {}".format(f_data))
        h5ad_mat = anndata.read_h5ad(f_data)
        meta, gxc_raw = basic_utils.h5ad_to_scf_rna_format(h5ad_mat, gid_col, cid_col)
        
        # check meta cells agree with gxc cells
        assert np.all(meta.index.values == gxc_raw.cell)
        # check genes are uniq 
        assert len(gxc_raw.gene) == len(np.unique(gxc_raw.gene)) 
    
        # get hvftrs
        logging.info("Preproc and get highly variable genes {}".format(f_data))
        if normalization_option == 'cpm':
            gxc_hvftr = preproc_utils.preproc_rna_cpm_based(
                                            gxc_raw, 
                                            sufficient_cell_coverage=0.01, 
                                            hv_percentile=30, hv_ncut=10)
        elif normalization_option == 'tpm':
            gene_lengths = gene_lengths_base.reindex(gxc_raw.gene)
            gxc_hvftr = preproc_utils.preproc_rna_tpm_based(
                                            gxc_raw, gene_lengths, impute_gene_lengths=True, 
                                            sufficient_cell_coverage=0.01, 
                                            hv_percentile=30, hv_ncut=10)
    
        # save
        logging.info("Saving to file {}".format(f_hvftr_data))
        h5ad_mat_hvftr = basic_utils.scf_rna_format_to_h5ad(meta, gxc_hvftr)
        h5ad_mat_hvftr.write(f_hvftr_data, compression='gzip')
    
    print("Done: {} -> {}".format(f_data, f_hvftr_data))
    return 


        
if __name__ == "__main__":
    log = basic_utils.create_logger()

    parser = cli_parser.create_parser_preproc()
    args = parser.parse_args()
    logging.info('* Parsing Command Line Arguments')

    # get input files
    data_files = args.input_datasets
    mods_selected = [cli_parser.parse_filename(data_file) for data_file in data_files]

    # specify output files
    outdir = args.output_dir
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    outprefix = args.output_prefix

    output_files = [
        os.path.join(outdir, "{}_{}".format(outprefix, os.path.basename(input_file)))
        for input_file in data_files
    ]

    # get dataset metadata
    input_normalizations = args.input_normalizations
    for option in input_normalizations:
        assert (option in ['mc', 'cpm', 'tpm'])

    # parameters
    gene_annotation_file = args.gene_annotation_file

    # # Settings
    df_genes = get_gene_annotation(gene_annotation_file).set_index('ensid')
    gene_lengths_base = (df_genes['end'] - df_genes['start'])

    for data_file, output_file, norm_option in zip(
        data_files, output_files, input_normalizations):

        preproc_dense(
            data_file, 
            output_file, 
            norm_option, 
            gene_lengths_base=gene_lengths_base, # required if normalization option == "tpm"
            gid_col='ensid', 
            cid_col='',
        )
        break

