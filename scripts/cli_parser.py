"""Command line interface is defined here.
"""

DESCRIPTION="""
SingleCellFusion is a computational tool to integrate single-cell transcriptome and epigenome datasets. 
"""

DESCRIPTION_preproc="""
SingleCellFusion is a computational tool to integrate single-cell transcriptome and epigenome datasets. 
This is the CLI for its preprocessing module 
(from count matrices to normalized HVG feature matrices).
"""

EPILOG="""
Contributors: Fangming Xie, Aditya Chandrasekar, Wayne I. Doyle, Ethan J. Armand, Eran Mukamel
Contact: Eran Mukamel (emukamel@ucsd.edu)
"""
import argparse
import os

def create_parser_preproc():
    """
    """
    parser = argparse.ArgumentParser(
        prog="SingleCellFusion_pre",
        description=DESCRIPTION_preproc,
        epilog=EPILOG,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    
    # Input/Output Dataset Settings
    parser.add_argument(
        "-i", "--input_datasets", 
        type=str,
        nargs="+",
        required=True,
        help='''(list of str; required) 
             Paths to .h5ad files, each containing a cell-by-gene feature matrix, 
             cell IDs and gene IDs. Cell IDs should be unique within each .h5ad file, 
             Gene IDs should be shared or partially shared across files. 
             Multiple inputs should be listed as a space seperated list of filenames. 
             '''
    )
    parser.add_argument(
        "-inorm", "--input_normalizations", 
        type=str,
        nargs="+",
        required=True,
        help='''(list of str; required) 
             Data modalities chosen from 'mc', 'cpm', or 'tpm'. This should be 
             listed in the same order as input_datasets 
             ''',
    )
    parser.add_argument(
        "-o", "--output_dir", 
        type=str,
        default="./results",
        help='''(str) 
             Directory to store output files
             '''
    )
    parser.add_argument(
        "-op", "--output_prefix", 
        metavar="OUT_PREFIX", 
        type=str,
        default="test_run",
        help='''(str) 
             The output files will contain this prefix
             '''
    )
    parser.add_argument(
        "-ga", "--gene_annotation_file", 
        type=str,
        default="",
        help='''(str) 
             Gene annotation file (bed format: chr, start, end, gene_id, gene_name)
             '''
    )
    return parser

def create_parser():
    """
    """
    parser = argparse.ArgumentParser(
        prog="SingleCellFusion",
        description=DESCRIPTION,
        epilog=EPILOG,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    
    ## ARGUMENTS DIRECTLY FED INTO SingleCellFusion CLI 
    # Input/Output Dataset Settings
    parser.add_argument(
        "-i", "--input_datasets", 
        metavar="xx.h5ad",
        type=str,
        nargs="+",
        required=True,
        help='''(list of str; required) 
             Paths to .h5ad files, each containing a cell-by-gene feature matrix, 
             cell IDs and gene IDs. Cell IDs should be unique within each .h5ad file, 
             Gene IDs should be shared or partially shared across files. 
             Multiple inputs should be listed as a space seperated list of filenames.
             '''
    )
    parser.add_argument(
        "-im", "--input_modalities", 
        metavar="rna/atac/mc",
        type=str,
        nargs="+",
        required=True,
        help='''(list of str; required)
             Data modalities chosen from 'rna', 'atac', or 'mc'. This should be 
             listed in the same order as input_datasets. 
             '''
    )
    # may need this in the future
    # parser.add_argument(
    #     "-im", "--input_meta", 
    #     type=str,
    #     required=True,
    #     help="(list of str) Input metadata csv file",
    # )

    parser.add_argument(
        "-f", "--feature_datasets", 
        metavar="xx.h5ad",
        type=str,
        nargs="+",
        required=True,
        help='''(list of str; required)
             Dataset(s) whose features all other datasets will impute into.
             Enter multiple datasets as a space-separated list of filenames.
             The features of these datasets will
             be the features kept in the output imputed data table.",
             '''
    )
    parser.add_argument(
        "-o", "--output_dir", 
        metavar="DIR",
        type=str,
        default="./results",
        help='''(str)
             Directory to store output files
             '''
    )
    parser.add_argument(
        "-op", "--output_prefix", 
        type=str,
        default="SingleCellFusion",
        help='''(str)
             The output files will contain this prefix.
             '''
    )

    # within modality smoothing 
    parser.add_argument(
        "--num_pcs", 
        type=int,
        default=50,
        help='''(integer)
             Number of Principal Components to keep for each dataset 
             for smoothing and for clustering/embedding after imputation.
             '''
    )
    parser.add_argument(
        "--smoothing_fractions", 
        default=[0.7, 0.1, 0.9],
        help='''(list of floats) 
             A list of three values between 0 to 1 that controls the relative contribution
             from the cell itself vs. its neighbors in within-dataset smoothing, 
             specified for 'rna', 'atac', 'mc' data, respectively.
             '''
    )

    # constraint kNN across modalities
    parser.add_argument(
        "--nearest_neighbors", 
        type=int,
        default=20,
        help='''(integer)
             Number of nearest neighbors used to impute data
             '''
    )
    parser.add_argument(
        "--relaxation", 
        type=int,
        default=3,
        help='''(integer)
             A value between 1 to infinity. It specifies how much the datasets should be
             assumed to be structurally similar.
             relaxation=1 enforces a hard limit that every cell receives equal number of nearest neighbors
             relaaxation=infinity approaches traditional kNN
             '''
    )

    # Arguments for Clustering
    parser.add_argument(
        "--leiden_n_neighbors", 
        type=int,
        default=30,
        help='''(integer) 
             Number of nearest neighbors to form in the integrated space, 
             the resulting nearest neighbor graph is used for Leiden clustering.
             It is passed into the python package leidenalg.
             '''
    )
    parser.add_argument(
        "--leiden_resolutions", 
        type=list,
        default=[0.1, 0.2, 0.4, 0.8],
        help='''(list of floats) 
             A list of resolutions to be used for Leiden Clustering.
             It is passed into the python package leidenalg.
             '''
    )

    # Arguments for UMAP
    parser.add_argument(
        "--umap_n_neighbors", 
        type=int,
		default=60,
        help='''(integer)
             Number of neighbors for UMAP. It is passed into the python package umap.UMAP(n_neighbors).
             '''
    )
    parser.add_argument(
        "--umap_min_dist", 
        type=float,
        default=0.5,
        help='''(float)
             Minimum distance for UMAP. It is passed into the python package umap.UMAP(min_dist).
             '''
    )
    return parser

def parse_filename(data_file):
    """turn a xxx/xxx/XXXX.h5ad into XXXX 
    """
    dataset_name = os.path.basename(data_file)
    if dataset_name.endswith('.h5ad'):
        dataset_name = dataset_name[:-len('.h5ad')]
    else:
        raise ValueError("filenames don't have the format xxxx.h5ad")
    return dataset_name

def modality_default_options(mod):
    """
    """
    if mod == 'mc':
        mod_direction = -1
        # norm_option = 'mc'
    elif mod == 'rna':
        mod_direction = 1
        # norm_option = 'cpm'
    elif mod == 'atac':
        mod_direction = 1
        # norm_option = 'tpm'
    else:
        raise ValueError("choose from ['mc', 'rna', 'atac']")
    return mod_direction 
