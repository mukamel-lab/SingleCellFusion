"""Command line interface is defined here.
"""

DESCRIPTION="""
SCF is a computational tool to integrate single-cell transcriptome and epigenome datasets. 
"""

EPILOG="""
Contributors: Fangming Xie, Aditya Chandrasekar, Wayne I. Doyle, Ethan J. Armand, Eran Mukamel
Contact: Eran Mukamel (emukamel@ucsd.edu)
"""
import argparse
import os

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
        "-i", "--inputs", 
        metavar="INPUT_FILES", 
        type=str,
        nargs="+",
        required=True,
        help="(list of str) Datasets inside the data directory to be integrated, \
              entered as a space seperated list of filenames. \
              './10x_cells_v2.h5ad' \
             ",
    )
    parser.add_argument(
        "-ic", "--input_file_categories", 
        type=str,
        nargs="+",
        required=True,
        # choices=['mc', 'rna', 'atac'],
        help="(list of str) Choose from ['mc', 'rna', 'atac']",
    )
    # may need this in the future
    # parser.add_argument(
    #     "-im", "--input_meta", 
    #     type=str,
    #     required=True,
    #     help="(list of str) Input metadata csv file",
    # )

    parser.add_argument(
        "-f", "--features", 
        metavar="FEATURE", 
        type=str,
        nargs="+",
        required=True,
        help="(list of str) Datasets to impute into, entered as a space-separated \
              list of filenames. The features of these datasets will \
              be the features kept in the output imputed data table.",
    )
    parser.add_argument(
        "-o", "--output", 
        metavar="OUT_DIR", 
        type=str,
        default="./results",
        help="(str) Directory to store output files",
    )
    parser.add_argument(
        "-op", "--output_prefix", 
        metavar="OUT_PREFIX", 
        type=str,
        default="test_run",
        help="(str) The output files will contain this prefix",
    )

    # within modality smoothing 
    parser.add_argument(
        "--num_pcs", 
        metavar="PCs", 
        type=int,
        default=50,
        help="(integer) Number of Principal Components to keep for each dataset \
              for smoothing and for clustering/embedding after imputation",
    )
    parser.add_argument(
        "--smoothing_fractions", 
        metavar="SMOOTHING_FRACTIONS",
        default=[0.7, 0.1, 0.9],
        help="(list of floats) values of 0 ~ 1, controlling the relative contribution from the cell itself vs. its neighbors. \
              Specified list of PS for [rna, atac, mc] data respectively.",
    )

    # constraint kNN across modalities
    parser.add_argument(
        "--nearest_neighbors", 
        metavar="kNN", 
        type=int,
        default=20,
        help="(integer) Number of nearest neighbors used to impute data",
    )
    parser.add_argument(
        "--relaxation", 
        type=int,
        default=3,
        help="(integer) Between 1 to infinity. It specifies how much the datasets should be \
              assumed to be structurally similar. \
              relaxation=1 enforces a hard limit that every cell receives equal number of nearest neighbors\
              relaaxation=infinity approaches traditional kNN",
    )

    # Arguments for Clustering
    parser.add_argument(
        "--leiden_n_neighbors", 
        type=int,
        default=30,
        help="(integer) Number of nearest neighbors to form in the integrated space, \
              the resulting kNN graph is used for Leiden clustering.",
    )
    parser.add_argument(
        "--leiden_resolutions", 
        type=list,
        default=[0.1, 0.2, 0.4, 0.8],
        help="(list of floats) A list of resolutions to be used for Leiden Clustering.",
    )

    # Arguments for UMAP
    parser.add_argument(
        "--umap_n_neighbors", 
        type=int,
		default=60,
        help="(integer) Number of neighbors for UMAP. It is passed into UMAP n_neighbors.",
    )
    parser.add_argument(
        "--umap_min_dist", 
        type=float,
        default=0.5,
        help="(float) Minimum distance for UMAP. It is passed into UMAP min_dist.",
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