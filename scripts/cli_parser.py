"""Command line interface is defined here.
"""

DESCRIPTION="""
SCF is a computational tool to integrate single-cell transcriptome and epigenome datasets. 

Contributors: Fangming Xie, Aditya Chandrasekar, Wayne I. Doyle, Ethan J. Armand, Eran Mukamel
Contact: Eran Mukamel (emukamel@ucsd.edu)
"""
import argparse

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
        default="test_run",
        help="Name of current configuration (prefix used in the output)",
    )

    # Input/Output Dataset Settings
    # what's required false/true?
    parser.add_argument(
        "-d", "--dataset", 
        metavar="DATASET_DIR", 
        type=str,
        default="./datasets",
        required=False,
        help="Directory containing input gene x cell tables in h5ad format",
    )

    parser.add_argument(
        "-m", "--files", 
        metavar="MODALITY", 
        type=str,
        nargs="+",
        required=True,
        help="Modalities inside the data directory to be integrated, \
            entered as a space seperated list of filenames. \
                '10x_cells_v2' \
                ",
    )

    # name of the cell col
    # name of the cluster col [option]
    # [option]
    # [option]

    parser.add_argument(
        "-m", "--", 
        metavar="MODALITY", 
        type=str,
        nargs="+",
        required=True,
        help="Modalities inside the data directory to be integrated, \
            entered as a space seperated list of filenames. \
                '10x_cells_v2' \
                ",
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

    # 
    parser.add_argument(
        "--smoothing_parameter", 
        metavar="PS",
        default=[0.9, 0.7, 0.1],
        help="a value from 0 ~ 1, controlling the relative contribution of the cell itself vs. its neighbors \
              TODO: Specified list of PS for [mc, rna, atac] \
              data respectively",
    )


    # enumerate options
    # make types explicit
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

    # parser.add_argument(
    #     "--n_cca", 
    #     metavar="NUMBER_CCA", 
    #     type=int,
    #     default=0,
    #     help="TODO",
    # )

    # # Number of PCs
    # parser.add_argument(
    #     "--n_dropped_pcs", 
    #     metavar="NUM_DROPPED_PCS", 
    #     type=list,
    #     default=[0, 0, 0],
    #     help="Specified list of Principal Components to drop for \
    #             [mc, rna, atac] respectively",
    # )
   
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