# SingleCellFusion

SingleCellFusion is a computational tool to integrate single-cell transcriptome and epigenome datasets. Code in this repository is used in [Luo et al., (2019) *BioRxiv*](https://www.biorxiv.org/content/10.1101/2019.12.11.873398v1) and in [Yao et al., (2020) *BioRxiv*](https://www.biorxiv.org/content/10.1101/2020.02.29.970558v2).

For more information and to cite this work:
- Luo, C. et al. Single nucleus multi-omics links human cortical cell regulatory genome diversity to disease risk variants. bioRxiv 2019.12.11.873398 (2019) [doi:10.1101/2019.12.11.873398](https://www.biorxiv.org/content/10.1101/2019.12.11.873398v1)
- Yao, Z. et al. An integrated transcriptomic and epigenomic atlas of mouse primary motor cortex cell types. bioRxiv 2020.02.29.970558 (2020) [doi:10.1101/2020.02.29.970558](https://www.biorxiv.org/content/10.1101/2020.02.29.970558v2)
- BRAIN Initiative Cell Census Network (BICCN) et al. A multimodal cell census and atlas of the mammalian primary motor cortex. bioRxiv 2020.10.19.343129 (2020) [doi:10.1101/2020.10.19.343129](https://www.biorxiv.org/content/10.1101/2020.10.19.343129v1)
- A [brief description](docs/scf_description.rst) of how SingleCellFusion works

Code contributors: [Fangming Xie](mailto:f7xie@ucsd.edu), Aditya Chandrasekar,  Wayne I. Doyle, [Ethan Armand](mailto:ejarmand@ucsd.edu)

Contact: [Eran Mukamel](mailto:emukamel@ucsd.edu)

## Installation
Step 1: Clone this repo.
```bash
git clone https://github.com/mukamel-lab/SingleCellFusion.git
cd SingleCellFusion
```

**TBD**
Step 2: Set up a conda environment and install dependent packages. (Skip this step if not needed.)
```bash
conda env create -f environment.yml # create an env named scf_dev
source activate scf_dev
```

## Usage
```
usage: SingleCellFusion [-h] -i xx.h5ad [xx.h5ad ...] -im rna/atac/mc [rna/atac/mc ...] -f xx.h5ad [xx.h5ad ...] [-o DIR] [-op OUTPUT_PREFIX]
                        [--num_pcs NUM_PCS] [--smoothing_fractions SMOOTHING_FRACTIONS] [--nearest_neighbors NEAREST_NEIGHBORS]
                        [--relaxation RELAXATION] [--leiden_n_neighbors LEIDEN_N_NEIGHBORS] [--leiden_resolutions LEIDEN_RESOLUTIONS]
                        [--umap_n_neighbors UMAP_N_NEIGHBORS] [--umap_min_dist UMAP_MIN_DIST]

SingleCellFusion is a computational tool to integrate single-cell transcriptome and epigenome datasets.

optional arguments:
  -h, --help            show this help message and exit
  -i xx.h5ad [xx.h5ad ...], --input_datasets xx.h5ad [xx.h5ad ...]
                        (list of str; required) Paths to .h5ad files, each containing a cell-by-gene feature matrix, cell IDs and gene IDs. Cell
                        IDs should be unique within each .h5ad file, Gene IDs should be shared or partially shared across files. Multiple inputs
                        should be listed as a space seperated list of filenames. (default: None)
  -im rna/atac/mc [rna/atac/mc ...], --input_modalities rna/atac/mc [rna/atac/mc ...]
                        (list of str; required) Data modalities chosen from 'rna', 'atac', or 'mc'. This should be listed in the same order as
                        input_datasets. (default: None)
  -f xx.h5ad [xx.h5ad ...], --feature_datasets xx.h5ad [xx.h5ad ...]
                        (list of str; required) Dataset(s) whose features all other datasets will impute into. Enter multiple datasets as a
                        space-separated list of filenames. The features of these datasets will be the features kept in the output imputed data
                        table.", (default: None)
  -o DIR, --output_dir DIR
                        (str) Directory to store output files (default: ./results)
  -op OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                        (str) The output files will contain this prefix. (default: SingleCellFusion)
  --num_pcs NUM_PCS     (integer) Number of Principal Components to keep for each dataset for smoothing and for clustering/embedding after
                        imputation. (default: 50)
  --smoothing_fractions SMOOTHING_FRACTIONS
                        (list of floats) A list of three values between 0 to 1 that controls the relative contribution from the cell itself vs.
                        its neighbors in within-dataset smoothing, specified for 'rna', 'atac', 'mc' data, respectively. (default: [0.7, 0.1,
                        0.9])
  --nearest_neighbors NEAREST_NEIGHBORS
                        (integer) Number of nearest neighbors used to impute data (default: 20)
  --relaxation RELAXATION
                        (integer) A value between 1 to infinity. It specifies how much the datasets should be assumed to be structurally
                        similar. relaxation=1 enforces a hard limit that every cell receives equal number of nearest neighbors
                        relaaxation=infinity approaches traditional kNN (default: 3)
  --leiden_n_neighbors LEIDEN_N_NEIGHBORS
                        (integer) Number of nearest neighbors to form in the integrated space, the resulting nearest neighbor graph is used for
                        Leiden clustering. It is passed into the python package leidenalg. (default: 30)
  --leiden_resolutions LEIDEN_RESOLUTIONS
                        (list of floats) A list of resolutions to be used for Leiden Clustering. It is passed into the python package leidenalg.
                        (default: [0.1, 0.2, 0.4, 0.8])
  --umap_n_neighbors UMAP_N_NEIGHBORS
                        (integer) Number of neighbors for UMAP. It is passed into the python package umap.UMAP(n_neighbors). (default: 60)
  --umap_min_dist UMAP_MIN_DIST
                        (float) Minimum distance for UMAP. It is passed into the python package umap.UMAP(min_dist). (default: 0.5)
```

### Example: integrate L5 ET cells from four data modalities from the mouse primary motor cortex
contains an example of integrating the layer 5 ET (L5 ET) neurons from 4 different datasets from the mouse primary motor cortex. The example includes the organized datasets, code, and results, which could be used as a template for other similar tasks.

```
cd ./example-MOp_L5ET
# run SingleCellFusion
./run_scf.sh
# visualize results
jupyter notebook visualize_results.ipynb
```
![]()

