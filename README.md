# SingleCellFusion

SingleCellFusion is a computational tool to integrate single-cell transcriptome and epigenome datasets. Code in this repository is used in [Luo et al 2019 BioRxiv](https://www.biorxiv.org/content/10.1101/2019.12.11.873398v1) and in [manuscript in preparation].

For more information and to cite this work:
- [mukamel-lab/SingleCellFusion_Loompy](https://github.com/mukamel-lab/SingleCellFusion_Loompy)--the github repository of an implementation of SingleCellFusion based on [Loompy](http://loompy.org/). (This repository--now called SingleCellFusion_Loompy--is previously called SingleCellFusion.)
- [Luo, C. et al. Single nucleus multi-omics links human cortical cell regulatory genome diversity to disease risk variants. bioRxiv 2019.12.11.873398 (2019) doi:10.1101/2019.12.11.873398](https://www.biorxiv.org/content/10.1101/2019.12.11.873398v1)


Code contributors: [Fangming Xie](mailto:f7xie@ucsd.edu), Wayne I. Doyle, [Ethan Armand](mailto:ejarmand@ucsd.edu)

Contact: [Eran Mukamel](mailto:emukamel@ucsd.edu)

## Installation
Step 1: Clone this repo.
```bash
git clone https://github.com/mukamel-lab/SingleCellFusion.git
cd SingleCellFusion
```

Step 2: Set up a conda environment and install dependent packages. (Skip this step if not needed.)
```bash
conda env create -f environment.yml # create an env named scf_dev
source activate scf_dev
```

## Usage
```./scripts``` contains the main code.

```./example_l5pt``` contains an example of integrating the layer 5 projection-track (L5 PT) neurons from 4 different datasets from the mouse primary motor cortex ([manuscript in preparation]). The example includes the organized datasets, code, and results, which could be used as a template for other similar tasks.
