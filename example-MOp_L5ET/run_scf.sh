#!/bin/bash

../scripts/SingleCellFusion \
    -i "./datasets/10x_cells_v2.h5ad" \
       "./datasets/smarter_cells.h5ad" \
       "./datasets/smarter_nuclei.h5ad" \
       "./datasets/snmcseq_gene.h5ad" \
    -im "rna" "rna" "rna" "mc" \
    -f "./datasets/10x_cells_v2.h5ad" \
    -op "test_may4" \
    -o "./results"
