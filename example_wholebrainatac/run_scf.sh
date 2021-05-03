#!/bin/bash

../scripts/SCF_main.py \
    -i "./datasets/10x_cells_v2.h5ad" "./datasets/snatac.h5ad" \
    -im "rna" "atac" \
    -f "./datasets/10x_cells_v2.h5ad" \
    -op "test_april27" \
    -o "./results"