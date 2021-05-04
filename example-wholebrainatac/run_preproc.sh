#!/bin/bash

ga="/cndd/Public_Datasets/BICCN/BICCN2.0_whole_mouse_brain/references/refdata-gex-mm10-2020-A/genes/genes_promoter_2kb_biccn2.0.bed"

../scripts/normalize_and_select_features.py \
    -i "./datasets_pre/CEMBA171206_3C_genes_promo2kb.h5ad" \
    -inorm "tpm" \
    -ga $ga \
    -op "test_preproc_may3" \
    -o "./datasets_processed"