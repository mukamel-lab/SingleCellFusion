FAQs
================
SingleCellFusion is under active development and function names and parameters will continue to be
changed until a stable release is reached. In the interim, we have provided some answers to common
questions and problems that can occur when using SingleCellFusion.

Why do you use loom files and how do I make one?
-------------------------------------------------
The loom file format allows SingleCellFusion to have a low memory footprint when analyzing large data
sets (such as 10x Genomics scRNA-seq data), and keep all of the meta-data in one centralized location.
The loompy package was developed by the Sten Linnarsson group and has excellent documentation at
`loompy.org <http://loompy.org/>`_.

Within a loom file features are stored in rows, and cells in columns. As an example to create a loom file,
say you have a pandas dataframe (df) in which the features are in rows and cells are in columns. The index of
this dataframe contains the unique feature IDs and the column header contains unique cell IDs. A loom file
can be generated with the following code::

    import loompy
    loompy.create(filename=filename,
                  layers={'':df.values},
                  row_attrs={'Accession:df.index.values},
                  col_attrs={'CellID:df.columns.values})

Why is my code using so much memory, even with the low_mem flag?
-----------------------------------------------------------------
Access of loom files is performed in batches to reduce the memory overload. In the basic recipe for
SingleCellFusion (pairwise_impute) the size of these batches is controlled by the parameter batch_x and
batch_y. If you are having memory issues, try reducing the size of these values to reduce your memory
overhead.

Why is my code using so many threads?
--------------------------------------
This is a problem with Numpy (see this `issue <https://github.com/numpy/numpy/issues/11826>`_). You can solve this
issue in two ways.

The easiest way is to run the following lines on your command line before running any Python scripts or notebooks::

    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export VECLIB_MAXIMUM_THREADS=1
    export NUMEXPR_NUM_THREADS=1

You can also run the below code in the first cell of a Python notebook or the beginning or a Python script. It must
be run before importing any other packages (including MoP)::

    import os
    os.environ["OMP_NUM_THREADS"] = '1'
    os.environ["OPENBLAS_NUM_THREADS"] = '1'
    os.environ["MKL_NUM_THREADS"] = '1'
    os.environ["VECLIB_MAXIMUM_THREADS"] = '1'
    os.environ["NUMEXPR_NUM_THREADS"] = '1'

You can specify the maximum number of threads that you want to use in that script or notebook by changing the value
from 1 to your desired integer.

This information came from this `StackOverflow question
<https://stackoverflow.com/questions/30791550/limit-number-of-threads-in-numpy>`_.

Why is my code running slow with the low_mem flag?
--------------------------------------------------
Although the loom file format has a number of benefits, the access and processing of data in the file
will get progressively slower as more data is added to the file. If you are finding that your code is
running too slow it can be helpful to make a second loom file containing just the relevant data for running
SingleCellFusion.

Another cause of slow code is that the batch size for processing code (see "Why is my code using so much
memory?" above) is too small. If you are not having memory issues, we recommend increasing the batch size
to speed up the code.

Why am I not finding many neighbors?
-------------------------------------
If you expect that similar cell types should be present in both data sets, this could be due to
the sparseness of your data. We have found that if you first smooth your data (we highly
recommend using `MAGIC <https://github.com/KrishnaswamyLab/MAGIC>`_. You can then use the
smoothed data to find nearest neighbors, and impute on the observed data. A tutorial using our
loom-based method of smoothing will be uploaded soon.


Is SingleCellFusion just for integrating data from different sequencing modalities?
-----------------------------------------------------------------------------------
No, theoretically this pipeline could be applied to integration across species or to find common cell
types across different research studies using the same sequencing technology. This is an active area
of development.

What happens if a cell type is present in only one modality?
-------------------------------------------------------------
In our experience, this situation is easily detectable. If the analysis is only performed on direct
mutual nearest neighbors, these cells will not make nearest neighbors and will be dropped from the analysis.
If the imputation is performed with the rescue, these cells will still not make mutual nearest neighbors.
Their imputed counts will then come from their mutual nearest neighbors within their own data set. These
imputed counts will not be similar to any observed counts, and these cells will self-segregate into their
own clusters and will be visually separate on a tSNE or uMAP embedding. For the kNN method, these cells will
make weak connections with a number of different cell types. This will lead to the imputation of counts that
are not similar to any observed counts, also leading to segregation into unique clusters.


