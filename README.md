# DeMixT (v 1.14.0)
DeMixT is a frequentist-based method and fast in yielding accurate estimates of cell proportions and compartment-specific expression profiles for two- and three-component deconvolution from heterogeneous tumor samples. 

# Updates

(10/04/2022) Added preprocessing functions and modified vignettes accordingly

(10/04/2022) Included R.h in DeMixT.c to fix compilation error caused by Rprintf.

(06/18/2022) Fixed a bug in ``R/DeMixT_GS.R:130``: changed ``return(class(try(solve(m), silent = T)) == "matrix")`` to ``return(class(try(solve(m), silent = T))[1] == "matrix")`` to make it work properly under both ``R 3.x`` and ``R 4.x``; since under ``R 4.x``, ``class(matrix object)`` returns ``"matrix", "array"``, instead only ``"matrix"``  under ``R 3.x``. Please be aware of this bug for those who installed ``DeMixT``(v1.10.0) from Bioconductor. 


# Installation
The DeMixT package is compatible with Windows, Linux and MacOS. The user can install it from ``Bioconductor``:
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DeMixT")
```

For Linux and MacOS, the user can also install the latest DeMixT from GitHub:
```
if (!require("devtools", quietly = TRUE))
    install.packages('devtools')

devtools::install_github("wwylab/DeMixT")
```

Check if DeMixT is installed successfully:
```
# load package
library(DeMixT)
```

**Note**: DeMixT relies on OpenMP for parallel computing. Starting from R 4.00, R no longer supports OpenMP on MacOS, meaning the user can only run DeMixT with one core on MacOS. We therefore recommend the users to mainly use Linux system for running DeMixT to take advantage of the multi-core parallel computation.

# Use ``DeMixT``
A tutorial is available at [https://wwylab.github.io/DeMixT/](https://wwylab.github.io/DeMixT/).

# Cite ``DeMixT``

[1] Ahn, J. et al. DeMix: Deconvolution for mixed cancer transcriptomes using raw measured data. Bioinformatics 29, 1865–1871 (2013).

[2] Wang, Z. et al. Transcriptome Deconvolution of Heterogeneous Tumor Samples with Immune Infiltration. iScience 9, 451–460 (2018).