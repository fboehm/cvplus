# cvplus
C++ code to implement CV+ (after Barber et al., 2021) for use with polygenic scores


## Inputs

1. n_fold: the number of folds for CV+
1. dbslmm_output_file_prefix: DBSLMM output files, one per chromosome per fold for a single trait. These need to be in a common directory because `cvplus` input is the file prefix. 
1. 
1. plink_file_prefix: plink binary files, one set per chromosome. We assume that the plink files are named by chromosome, chr1.bed, chr1.bim, chr1.fam, etc., and in a common directory.
1. training & test & verification set indicator files - appropriately named and in a common subdirectory
1. output directory, where the true pheno, lower, and upper files text files will be saved.
1. 
