#!/bin/bash

#SBATCH --partition=mulan,main
#SBATCH --time=24:00:00
#SBATCH --job-name=cvpluscpp
#SBATCH --mem-per-cpu=16G
#SBATCH --output=/net/mulan/home/fredboe/research/ukb-intervals/cluster_outputs/cvpluscpp.out
#SBATCH --error=/net/mulan/home/fredboe/research/ukb-intervals/cluster_outputs/cvpluscpp.err



# path to the executable
cvplus_path=/net/mulan/home/fredboe/research/cvplus/src/cvplus

# cvplus command
${cvplus_path} --n_fold 5 --dbslmm_output_file_prefix ~/research/ukb-intervals/05_internal_c/pheno1/DBSLMM/summary_ukb_cross --plink_file_prefix /net/mulan/disk2/yasheng/predictionProject/plink_file/ukb/chr --alpha 0.1 --path_to_indicator_files ~/research/ukb-intervals/03_subsample/continuous/pheno1/indicator_files/ --path_to_true_pheno_file ~/research/ukb-intervals/03_subsample/continuous/pheno1/true_pheno.txt --outpath ~/research/ukb-intervals/05_internal_c/pheno1/cvplus_outputs/ 

