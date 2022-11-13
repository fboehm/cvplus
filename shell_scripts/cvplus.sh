#!/bin/bash

#SBATCH --partition=mulan,main
#SBATCH --time=48:00:00
#SBATCH --job-name=cvpluscpp
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=4
#SBATCH --array 1-110
#SBATCH --output=/net/mulan/home/fredboe/research/ukb-intervals/cluster_outputs/cvplus_%a.out
#SBATCH --error=/net/mulan/home/fredboe/research/ukb-intervals/cluster_outputs/cvplus_%a.err



# path to the executable
cvplus_path=/net/mulan/home/fredboe/research/cvplus/src/cvplus1

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

let k=0

for fold_num in `seq 1 5`; do
    for chr_num in `seq 1 22`; do
        let k=${k}+1
        if [ ${k}==${SLURM_ARRAY_TASK_ID} ]; then
            # cvplus1 command
            ${cvplus_path} --n_fold 5 --dbslmm_output_file_prefix ~/research/ukb-intervals/05_internal_c/pheno1/DBSLMM/summary_ukb_cross \
                            --plink_file_prefix /net/mulan/disk2/yasheng/predictionProject/plink_file/ukb/chr \
                            --alpha 0.1 \
                            --path_to_indicator_files ~/research/ukb-intervals/03_subsample/continuous/pheno1/indicator_files/ \
                            --path_to_true_pheno_file ~/research/ukb-intervals/03_subsample/continuous/pheno1/true_pheno.txt \
                            --outpath ~/research/ukb-intervals/05_internal_c/pheno1/cvplus_outputs/ \
                            --thread_num ${SLURM_CPUS_PER_TASK} \
                            --fold_num ${fold_num} \
                            --chr_num ${chr_num} \
                            --path_to_pgs_files  ~/research/ukb-intervals/05_internal_c/pheno1/cvplus1_outputs/
        fi
    done
done
    


