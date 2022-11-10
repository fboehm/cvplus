#define ARMA_DONT_USE_WRAPPER

#include <armadillo>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include <math.h>       /* floor */
#include <boost/lexical_cast.hpp>
#include <iterator>
#include <regex>
#include <omp.h>

#include "main.hpp"
#include "read_inputs.hpp"
#include "helpers.hpp"
#include "standardize.hpp"
#include "dtpr.hpp"

#pragma omp declare reduction(+ : arma::vec : \
                              omp_out += omp_in) \
                    initializer( omp_priv = arma::zeros<arma::vec>(omp_orig.n_rows))




int main(int argc, char *argv[])
{
    std::cout << "starting main function " << std::endl;
    // set up cPar
    PARAM cPar;
    parse_args(argc, argv, cPar);
    // read indicator files for verification set
    std::string verification_indicator_file = cPar.path_to_indicator_files + std::string("indicator_verification.txt");
    std::vector<std::string> verification_indic_string = read_one_column_file(verification_indicator_file);
    std::vector<int> verification_indic;
    castContainer(verification_indic_string, verification_indic);
    std::vector<int> verification_indices = get_indices(verification_indic);
    arma::uvec verification_indices_arma = arma::conv_to<arma::uvec>::from(verification_indices);
    
/*    std::vector<arma::vec> pgs(cPar.n_fold);
    std::vector<arma::vec> v_pgs(cPar.n_fold);
 */
    // specify max block size in number of SNPs
    uint max_block_size = 10000;
    std::cout << "Starting chr " << cPar.chr_num << std::endl;
    // write args to analyze_one_fold_one_chr with cPar contents
    // set up bed file stream
    std::string bed_file = cPar.plink_file_prefix + std::to_string(cPar.chr_num) + std::string(".bed");
    std::ifstream bed_file_stream(bed_file.c_str(), std::ios::binary);
    std::string bim_file_name = cPar.plink_file_prefix + std::to_string(cPar.chr_num) + std::string(".bim");
    std::vector<std::string> bim = read_bim_file(bim_file_name); // 1 vector, rs_id
    std::string training_indicator_file = cPar.path_to_indicator_files + std::string("indicator_training_fold") + std::to_string(cPar.fold_num) + std::string(".txt");
    std::string test_indicator_file = cPar.path_to_indicator_files + std::string("indicator_test_fold") + std::to_string(cPar.fold_num) + std::string(".txt");
    std::vector<std::string> training_indic_string = read_one_column_file(training_indicator_file);
    std::vector<std::string> test_indic_string = read_one_column_file(test_indicator_file);
    // convert to std::vector <int>
    std::vector<int> training_indic;
    castContainer(training_indic_string, training_indic);
    std::vector<int> test_indic;
    castContainer(test_indic_string, test_indic);
    std::vector<int> training_indices = get_indices(training_indic);
    std::vector<int> test_indices = get_indices(test_indic);
    arma::uvec training_indices_arma = arma::conv_to<arma::uvec>::from(training_indices);
    arma::uvec test_indices_arma = arma::conv_to<arma::uvec>::from(test_indices);
    /*test_indices_all_folds[fold] = test_indices_arma;
    training_indices_all_folds[fold] = training_indices_arma;
    */
    std::string dbslmm_output_fn = cPar.dbslmm_output_file_prefix + std::to_string(cPar.fold_num) + std::string("_chr") + std::to_string(cPar.chr_num) + std::string("_best.dbslmm.txt");
    std::vector<std::vector<std::string>> DBSLMM = read_DSBLMM_output(dbslmm_output_fn); // 3 vectors, rs_id, allele, effect
    // https://stackoverflow.com/questions/49441588/c-how-to-check-if-contents-of-vector-exist-in-another-vector
    uint n_blocks = ceil(DBSLMM[0].size() /(double) max_block_size);
    uint n_snp = DBSLMM[0].size();
    // make the indicator vector for bim snps being in the DBSLMM output file

    std::vector<bool> bim_snp_in_DBSLMM_output = is_in(DBSLMM[0], bim);
    // read one SNP's genotypes for all subjects
    // determine pos value for readSNP function
    // we only read SNPs that are in the DBSLMM output file

    int DBSLMM_snp = 0; // counter to progress through DBSLMM output file
    // make subject indicator to know which subjects to read genotypes of
    std::vector<int> subject_indicator_pre = add_two_integer_vectors(training_indic, test_indic);
    std::vector<int> subject_indicator = add_two_integer_vectors(subject_indicator_pre, verification_indic);
    // determine length of product_vec and v_product_vec
    arma::vec product_vec(sum_vec(test_indic));
    arma::vec v_product_vec(sum_vec(verification_indic));
    
    // https://stackoverflow.com/questions/28607912/sum-values-of-2-vectors
    //https://stackoverflow.com/questions/11773115/parallel-for-loop-in-openmp
    for (uint block_num = 0; block_num < n_blocks; block_num++){
        // set snp_block_size
        uint snp_block_size;
        if (block_num + 1 != n_blocks){
            snp_block_size = max_block_size;
        } else {
            snp_block_size = n_snp % max_block_size;
        }
        // initialize geno_mat
        
        arma::mat geno_mat = arma::zeros<mat>(sum_vec(subject_indicator), snp_block_size);
        arma::mat test_geno_mat = arma::zeros<mat>(sum_vec(test_indic), snp_block_size);
        arma::mat verif_geno_mat = arma::zeros<mat>(sum_vec(verification_indic), snp_block_size);
        arma::vec geno = arma::zeros<vec>(sum_vec(subject_indicator));
        uint bim_start_point = block_num * max_block_size;
        uint bim_end_point = bim_start_point + snp_block_size;
        //need to iterate bim_snp and DBSLMM_snp in this loop; however
        // want to put DBSLMM_snp in the for () and have bim_snp as an extra counter
        #pragma omp parallel for num_threads(cPar.thread_num) reduction(+:product_vec,v_product_vec)
        for (uint bim_snp = bim_start_point; bim_snp < bim_end_point; bim_snp++)
        {
            // check if SNP from bim is in DBSLMM file
            if (bim_snp_in_DBSLMM_output[bim_snp])
            {
                readSNP(bim_snp, subject_indicator, bed_file_stream, geno);
                // partition geno into training, test, and verif sets
                arma::vec training_geno = subset(geno, training_indices_arma);
                arma::vec test_geno = subset(geno, test_indices_arma);
                arma::vec verif_geno = subset(geno, verification_indices_arma);
                arma::vec effects(snp_block_size);
                // standardize verif_geno & test_geno
                for (uint col = 0; col < geno_mat.n_cols; col++){
                    //define training_geno, test_geno & verif_geno
                    arma::vec test_geno_std = standardize(training_geno, test_geno);
                    arma::vec verif_geno_std = standardize(training_geno, verif_geno);
                    double dd = std::stod(DBSLMM[2][DBSLMM_snp]);
                    // put back into matrix
                    test_geno_mat.col(col) = test_geno_std;
                    verif_geno_mat.col(col) = verif_geno_std;
                    effects(col) = dd; 
                }
                // multiply standardized genotypes by DBSLMM effect for that snp
                // multiply the standardized test set genotypes by effect for that snp
                product_vec += test_geno_mat * effects;
                v_product_vec += verif_geno_mat * effects;
                DBSLMM_snp++; // advance counter for snps in DBSLMM file
                // note that we assume that snps in DBSLMM file is a subset of snps in bim file
            }
        }
            
    } // end loop over blocks
    
    bed_file_stream.close();
    return 0;
}

