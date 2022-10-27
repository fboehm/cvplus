#define ARMA_DONT_USE_WRAPPER

#include <armadillo>
#include <iostream>
#include <string>
#include<algorithm>


#include "analyze_one.hpp"
#include "read_inputs.hpp"
#include "dtpr.hpp"
#include "helpers.hpp"
#include "standardize.hpp"

//GOAL: analyze_one reads in a single DBSLMM output file - for a single chromosome - fold pair
// then
// reads in the bed file for the corresponding chromosome and calculates the matrix product 
// genotypes * effects

//' Calculate test set residuals for one test set for a single fold
//' 
//' @param DBSLMM_output_file file path for a single DBSLMM output file for a single fold
//' @param bed_file plink bed file path
//' @param bim_file plink bim file path
//' @param training_indicator_file file path to a text file containing a single column of 1s and 0s to indicate training set membership for this fold
//' @param test_indicator_file  file path to a text file containing a single column of 1s and 0s to indicate test set membership for this fold
//' @param test_pheno_file file path to a text file containing a single column of true phenotype values for the test set for this fold
//' @return 

arma::vec analyze_one_fold(std::string DBSLMM_output_file, 
                std::string bed_file, 
                std::string bim_file,
                std::string training_indicator_file, 
                std::string test_indicator_file, 
                std::string test_pheno_file){

    // read indicator files
    std::vector <std::string> training_indic_string = read_one_column_file(training_indicator_file);
    std::vector <std::string> test_indic_string = read_one_column_file(test_indicator_file);
    // convert to std::vector <int>
    std::vector<int> training_indic;
    castContainer(training_indic_string, training_indic);
    std::vector<int> test_indic;
    castContainer(test_indic_string, test_indic);
    

    // get indices from indicator vectors
    std::vector <int> training_indices = get_indices(training_indic);
    std::vector <int> test_indices = get_indices(test_indic);
    arma::uvec training_indices_arma = arma::conv_to< arma::uvec >::from(training_indices);
    arma::uvec test_indices_arma = arma::conv_to< arma::uvec >::from(test_indices);
    // subset effects vector to have only snps in both DBSLMM output file & bim file
    // we'll also use the resulting indicator vector when reading the bed file
    std::vector<std::vector <std::string> > DBSLMM = read_DSBLMM_output(DBSLMM_output_file); // 3 vectors, rs_id, allele, effect
    std::vector<std::vector <std::string> > bim = read_bim_file(bim_file); // 2 vectors, rs_id and allele
    //https://stackoverflow.com/questions/49441588/c-how-to-check-if-contents-of-vector-exist-in-another-vector
    //make the indicator vector for bim snps being in the DBSLMM output file
    std::vector < bool > bim_snp_in_DBSLMM_output = is_in(DBSLMM[0], bim[0]);
    // read one SNP's genotypes for all subjects
    // determine pos value for readSNP function
    // we only read SNPs that are in the DBSLMM output file
    ifstream bed_file_stream(bed_file.c_str(), ios::binary);
    arma::vec geno;
    int DBSLMM_snp = 0;
    std::vector < int > subject_indicator = training_indic + test_indic;
    for (int bim_snp = 0; bim_snp < bim_snp_in_DBSLMM_output.size(); bim_snp++){
        //check if SNP from bim is in DBSLMM file
        if (bim_snp_in_DBSLMM_output[bim_snp]){
            readSNP(bim_snp, subject_indicator, bed_file_stream, geno);
            //partition geno into training, test, and verif sets
            arma::vec training_geno = subset(geno, training_indices_arma);
            arma::vec test_geno = subset(geno, test_indices_arma);
            // standardize verif_geno & test_geno
            arma::vec test_geno_std = standardize(training_geno, test_geno);
            // multiply standardized genotypes by DBSLMM effect for that snp
            


            DBSLMM_snp++;//advance counter for snps in DBSLMM file
            //note that we assume that snps in DBSLMM file is a subset of snps in bim file 
        }
    }
    // read in pheno values for test set for this fold
    std::vector <std::string> test_pheno_string = read_one_column_file(test_pheno_file.c_str());
    // convert string to numeric
    

    // Calculate residuals
    arma::vec residuals = test_pheno - pgs;
    return(residuals);
}
