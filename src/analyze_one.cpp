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

analyze_one_fold(std::string DBSLMM_output_file, 
                std::string bed_file, 
                std::string bim_file,
                std::string training_indicator_file, 
                std::string test_indicator_file, 
                std::string verification_indicator_file){
    // read indicator files
    arma::vec training_indic = read_one_column_file(training_indicator_file, "integer");
    arma::vec test_indic = read_one_column_file(test_indicator_file, "integer");
    arma::vec verif_indic = read_one_column_file(verification_indicator_file, "integer");
    arma::vec training_indices = get_indices(training_indic);
    arma::vec test_indices = get_indices(test_indic);
    arma::vec verif_indices = get_indices(verif_indic);
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
    ifstream bed_file_stream(bed_str.c_str(), ios::binary);
    arma::vec geno;
    int DBSLMM_snp = 0;

    for (int bim_snp = 0; bim_snp < bim_snp_in_DBSLMM_output.size(); bim_snp++){
        //check if SNP from bim is in DBSLMM file
        if (bim_snp_in_DBSLMM_output[bim_snp]){
            readSNP(bim_snp, subject_indicator, bed_file_stream, geno);
            //partition geno into training, test, and verif sets
            arma::vec training_geno = subset(geno, training_indices);
            arma::vec verif_geno = subset(geno, verif_indices);
            arma::vec test_geno = subset(geno, test_indices);
            // standardize verif_geno & test_geno
            arma::vec verif_geno_std = standardize(training_geno, verif_geno);
            arma::vec test_geno_std = standardize(training_geno, test_geno);
            // 
            DBSLMM_snp++;//advance counter for snps in DBSLMM file
            //note that we assume that snps in DBSLMM file is a subset of snps in bim file 
        }
    }
}
