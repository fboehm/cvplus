#define ARMA_DONT_USE_WRAPPER

#include <armadillo>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include "main.hpp"
#include "analyze_one.hpp"
#include "read_inputs.hpp"
#include "helpers.hpp"
#include "standardize.hpp"
#include "dtpr.hpp"


int main(int argc, char *argv[]){
	// set up cPar
    PARAM cPar;
    parse_args(argc, argv, cPar);
    // read indicator files for verification set
    std::string verification_indicator_file = cPar.path_to_indicator_files + std::string("indicator_verification.txt");
    std::vector <std::string> verification_indic_string = read_one_column_file(verification_indicator_file);
    std::vector<int> verification_indic;
    castContainer(verification_indic_string, verification_indic);
    std::vector <int> verification_indices = get_indices(verification_indic);
    arma::uvec verification_indices_arma = arma::conv_to< arma::uvec >::from(verification_indices);

    // define residuals_vv to hold 5 vectors of residuals, one per fold
    std::vector <arma::vec> residuals_vv;
    std::vector < arma::uvec > test_indices_all_folds; // holds 5 test index arma vectors. We need these vectors later when we assemble residuals vector for entire "training + test" set
    std::vector < arma::uvec > training_indices_all_folds;
    // read true phenos for subsequent use inside folds loop
    std::vector <std::string> true_pheno_string = read_one_column_file(cPar.path_to_true_pheno_files + std::string("true_pheno.txt"));
    std::vector <double> true_pheno_double;
    castContainer(true_pheno_double, true_pheno_string);
    // now convert to arma::vec
    arma::vec true_pheno = arma::conv_to<arma::vec>::from(true_pheno_double);
    
    std::vector< arma::vec > pgs;
    for (int chr = 1; chr <= 22; chr++){
        // write args to analyze_one_fold_one_chr with cPar contents
        //set up bed file stream
        std::string bed_file = cPar.plink_file_prefix + std::to_string(chr) + std::string(".bed");
        std::ifstream bed_file_stream(bed_file.c_str(), std::ios::binary);
        std::vector<std::vector <std::string> > bim = read_bim_file(cPar.plink_file_prefix + std::to_string(chr) + std::string(".bim")); // 2 vectors, rs_id and allele
        // get indices from indicator vectors
        for (int fold = 1; fold <= cPar.n_fold; fold++){
            std::string training_indicator_file = cPar.path_to_indicator_files + std::string("indicator_training_fold") + std::to_string(fold) + std::string(".txt");
            std::string test_indicator_file = cPar.path_to_indicator_files + std::string("indicator_test_fold") + std::to_string(fold) + std::string(".txt");
            std::vector <std::string> training_indic_string = read_one_column_file(training_indicator_file);
            std::vector <std::string> test_indic_string = read_one_column_file(test_indicator_file);
            // convert to std::vector <int>
            std::vector<int> training_indic;
            castContainer(training_indic_string, training_indic);
            std::vector<int> test_indic;
            castContainer(test_indic_string, test_indic);
            std::vector <int> training_indices = get_indices(training_indic);
            std::vector <int> test_indices = get_indices(test_indic);
            arma::uvec training_indices_arma = arma::conv_to< arma::uvec >::from(training_indices);
            arma::uvec test_indices_arma = arma::conv_to< arma::uvec >::from(test_indices);
            test_indices_all_folds[fold] = test_indices_arma;
            training_indices_all_folds[fold] = training_indices_arma;
            
            // subset effects vector to have only snps in both DBSLMM output file & bim file
            // we'll also use the resulting indicator vector when reading the bed file
            std::vector<std::vector <std::string> > DBSLMM = read_DSBLMM_output(cPar.dbslmm_output_file_prefix + std::to_string(fold) + std::string("_chr") + std::to_string(chr) + std::string("_best.dbslmm.txt")); // 3 vectors, rs_id, allele, effect
            //https://stackoverflow.com/questions/49441588/c-how-to-check-if-contents-of-vector-exist-in-another-vector
            //make the indicator vector for bim snps being in the DBSLMM output file
            std::vector < bool > bim_snp_in_DBSLMM_output = is_in(DBSLMM[0], bim[0]);
            // read one SNP's genotypes for all subjects
            // determine pos value for readSNP function
            // we only read SNPs that are in the DBSLMM output file
            
            arma::vec geno;
            int DBSLMM_snp = 0;// counter to progress through DBSLMM output file
            //make subject indicator to know which subjects to read genotypes of
            std::vector < int > subject_indicator = training_indic + test_indic + verification_indic;
            arma::vec product_vec;
            //https://stackoverflow.com/questions/28607912/sum-values-of-2-vectors
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
                    double dd = std::stod(DBSLMM[2][DBSLMM_snp]);
                    // multiply the standardized test set genotypes by effect for that snp
                    product_vec += test_geno_std *(double) dd;            
                    DBSLMM_snp++;//advance counter for snps in DBSLMM file
                    //note that we assume that snps in DBSLMM file is a subset of snps in bim file 
                }
            }
            //store product_vec
            pgs[fold] += product_vec;

        }// end loop over folds
    } // end loop over chr
    //loop over folds
    for (int fold = 1; fold <= cPar.n_fold; fold++){
        residuals_vv[fold] = abs(true_pheno.elem(test_indices_all_folds[fold]) - pgs[fold]);
    }
    //assemble residuals_vv into a single arma::vec for all "training + test" subjects
    //initialize resids vector with NaN values
    arma::vec resids(verification_indic.size()); // resids is a vector with one entry per subject in the fam file
    resids.fill(datum::nan);
    for (int fold = 1; fold <= cPar.n_fold; fold++){
        populate_vec(residuals_vv[fold], test_indices_all_folds[fold], resids);
    }
    //resids, from above, still has length equal to one entry per subject in teh fam file. 
    // Some entries, right now, are still datum::nan
    // we need to extract the entries that are not missing into a new vector
    // alternatively, we might remove the entries that ARE datum::nan
    

    return 0;
}

void parse_args(int argc, char *argv[], PARAM &cPar){
    std::string str;         
        


    for (int i = 0; i < argc; ++i){
        if (strcmp(argv[i], "--n_fold") == 0){
            ++i;
            str.clear();
            str.assign(argv[i]);
            cPar.n_fold = atoi(str.c_str());

        }
        else if (strcmp(argv[i], "--dbslmm_output_file_prefix") == 0){
            ++i;
            str.clear();
            str.assign(argv[i]);
            cPar.dbslmm_output_file_prefix = str;
        }
        else if (strcmp(argv[i], "--plink_file_prefix") == 0){
            ++i;
            str.clear();
            str.assign(argv[i]);
            cPar.plink_file_prefix = str;
        }
        else if (strcmp(argv[i], "--alpha") == 0){
            ++i;
            str.clear();
            str.assign(argv[i]);
            cPar.alpha = std::stod(str.c_str());
        } 
        else if (strcmp(argv[i], "--path_to_indicator_files") == 0){
            ++i;
            str.clear();
            str.assign(argv[i]);
            cPar.path_to_indicator_files = str;

        }
        else if (strcmp(argv[i], "--path_to_true_pheno_files") == 0){
            ++i;
            str.clear();
            str.assign(argv[i]);
            cPar.path_to_true_pheno_files = str;

        }
    }
    return; 
}

