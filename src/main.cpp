#define ARMA_DONT_USE_WRAPPER

#include <armadillo>
#include <iostream>
#include <string>

#include "main.hpp"
#include "analyze_one.hpp"
#include "read_inputs.hpp"
#include "helpers.hpp"


int main(int argc, char *argv[]){
	
    PARAM cPar;
    parse_args(argc, argv, cPar);
    // define residuals_vv to hold 5 vectors of residuals, one per fold
    std::vector <arma::vec> residuals_vv;
    std::vector <std::vector <std::string> > te_indic;
    for (int fold = 1; fold <= cPar.n_fold; fold++){
            //analyze_one_chr_fold_pair();
        std::vector< arma::vec > pgs;
        for (int chr = 1; chr <= 22; chr++){
            // write args to analyze_one_fold_one_chr with cPar contents
            pgs[fold] += analyze_one_fold_one_chr(cPar.dbslmm_output, 
                                            cPar.plink_file_prefix + std::to_string(chr) + std::string(".bed"), 
                                            cPar.plink_file_prefix + std::to_string(chr) + std::string(".bim"),
                                            cPar.path_to_indicator_files + std::string("indicator_training_fold") + std::to_string(fold) + std::string(".txt"), 
                                            cPar.path_to_indicator_files + std::string("indicator_test_fold") + std::to_string(fold) + std::string(".txt")
                                            );
        }
        // read in true pheno values, then get residuals
        
        std::vector <std::string> true_pheno_string = read_one_column_file(cPar.path_to_true_pheno_files + std::string("true_pheno_test_fold") + std::to_string(fold) + std::string(".txt"));
        std::vector <double> true_pheno_double;
        castContainer(true_pheno_double, true_pheno_string);
        // now convert to arma::vec
        arma::vec true_pheno = arma::conv_to<arma::vec>::from(true_pheno_double);
        
        residuals_vv[fold] = abs(true_pheno - pgs[fold]);
        // read test membership indicator file
        te_indic[fold] = read_one_column_file(cPar.path_to_indicator_files + std::string("indicator_test_fold") + std::to_string(fold) + std::string(".txt"));

    }
    //assemble residuals_vv into a single arma::vec for all "training + test" subjects
    
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
        else if (strcmp(argv[i], "--dbslmm_output") == 0){
            ++i;
            str.clear();
            str.assign(argv[i]);
            cPar.dbslmm_output = str;
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

