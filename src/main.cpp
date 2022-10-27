#define ARMA_DONT_USE_WRAPPER

#include <armadillo>
#include <iostream>
#include <string>

#include "main.hpp"

int main(int argc, char *argv[]){
	
    PARAM cPar;
    parse_args(argc, argv, cPar);


    // loop over chr 
    for (int chr = 1; chr <= 22; chr++){
        //loop over folds
        for (int fold = 1; fold <= cPar.n_fold; fold++){
            //analyze_one_chr_fold_pair();
        }
    }
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
    }
    return; 
}

