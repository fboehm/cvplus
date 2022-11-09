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



void parse_args(int argc, char *argv[], PARAM &cPar)
{
    std::string str;

    for (int i = 0; i < argc; ++i)
    {
        if (strcmp(argv[i], "--n_fold") == 0)
        {
            ++i;
            str.clear();
            str.assign(argv[i]);
            cPar.n_fold = atoi(str.c_str());
        }
        else if (strcmp(argv[i], "--dbslmm_output_file_prefix") == 0)
        {
            ++i;
            str.clear();
            str.assign(argv[i]);
            cPar.dbslmm_output_file_prefix = str;
        }
        else if (strcmp(argv[i], "--plink_file_prefix") == 0)
        {
            ++i;
            str.clear();
            str.assign(argv[i]);
            cPar.plink_file_prefix = str;
        }
        else if (strcmp(argv[i], "--alpha") == 0)
        {
            ++i;
            str.clear();
            str.assign(argv[i]);
            cPar.alpha = std::stod(str.c_str());
        }
        else if (strcmp(argv[i], "--path_to_indicator_files") == 0)
        {
            ++i;
            str.clear();
            str.assign(argv[i]);
            cPar.path_to_indicator_files = str;
        }
        else if (strcmp(argv[i], "--path_to_true_pheno_file") == 0)
        {
            ++i;
            str.clear();
            str.assign(argv[i]);
            cPar.path_to_true_pheno_file = str;
        }
        else if (strcmp(argv[i], "--outpath") == 0)
        {
            ++i;
            str.clear();
            str.assign(argv[i]);
            cPar.outpath = str;
        }
                else if (strcmp(argv[i], "--thread_num") == 0)
        {
            ++i;
            str.clear();
            str.assign(argv[i]);
            cPar.thread_num = std::stoi(str.c_str());
        }
        else if (strcmp(argv[i], "--fold_num") == 0)
        {
            ++i;
            str.clear();
            str.assign(argv[i]);
            cPar.fold_num = std::stoi(str.c_str());
        }
        else if (strcmp(argv[i], "--chr_num") == 0)
        {
            ++i;
            str.clear();
            str.assign(argv[i]);
            cPar.chr_num = std::stoi(str.c_str());
        }
    }
    return;
}
