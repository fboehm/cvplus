#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include "read_inputs.hpp"
#include "standardize.hpp"
#include "dtpr.hpp"


class PARAM {
public:
	//parameters
    int n_fold; 
	std::string dbslmm_output_file_prefix;
	std::string plink_file_prefix;
	double alpha;
	std::string path_to_indicator_files;
	std::string path_to_true_pheno_file; 
	std::string outpath;
};


void parse_args(int argc, char *argv[], PARAM &cPar);
