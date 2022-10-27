#include <string>


class PARAM {
public:
	//parameters
    int n_fold; 
	std::string dbslmm_output;
	std::string plink_file_prefix;
	double alpha;
	std::string path_to_indicator_files;
	std::string path_to_true_pheno_files; 
};
