#include <armadillo>
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>

arma::vec analyze_one_fold_one_chr(std::string DBSLMM_output_file, 
            std::string bed_file, 
            std::string bim_file,
            std::string training_indicator_file, 
            std::string test_indicator_file);
 