#define ARMA_DONT_USE_WRAPPER

#include <armadillo>
#include <iostream>
#include <string>
#include <regex>
#include <iterator>
#include <boost/lexical_cast.hpp>



std::vector <std::string> read_one_column_file(const std::string filepath);

std::vector<std::string> split(const std::string str, const std::string regex_str);

std::vector<std::vector <std::string> > read_DSBLMM_output(const std::string filepath);

std::vector <std::string> read_bim_file(const std::string filepath);