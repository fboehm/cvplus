#define ARMA_DONT_USE_WRAPPER

#include <armadillo>
#include <iostream>
#include <string>
#include <regex>
#include <iterator>
#include <boost/lexical_cast.hpp>



arma::vec read_one_column_file(const std::string filepath, std::string output_type);

std::vector<std::string> split(const std::string str, const std::string regex_str);

arma::vec convert_string_to_indices(std::vector <std::string> in_string);

arma::vec convert_string_to_indices(std::vector <std::string> in_string);

template<typename C1, typename C2> 
void castContainer(const C1& source, C2& destination);

template<typename T, typename T2>
std::vector<T>& operator<<(std::vector<T>& v, T2 t);

std::vector<std::vector <std::string> > read_DSBLMM_output(const std::string filepath);