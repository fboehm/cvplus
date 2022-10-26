#include <vector>
#include <string>

unsigned int sum_vec(std::vector<int> vv);

std::vector <int> convert_bool_to_int(std::vector < bool > boolean_vector);

std::vector <bool> is_in(const std::vector<std::string>& v1, const std::vector<std::string>& v2);

bool is_one(int x);

arma::vec get_indices(arma::vec A);

vector<int> make_ones_and_zeroes_vec(arma::uvec ones_positions, unsigned int length);

arma::mat subset(arma::mat matrix, arma::uvec indices);

arma::vec subset(arma::vec vector, arma::uvec indices);

std::vector<int> make_integer_vector(int start, int end);

arma::uvec convert_string_to_indices(std::vector <std::string> in_string);

