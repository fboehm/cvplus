#include <vector>
#include <armadillo>
#include <math.h>       /* floor */
#include <algorithm> /* std::sort, std::set_difference */
#include <fstream> //std::ifstream
#include <string>
#include <boost/lexical_cast.hpp>
#include <iterator>
#include <regex>

unsigned int sum_vec(std::vector<int> vv);

std::vector <int> convert_bool_to_int(std::vector < bool > boolean_vector);

std::vector <bool> is_in(const std::vector<std::string>& v1, const std::vector<std::string>& v2);

bool is_one(int x);

std::vector <int> get_indices(std::vector<int> A);

std::vector<int> make_ones_and_zeroes_vec(arma::uvec ones_positions, unsigned int length);

arma::mat subset(arma::mat matrix, arma::uvec indices);

arma::vec subset(arma::vec vector, arma::uvec indices);

std::vector<int> make_integer_vector(int start, int end);

std::vector <int> convert_string_to_indices(std::vector <std::string> in_string);

template<typename C1, typename C2>
void castContainer(const C1& source, C2& destination);

/*template<typename T, typename T2>
std::vector<T>& operator<<(std::vector<T>& v, T2 t);
*/

/*template <typename T>
std::vector<T> operator+(std::vector<T> lhs, const std::vector<T> & rhs);
*/

arma::vec populate_vec(arma::vec small_v, arma::uvec indices, arma::vec big_v);


std::vector <int> add_two_integer_vectors(std::vector <int> a, std::vector <int> b);
