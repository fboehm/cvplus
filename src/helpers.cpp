#include <armadillo>
#include <math.h>       /* floor */
#include <algorithm> /* std::sort, std::set_difference */
#include <fstream> //std::ifstream
#include <string>
#include <boost/lexical_cast.hpp>
#include <iterator>
#include <regex>


#include "helpers.hpp"



//' Sum entries in a vector with only nonnegative integer entries 
//' 
//' @param vv vector<int> to be summed
//' @return unsigned int that is the sum of the entries in the vector

unsigned int sum_vec(std::vector<int> vv){
  unsigned int result = std::accumulate(vv.begin(), vv.end(),
                  decltype(vv)::value_type(0));
  return result;
}

//' Convert a boolean vector to an integer vector of 1s and 0s
std::vector <int> convert_bool_to_int(std::vector < bool > boolean_vector){
    std::vector <int> output;
    for (int i = 0; i < boolean_vector.size(); i++){
        output.push_back(int(boolean_vector[i]));
    }
    return (output);
}


//' For every element in v2, is it in v1?
//'
//' @param v1 reference vector
//' @param v2 vector to be examined
//' @return a vector of booleans with same size as v2. 
//' @details if ith element of v2 is in v1 then ith element of output is true

std::vector <bool> is_in(const std::vector<std::string>& v1, const std::vector<std::string>& v2) {
  std::vector <bool> output(v2.size(), false);
  for (size_t i = 0; i < v2.size(); i++){
      for (size_t j = 0; j < v1.size(); j++){
          if (v2[i] == v1[j]){
              output[i] = true;
              break;
          }
      }
  }
  return(output); 
}



bool is_one(int x){
    return x == 1;
}

//' For A a vector of 0s and 1s, find the indices for entries equal to 1
//'
//' @param A Armadillo vector of 1s and 0s
//' @reference https://stackoverflow.com/questions/25846235/finding-the-indexes-of-all-occurrences-of-an-element-in-a-vector
//' @return Armadillo vector of positive indices where 1s are found in A


arma::vec get_indices(arma::vec A){
    arma::vec output;
    std::vector<int>::iterator iter = A.begin();
    while ((iter = std::find_if(iter, A.end(), is_one)) != A.end())
    {
        // Do something with iter
        output.push_back(iter);    
        iter++;
    }
    return (output);
}


//' Make a vector of zeros and ones with ones at positions indicated
//' 
//' @param ones_positions an arma::uvec specifying the positions where the 1's go
//' @param length the total length of the outputted vector
//' @return an unsigned integer vector of zeros and ones
//' @reference https://gallery.rcpp.org/articles/armadillo-subsetting/

vector<int> make_ones_and_zeroes_vec(arma::uvec ones_positions, unsigned int length){
  cout << "starting make_ones_and_zeroes_vec"<<endl;
  cout <<"result length is: " << length << endl;
  arma::vec result;
  result.zeros(length); //fill vector with all zeros
  //construct a vector for replacing zeroes with ones
  arma::vec ones_vector;
  ones_vector.ones(ones_positions.n_elem);
  cout << "ones_vector has length: " << ones_positions.n_elem << endl;
  // replace zeroes with ones
  //result.elem(ones_positions) = ones_vector;
  result.elem(ones_positions) = ones_vector;
  //convert to vector<int>
  vector<int> out = conv_to<vector<int> >::from(result);
  return out;
}


//' Subset a matrix's rows by indices 
//' 
//' @param mat a matrix, eg., of genotypes, for the entire cohort, with one subject per row
//' @param test_indices vector with subject indices to go into test set
//' @return matrix of genotypes for the subsetted collection of subjects

arma::mat subset(arma::mat matrix, arma::uvec indices){
  arma::mat result = matrix.rows(indices);
  return(result);
}

//' Subset a vector by indices
//' 
//' @param vector a vector, arma::vec
//' @param indices arma::vec of indices to indicate which entries to extract 
//' @return vector of values for the subsetted collection of subjects

arma::vec subset(arma::vec vector, arma::uvec indices){
  arma::vec result = vector.elem(indices);
  return(result);
}

//' Construct an integer vector from start to end, for integers start and end
//' 
//' @param start smallest and first integer value
//' @param end largest and last integer value
//' @return integer vector, start, start + 1, ..., end

std::vector<int> make_integer_vector(int start, int end){
  std::vector<int> myVec;
  for( int i = start; i <= end; i++ ) //spacing between value is 1.
    myVec.push_back( i );
  return myVec;
}


//' Convert std::vector <string> to indices
//' 
//' @param string a string vector
//' @return arma::uvec vector, for use as indices in subsetting armadillo matrices or vectors

arma::uvec convert_string_to_indices(std::vector <std::string> in_string){
  vector<int> vectorOfIntegers;
  castContainer(in_string, vectorOfIntegers);
  arma::uvec result = conv_to< arma::uvec >::from(vectorOfIntegers);
  return (result);
} 





