#include <armadillo>
#include <vector>


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