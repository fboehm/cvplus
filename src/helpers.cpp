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

