#define ARMA_DONT_USE_WRAPPER

#include <armadillo>
#include <iostream>
#include <string>
#include <regex>
#include <iterator>
#include <boost/lexical_cast.hpp>

#include "read_inputs.hpp"

//' Read a one-column text file into an armadillo vector
//'
//' @param filepath a string containing the full file path to text file with exactly one column
//' @param output_type a length-one string vector containing the output type, either "integer" or "double"
//' @return Armadillo vector with file contents 

arma::vec read_one_column_file(const std::string filepath, std::string output_type){
    std::ifstream infile;
    infile.open(filepath.c_str()); //read mode
    if(infile.fail()){ // checks to see if file opened  
        std::cout << "error - file didn't open" << std::endl; 
        return; // no point continuing if the file didn't open...
    }
    std::string line;
    std::vector<std::string> result;
    while(std::getline(infile, line)){ 
        std::vector<std::string> l0 = split(line, " "); //split line with space delimiter 
        result.push_back(l0[0]); // append only the first entry in line
    } 
    infile.close(); 
    arma::vec out;
    if (output_type == "integer"){
        out = convert_string_to_indices(result);  
    } else {
        out = convert_string_to_doubles(result); 
    }
    return (out); 
}

//' Split a string vector with a regex
//' 
//' @references https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c

std::vector<std::string> split(const std::string str, const std::string regex_str)
{
  std::regex regexz(regex_str);
  std::vector<std::string> list(std::sregex_token_iterator(str.begin(), str.end(), regexz, -1),
                                std::sregex_token_iterator());
  return list;
}

//' Convert std::vector <string> to indices
//' 
//' @param string a string vector
//' @return arma::uvec vector, for use as indices in subsetting armadillo matrices or vectors

arma::vec convert_string_to_indices(std::vector <std::string> in_string){
  std::vector<int> vectorOfIntegers;
  castContainer(in_string, vectorOfIntegers);
  arma::vec result = arma::conv_to< arma::vec >::from(vectorOfIntegers);
  return (result);
} 

//' Convert string vector to doubles armadillo vector 
//'
//' @param stringVector a string vector where each entry is a double, quoted.
//' @return armadillo vector of doubles
arma::vec convert_string_to_doubles(std::vector <std::string> stringVector){
    std::vector<double> doubleVector(stringVector.size());
    std::transform(stringVector.begin(), stringVector.end(), doubleVector.begin(), [](const std::string& val)
    {
        return std::stod(val);
    });
    arma::vec y = arma::conv_to< arma::vec >::from(doubleVector);
    return (y);
}


//' a single string to single integer function
//' @references https://www.py4u.net/discuss/90965

template<typename C1, typename C2>
void castContainer(const C1& source, C2& destination)
{
  typedef typename C1::value_type source_type;
  typedef typename C2::value_type destination_type;
  destination.resize(source.size());
  std::transform(source.begin(), source.end(), destination.begin(), boost::lexical_cast<destination_type, source_type>);
}

template<typename T, typename T2>
std::vector<T>& operator<<(std::vector<T>& v, T2 t)
{
  v.push_back(T(t));
  return v;
}


//' Read a single DBSLMM output file for a single chromosome and a single fold
std::vector<std::vector <std::string> > read_DSBLMM_output(const std::string filepath){
    std::ifstream infile;
    infile.open(filepath.c_str()); //read mode
    if(infile.fail()){ // checks to see if file opened  
        std::cout << "error - file didn't open" << std::endl; 
        return; // no point continuing if the file didn't open...
    }
    std::string line;
    std::vector<std::string> rs_id;
    std::vector<std::string> allele;
    std::vector<std::string> effect;
    while(std::getline(infile, line)){ 
        std::vector<std::string> l0 = split(line, " "); //split line with space delimiter 
        rs_id.push_back(l0[0]); // append the first entry in line
        allele.push_back(l0[1]);
        effect.push_back(l0[3]);// use fourth column as effect
    } 
    infile.close(); 
    std::vector<std::vector <std::string> > out = {rs_id, allele, effect};
    return(out);    
}


//' Read a single bim file 
std::vector<std::vector <std::string> > read_bim_file(const std::string filepath){
    std::ifstream infile;
    infile.open(filepath.c_str()); //read mode
    if(infile.fail()){ // checks to see if file opened  
        std::cout << "error - file didn't open" << std::endl; 
        return; // no point continuing if the file didn't open...
    }
    std::string line;
    std::vector<std::string> rs_id;
    std::vector<std::string> allele;
    while(std::getline(infile, line)){ 
        std::vector<std::string> l0 = split(line, " "); //split line with space delimiter 
        rs_id.push_back(l0[1]); 
        allele.push_back(l0[4]);
    } 
    infile.close(); 
    std::vector<std::vector <std::string> > out = {rs_id, allele};
    return(out); 
}

