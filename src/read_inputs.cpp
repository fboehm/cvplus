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
//' @return standard string vector  

std::vector <std::string> read_one_column_file(const std::string filepath){
    std::ifstream infile;
    infile.open(filepath.c_str()); //read mode
    std::string line;
    std::vector<std::string> result;
    while(std::getline(infile, line)){ 
        std::vector<std::string> l0 = split(line, " "); //split line with space delimiter 
        result.push_back(l0[0]); // append only the first entry in line
    } 
    infile.close(); 
    return (result); 
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







//' Read a single DBSLMM output file for a single chromosome and a single fold
std::vector<std::vector <std::string> > read_DSBLMM_output(const std::string filepath){
    std::ifstream infile;
    infile.open(filepath.c_str()); //read mode
    std::string line;
    std::vector<std::string> rs_id;
    std::vector<std::string> allele;
    std::vector<std::string> effect;
    while(std::getline(infile, line)){ 
        std::vector<std::string> l0 = split(line, " "); //split line with tab delimiter 
        rs_id.push_back(l0[0]); // append the first entry in line
        allele.push_back(l0[1]);
        effect.push_back(l0[3]);// use fourth column as effect
    } 
    infile.close(); 
    std::vector<std::vector <std::string> > out = {rs_id, allele, effect};
    return(out);    
}


//' Read a single bim file 
std::vector <std::string> read_bim_file(const std::string filepath){
    std::ifstream infile;
    infile.open(filepath.c_str()); //read mode
    std::string line;
    std::vector<std::string> rs_id;
    //std::vector<std::string> allele;
    while(std::getline(infile, line)){ 
        std::vector<std::string> l0 = split(line, "\t"); //split line with tab delimiter
        //std::cout << "l0 has length: " << l0.size() << std::endl; 
        rs_id.push_back(l0[1]); 
        //std::cout << "rs_id is: " << l0[1] << std::endl;
        //allele.push_back(l0[4]);
    } 
    infile.close(); 
    return(rs_id); 
}

