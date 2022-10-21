#define ARMA_DONT_USE_WRAPPER

#include <armadillo>
#include <iostream>
#include <string>


arma::vec read_one_column_file(const std::string filepath){
    std::ifstream infile;
    infile.open(filepath.c_str()); //read mode
    if(infile.fail()){ // checks to see if file opened  
        std::cout << "error - file didn't open" << std::endl; 
        return 1; // no point continuing if the file didn't open...
    }
    std::string line;
    std::vector<std::string> result;
    while(std::getline(infile, line)){ 
        std::vector<std::string> l0 = std::split(line, " "); //split line with space delimiter 
        result.push_back(l0[0]); // append only the first entry in line
    } 
    infile.close(); 
    arma::vec out = convert_string_to_indices(result);  
    return (out); 
}

