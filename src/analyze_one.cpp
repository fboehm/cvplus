#define ARMA_DONT_USE_WRAPPER

#include <armadillo>
#include <iostream>
#include <string>


#include "analyze_one.hpp"
#include "read_inputs.hpp"

//GOAL: analyze_one reads in a single DBSLMM output file - for a single chromosome - fold pair
// then
// reads in the bed file for the corresponding chromosome and calculates the matrix product 
// genotypes * effects

analyze_one(std::string DBSLMM_output_file, std::string bed_file, std::string bim_file){
    // subset effects vector to have only snps in both DBSLMM output file & bim file
    // we'll also use the resulting indicator vector when reading the bed file
    std::vector<std::vector <std::string> > DBSLMM = read_DSBLMM_output(DBSLMM_output_file); // 3 vectors, rs_id, allele, effect
    std::vector<std::vector <std::string> > bim = read_bim_file(bim_file); // 2 vectors, rs_id and allele
    //https://stackoverflow.com/questions/49441588/c-how-to-check-if-contents-of-vector-exist-in-another-vector
    //https://www.geeksforgeeks.org/stdfind_first_of-in-cpp/
    //make the indicator vector for bim snps being in the DBSLMM output file
    std::vector < bool > bim_snp_in_DBSLMM_output = std::find_first_of(DBSLMM[0].begin(), 
                                                                        DBSLMM[0].end(), 
                                                                        bim[0].begin(), 
                                                                        bim[0].end()) != DBSLMM[0].end()
    
}