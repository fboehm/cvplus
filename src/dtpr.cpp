#include <algorithm>
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <sstream>
#include <fstream>
#include <bitset>
#include <numeric>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <cctype>
#include <boost/math/distributions/students_t.hpp>

#include <armadillo>

#include "dtpr.hpp"

using namespace std;
using namespace arma;

// input genotype data
// modify from GEMMA, Xiang Zhou et al.
//' Read genotype data from plink bed file
//' 
//' @param pos position in bim file
//' @param indicator_idv indicator vector of which individuals to include in analysis
//' @param infile ifstream object for reading bim file
//' @param geno_mat genotype matrix
//' @param n_snp number of SNPs in the block
//' @return void
//' @details indicator_idv contains one entry per subject in the plink files. Each entry is either 0 or 1. Those subjects with 1 get their SNP genotype read, while those with zero don't.

void readSNP(const int pos, //integer position (starting from zero) within the bim file, for the first snp in the block
                   const vector<int> &indicator_idv, 
                   ifstream &infile, 
                   mat &geno_mat,
				   const int n_snp){
				   //,
				   //double &maf
				   //) {

	// debug_msg("entered");
	size_t ni_total = indicator_idv.size(), n_bit;// set ni_total to length of indicator_idv 
	if (ni_total % 4 == 0) {
		n_bit = ni_total / 4;
	}
	else {
		n_bit = ni_total / 4 + 1;
	}

	// n_bit, and 3 is the number of magic numbers.
	for (uint snp = 0; snp < n_snp; snp++){
		vec geno = geno_mat.row(snp);
		infile.seekg((pos + snp) * n_bit + 3);//what is seekg? Sets the position of the next character to be extracted from the input stream.
		//seekg takes an object of ifstream class

		// Read genotypes.
		char ch[1];
		bitset<8> b;
		
		double geno_mean = 0.0;
		size_t c = 0, c_idv = 0;
		vector<size_t> geno_miss;
		int freq[3];
		freq[0] = freq[1] = freq[2] = 0;
		for (size_t i = 0; i < n_bit; ++i) {
			infile.read(ch, 1);
			b = ch[0];

			// Minor allele homozygous: 2.0; major: 0.0.
			for (size_t j = 0; j < 4; ++j) {
				if ((i == (n_bit - 1)) && c == ni_total) {
					break;
				}
				if (indicator_idv[c] == 0) {
					c++;
					continue;
				}
				c++;

				if (b[2 * j] == 0) {
					if (b[2 * j + 1] == 0) {
						geno(c_idv) = 2.0;
						geno_mean += 2.0;
						freq[2]++;
					}
					else {
						geno(c_idv) = 1.0;
						geno_mean += 1.0;
						freq[1]++;
					}
				}
				else {
					if (b[2 * j + 1] == 1) {
						geno(c_idv) = 0.0;
						geno_mean += 0.0;
						freq[0]++; 
					}
					else {
						geno_miss.push_back(c_idv);
					}
				}

				c_idv++;
			}
		}
			// max imputation
		// int imp_val = distance(freq, max_element(freq, freq + 3));
		// mean imputation
		geno_mean /= (double)(c_idv - geno_miss.size());
		for (size_t i = 0; i < geno_miss.size(); ++i) 
			geno(geno_miss[i]) = geno_mean;
		double af = 0.5 * sum(geno) / geno.n_elem; 
		geno_mat.row(snp) = geno;
	}
	//maf = min(af, 1.0 - af);

	return ;
}

