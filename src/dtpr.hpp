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

using namespace std;
using namespace arma;



void readSNP(const int pos, 
                   const vector<int> &indicator_idv, 
                   ifstream &infile, 
                   mat &geno_mat, 
                   const int n_snp);
                   