#define ARMA_DONT_USE_WRAPPER

#include <armadillo>
#include <iostream>

#include "standardize.hpp"

arma::vec standardize(arma::vec training, arma::vec test){
    test -= mean(training);
    test /= stddev(training);
    return(test);
}
