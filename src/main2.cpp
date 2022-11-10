#define ARMA_DONT_USE_WRAPPER

#include <armadillo>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include <math.h>       /* floor */
#include <boost/lexical_cast.hpp>
#include <iterator>
#include <regex>
#include <omp.h>

#include "main.hpp"
#include "read_inputs.hpp"
#include "helpers.hpp"
#include "standardize.hpp"
#include "dtpr.hpp"

int main(int argc, char *argv[]){
    // set up cPar
    PARAM cPar;
    parse_args(argc, argv, cPar);
    // define residuals_vv to hold 5 vectors of residuals, one per fold
    std::vector<arma::vec> residuals_vv(cPar.n_fold);
    std::vector<arma::uvec> test_indices_all_folds(cPar.n_fold); // holds 5 test index arma vectors. We need these vectors later when we assemble residuals vector for entire "training + test" set
    std::vector<arma::uvec> training_indices_all_folds(cPar.n_fold);
    std::vector<arma::uvec> short_test_indices_all_folds(cPar.n_fold);
    std::vector<arma::vec> v_pgs(cPar.n_fold);
    // read true phenos for subsequent use inside folds loop
    std::vector<std::string> true_pheno_string = read_one_column_file(cPar.path_to_true_pheno_file);
    //make a vector to indicate missing pheno values, ie "NA", from true_pheno_string
    std::vector <double> true_pheno_double;
    std::vector <int> true_pheno_missingness_indicator;
    for (uint i = 0; i < true_pheno_string.size(); i++){
        if (true_pheno_string[i] != "NA"){
            true_pheno_double.push_back(std::stod(true_pheno_string[i]));
            true_pheno_missingness_indicator.push_back(0);
        } else {
            true_pheno_missingness_indicator.push_back(1);
        }
    }
    arma::vec true_pheno = arma::conv_to<arma::vec>::from(true_pheno_double);
    std::vector <arma::vec> pgs_all_folds(cPar.n_fold);
    std::vector <arma::vec> v_pgs_all_folds(cPar.n_fold);
    arma::vec resids; // resids is a vector with one entry per subject in the fam file
        
    for (uint fold = 0; fold < cPar.n_fold; fold++){
        // read pgs for each chr-fold pair
        arma::vec pgs_arma;
        arma::vec v_pgs_arma;
        
        for (uint chr = 1; chr <= 22; chr++){
            // assemble file name
            std::string pgs_file = cPar.path_to_pgs_files + std::string("pgs_chr") + std::to_string(chr) + std::string("_fold") + std::to_string(fold + 1) + std::string(".txt");
            arma::vec bar;
            bar.load(pgs_file.c_str());
            if (chr == 1){
                pgs_arma = arma::zeros<vec>(bar.n_elem);
            }        
            pgs_arma += bar;
            // verif set
            std::string v_pgs_file = cPar.path_to_pgs_files + std::string("verif_pgs_chr") + std::to_string(chr) + std::string("_fold") + std::to_string(fold + 1) + std::string(".txt");
            arma::vec v_bar;
            v_bar.load(v_pgs_file.c_str());
            if (chr == 1){
                v_pgs_arma = arma::zeros<vec>(v_bar.n_elem);
            }        
            v_pgs_arma += v_bar;
        }
        pgs_all_folds.push_back(pgs_arma);
        v_pgs_all_folds.push_back(v_pgs_arma);
        // read test indicators
        std::string test_indicator_file = cPar.path_to_indicator_files + std::string("indicator_test_fold") + std::to_string(fold + 1) + std::string(".txt");
        std::vector<std::string> test_indic_string = read_one_column_file(test_indicator_file);
        // convert to std::vector <int>
        std::vector<int> test_indic;
        castContainer(test_indic_string, test_indic);
        //get indices from indicator vector
        std::vector<int> test_indices = get_indices(test_indic);
        arma::uvec test_indices_arma = arma::conv_to<arma::uvec>::from(test_indices);
        //
        if (fold == 0){
            resids.set_size(test_indic.size());
            resids.fill(datum::nan);
        }
        //make true_pheno_double as a std::vector <double> without missing values
        // now convert to arma::vec
        arma::vec true_pheno = arma::conv_to<arma::vec>::from(true_pheno_double);
        std::vector <int> short_test_indic;
        for (uint subject = 0; subject < test_indic.size(); subject++){
            if (true_pheno_missingness_indicator[subject] == 0){
                short_test_indic.push_back(test_indic[subject]);
            }
        }
        std::vector <int> short_test_indices = get_indices(short_test_indic);
        arma::uvec short_test_indices_arma = arma::conv_to<arma::uvec>::from(short_test_indices);
        short_test_indices_all_folds.push_back(short_test_indices_arma);
        
        arma::vec foo = abs(true_pheno.elem(short_test_indices_all_folds[fold]) - pgs_all_folds[fold]);
        residuals_vv.push_back(foo);
        // assemble residuals_vv into a single arma::vec for all "training + test" subjects
        // initialize resids vector with NaN values
        populate_vec(residuals_vv[fold], test_indices_all_folds[fold], resids);
        // resids, from above, still has length equal to one entry per subject in teh fam file.
        //  Some entries, right now, are still datum::nan
        //  we extract the entries that are not datum::nan
        // resids2 then gets added to vector of fitted values
        //// we already have the five fitted values for each verification set subject ////
        // need to loop over verif subjects
        // for each subject, we construct a vector with length equal to the number of "test + training" subjects
    }//end loop over folds
    arma::vec resids2 = resids.elem(arma::find_finite(resids));
    std::string verification_indicator_file = cPar.path_to_indicator_files + std::string("indicator_verification.txt");
    std::vector<std::string> verification_indic_string = read_one_column_file(verification_indicator_file);
    std::vector<int> verification_indic;
    castContainer(verification_indic_string, verification_indic);
    std::vector<int> verification_indices = get_indices(verification_indic);
    arma::uvec verification_indices_arma = arma::conv_to<arma::uvec>::from(verification_indices);
 
    arma::vec upper(verification_indices_arma.n_elem);
    arma::vec lower(verification_indices_arma.n_elem);
    for (int v_subject = 0; v_subject < verification_indices_arma.n_elem; v_subject++)
    {
        arma::vec v_fitted(verification_indic.size());
        v_fitted.fill(datum::nan);
        for (int fold = 0; fold < cPar.n_fold; fold++)
        {
            double foo_val = v_pgs[fold](v_subject);
            arma::vec foo(test_indices_all_folds[fold].n_elem); // resids is a vector with one entry per subject in the fam file
            foo.fill(foo_val);
            populate_vec(foo, test_indices_all_folds[fold], v_fitted);
        } // v_fitted now has nan's in it
        arma::vec v_fitted2 = v_fitted.elem(arma::find_finite(v_fitted));
        // calculate fitted + resids vector
        arma::vec arg_plus = v_fitted2 + resids2;
        arma::vec ap_sorted = arma::sort(arg_plus, "ascend");
        arma::vec arg_minus = v_fitted2 - resids2;
        arma::vec am_sorted = arma::sort(arg_minus, "ascend");
        // determine quantile using the sorted vector
        int p_index = ceil((1 - cPar.alpha) * (double) (v_fitted2.n_elem + 1));
        int m_index = floor(cPar.alpha * (double) (v_fitted2.n_elem + 1));
        // since c++ indexes start with zero, we subtract 1 when from p_index and m_index to get the correct entries
        // eg, for the fifth smallest value, we need index 4
        upper(v_subject) = ap_sorted(p_index - 1);
        lower(v_subject) = am_sorted(m_index - 1);
    }
    // specify filename for true pheno vector for verif subjects
    std::string true_pheno_outfile = cPar.outpath + std::string("true_pheno.txt");
    std::string upper_outfile = cPar.outpath + std::string("upper.txt");
    std::string lower_outfile = cPar.outpath + std::string("lower.txt");

    // subset true phenotype vector to verif subjects
    arma::vec true_pheno_subset = subset(true_pheno, verification_indices_arma);
    true_pheno_subset.save(true_pheno_outfile.c_str(), arma::arma_ascii);
    upper.save(upper_outfile.c_str(), arma::arma_ascii);
    lower.save(lower_outfile.c_str(), arma::arma_ascii);
    return 0;
}

