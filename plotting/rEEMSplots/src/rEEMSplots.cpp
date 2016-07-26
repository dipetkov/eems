
#include <RcppEigen.h>
#include <Rcpp.h>
using namespace Rcpp;

void rcppstandardize_rates(const Eigen::VectorXd &tiles, const Eigen::VectorXd &rates,
                           const Eigen::MatrixXd &seeds, const Eigen::MatrixXd &marks,
                           const std::string &distm, Eigen::VectorXd &zvals, 
                           Eigen::VectorXd &prgt0, Eigen::VectorXd &prlt0);
void rcppdont_standardize_rates(const Eigen::VectorXd &tiles, const Eigen::VectorXd &rates,
                                const Eigen::MatrixXd &seeds, const Eigen::MatrixXd &marks,
                                const std::string &distm, Eigen::VectorXd &zvals);
Eigen::MatrixXd euclidean_dist(const Eigen::MatrixXd &X, const Eigen::MatrixXd &Y);
Eigen::MatrixXd greatcirc_dist(const Eigen::MatrixXd &X, const Eigen::MatrixXd &Y);


RcppExport SEXP rEEMSplots__rcppstandardize_rates(SEXP tiles_, SEXP rates_, SEXP seeds_, SEXP marks_, SEXP distm_) {
    BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type tiles(tiles_);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type rates(rates_);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type seeds(seeds_);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type marks(marks_);
    Rcpp::traits::input_parameter< const std::string& >::type distm(distm_);
    Eigen::VectorXd zvals, prgt0, prlt0;
    rcppstandardize_rates(tiles, rates, seeds, marks, distm, zvals, prgt0, prlt0);
    return Rcpp::List::create(Rcpp::Named("zvals") = zvals,
                              Rcpp::Named("prgt0") = prgt0,
                              Rcpp::Named("prlt0") = prlt0);
    END_RCPP
}
RcppExport SEXP rEEMSplots__rcppdont_standardize_rates( SEXP tiles_, SEXP rates_, SEXP seeds_, SEXP marks_, SEXP distm_) {
    BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type tiles(tiles_);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type rates(rates_);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type seeds(seeds_);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type marks(marks_);
    Rcpp::traits::input_parameter< const std::string& >::type distm(distm_);
    Eigen::VectorXd zvals;
    rcppdont_standardize_rates(tiles, rates, seeds, marks, distm, zvals);
    return Rcpp::List::create(Rcpp::Named("zvals") = zvals);
    END_RCPP
}
RcppExport SEXP rEEMSplots__euclidean_dist( SEXP X_, SEXP Y_) {
    BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X(X_);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Y(Y_);
    __result = Rcpp::wrap(euclidean_dist(X, Y));
    return __result;
    END_RCPP    
}
RcppExport SEXP rEEMSplots__greatcirc_dist( SEXP X_, SEXP Y_) {
    BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X(X_);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Y(Y_);
    __result = Rcpp::wrap(greatcirc_dist(X, Y));
    return __result;
    END_RCPP    
}
