
#include <RcppEigen.h>
#include <Rcpp.h>
using namespace Rcpp;

// rcppstandardize_rates
Eigen::MatrixXd rcppstandardize_rates(const Eigen::VectorXd &tilex, const Eigen::VectorXd &rates,
				      const Eigen::VectorXd &xseed, const Eigen::VectorXd &yseed,
				      const Eigen::MatrixXd &marks, const Eigen::VectorXd &nmrks,
				      const std::string &distm);
Eigen::MatrixXd euclidean_dist(const Eigen::MatrixXd &X, const Eigen::MatrixXd &Y);



RcppExport SEXP rEEMSplots__rcppstandardize_rates( SEXP tiles_, SEXP rates_, SEXP xseed_, SEXP yseed_,
						   SEXP marks_, SEXP nmrks_, SEXP distm_) {
BEGIN_RCPP
  SEXP Zmrks_;
  {
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type tiles(tiles_);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type rates(rates_);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type xseed(xseed_);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type yseed(yseed_);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type nmrks(nmrks_);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type marks(marks_);
    Rcpp::traits::input_parameter< const std::string& >::type distm(distm_);
    
    Eigen::MatrixXd Zmrks = rcppstandardize_rates(tiles,rates,xseed,yseed,marks,nmrks,distm);
    
    PROTECT(Zmrks_ = Rcpp::wrap(Zmrks));
  }
  UNPROTECT(1);
  return Zmrks_;
END_RCPP
}
RcppExport SEXP rEEMSplots__euclidean_dist( SEXP X_, SEXP Y_) {
BEGIN_RCPP
  SEXP D_;
  {
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X(X_);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Y(Y_);

    Eigen::MatrixXd D = euclidean_dist(X,Y);

    PROTECT(D_ = Rcpp::wrap(D));
  }
  UNPROTECT(1);
  return D_;
END_RCPP    
}
