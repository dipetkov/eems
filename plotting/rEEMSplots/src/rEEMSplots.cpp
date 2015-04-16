
#include <RcppEigen.h>
#include <Rcpp.h>
using namespace Rcpp;

using std::string;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;
using Eigen::VectorXi;

// rcppstandardize_rates
MatrixXd rcppstandardize_rates(const VectorXd &tilex, const VectorXd &rates,
			       const VectorXd &xseed, const VectorXd &yseed,
			       const MatrixXd &marks, const VectorXd &nmrks,
			       const std::string &distm);
MatrixXd rcppnotstandardize_rates(const VectorXd &tilex, const VectorXd &rates,
				  const VectorXd &xseed, const VectorXd &yseed,
				  const MatrixXd &marks, const VectorXd &nmrks,
				  const std::string &distm);
RcppExport SEXP rEEMSplots__rcppstandardize_rates( SEXP tiles_, SEXP rates_, SEXP xseed_, SEXP yseed_,
						   SEXP marks_, SEXP nmrks_, SEXP distm_) {
BEGIN_RCPP
  SEXP Zmrks_;
  {
    Rcpp::traits::input_parameter< const VectorXd& >::type tiles(tiles_);
    Rcpp::traits::input_parameter< const VectorXd& >::type rates(rates_);
    Rcpp::traits::input_parameter< const VectorXd& >::type xseed(xseed_);
    Rcpp::traits::input_parameter< const VectorXd& >::type yseed(yseed_);
    Rcpp::traits::input_parameter< const VectorXd& >::type nmrks(nmrks_);
    Rcpp::traits::input_parameter< const MatrixXd& >::type marks(marks_);
    Rcpp::traits::input_parameter< const string& >::type distm(distm_);
    
    Eigen::MatrixXd Zmrks = rcppstandardize_rates(tiles,rates,xseed,yseed,marks,nmrks,distm);
    
    PROTECT(Zmrks_ = Rcpp::wrap(Zmrks));
  }
  UNPROTECT(1);
  return Zmrks_;
END_RCPP
}
RcppExport SEXP rEEMSplots__rcppnotstandardize_rates( SEXP tiles_, SEXP rates_, SEXP xseed_, SEXP yseed_,
						      SEXP marks_, SEXP nmrks_, SEXP distm_) {
BEGIN_RCPP
  SEXP Zmrks_;
  {
    Rcpp::traits::input_parameter< const VectorXd& >::type tiles(tiles_);
    Rcpp::traits::input_parameter< const VectorXd& >::type rates(rates_);
    Rcpp::traits::input_parameter< const VectorXd& >::type xseed(xseed_);
    Rcpp::traits::input_parameter< const VectorXd& >::type yseed(yseed_);
    Rcpp::traits::input_parameter< const VectorXd& >::type nmrks(nmrks_);
    Rcpp::traits::input_parameter< const MatrixXd& >::type marks(marks_);
    Rcpp::traits::input_parameter< const string& >::type distm(distm_);
    
    Eigen::MatrixXd Zmrks = rcppnotstandardize_rates(tiles,rates,xseed,yseed,marks,nmrks,distm);
    
    PROTECT(Zmrks_ = Rcpp::wrap(Zmrks));
  }
  UNPROTECT(1);
  return Zmrks_;
END_RCPP
}
