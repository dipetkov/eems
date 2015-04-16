// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;
using Eigen::VectorXi;

const double pi_180 = M_PI/180.0;

// Squared Euclidean distance: Taking the square root is not necessary
// because EEMS uses the distances to find closest points. For example,
//   pairwise_distance(X,Y).col(0).minCoeff( &closest )
// finds the row/point in X that is closest to the first row/point in Y
//
// There is no need to make this function available from R
MatrixXd euclidean_dist(const MatrixXd &X, const MatrixXd &Y) {
  return(  X.rowwise().squaredNorm().eval().replicate(1,Y.rows())
	 + Y.rowwise().squaredNorm().eval().transpose().replicate(X.rows(),1)
	 - 2.0*X*Y.transpose() );
}
// Great circle distance, up to a constant of proportionality equal to 2*R
// where R is the earth's radius
//
// There is no need to make this function available from R
MatrixXd greatcirc_dist(const MatrixXd &X, const MatrixXd &Y) {
  int nr = X.rows();
  int nc = Y.rows();
  MatrixXd lon1 = X.col(0).replicate(1,nc) * pi_180;
  MatrixXd lat1 = X.col(1).replicate(1,nc) * pi_180;
  MatrixXd lon2 = Y.col(0).transpose().replicate(nr,1) * pi_180;
  MatrixXd lat2 = Y.col(1).transpose().replicate(nr,1) * pi_180;
  MatrixXd dlon = 0.5*(lon2 - lon1);
  MatrixXd dlat = 0.5*(lat2 - lat1);
  MatrixXd a = dlat.array().sin().square().matrix() +
    (dlon.array().sin().square() * lat1.array().cos() * lat2.array().cos()).matrix();
  MatrixXd c = (a.array()<1.0).select(a.array().sqrt(),MatrixXd::Ones(nr,nc)).array().asin();
  return (c); // Instead of (2*R*c) where R = 6378137 is the Earth's radius.
}
// Compute pairwise distances between the rows of X and the rows of Y.
// Choose either Euclidean or great circle distance
//
// There is no need to make this function available from R
MatrixXd pairwise_dist(const MatrixXd &X, const MatrixXd &Y, const std::string &distm) {
  if (!distm.compare("greatcirc")) {
    return (greatcirc_dist(X,Y));
  } else {
    return (euclidean_dist(X,Y));
  }
}
// Compute one contour, by filling in each of the pixels/marks
//
// There is no need to make this function available from R
void compute_contour_vals(MatrixXd &zvals, const MatrixXd &marks,
			  const VectorXd &now_rates, const MatrixXd &now_seeds,
			  const std::string &distm) {
  int nxmrks = zvals.cols();
  int nymrks = zvals.rows();
  MatrixXd distances = pairwise_dist(marks,now_seeds,distm);
  int closest = 0;
  for ( int j = 0 ; j < nxmrks*nymrks ; j++ ) {
    int r = j % nxmrks;
    int c = j / nxmrks;
    distances.row(j).minCoeff( &closest );
    zvals(c,r) = now_rates( closest );
  }
}
// Compute the average contour, by calling compute_contour_vals repeatedly
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::export]]
MatrixXd rcppstandardize_rates(const VectorXd &tiles, const VectorXd &rates,
                               const VectorXd &xseed, const VectorXd &yseed,
                               const MatrixXd &marks, const VectorXd &nmrks,
                               const std::string &distm) {
  bool use_weighted_mean = true;
  int nxmrks = nmrks(0);
  int nymrks = nmrks(1);
  MatrixXd Zvals = MatrixXd::Zero(nymrks,nxmrks);
  MatrixXd zvals = MatrixXd::Zero(nymrks,nxmrks);
  int niters = tiles.size();
  for ( int i = 0, pos = 0 ; i < niters ; i++ ) {
    int now_tiles = (int)tiles(i);
    VectorXd now_rates = rates.segment(pos,now_tiles);
    VectorXd now_xseed = xseed.segment(pos,now_tiles);
    VectorXd now_yseed = yseed.segment(pos,now_tiles);
    MatrixXd now_seeds(now_tiles, 2);
    now_seeds << now_xseed,now_yseed;
    if (use_weighted_mean) {
      compute_contour_vals(zvals,marks,now_rates,now_seeds,distm);
      zvals = zvals.array() - zvals.mean();
    } else {
      now_rates = now_rates.array() - now_rates.mean();
      compute_contour_vals(zvals,marks,now_rates,now_seeds,distm);
    }
    Zvals += zvals; pos += now_tiles;
  }
  ///Zvals = Zvals.array() / niters;
  return Zvals.transpose();
}
MatrixXd rcppnotstandardize_rates(const VectorXd &tiles, const VectorXd &rates,
                                  const VectorXd &xseed, const VectorXd &yseed,
                                  const MatrixXd &marks, const VectorXd &nmrks,
                                  const std::string &distm) {
  bool use_weighted_mean = true;
  int nxmrks = nmrks(0);
  int nymrks = nmrks(1);
  MatrixXd Zvals = MatrixXd::Zero(nymrks,nxmrks);
  MatrixXd zvals = MatrixXd::Zero(nymrks,nxmrks);
  int niters = tiles.size();
  for ( int i = 0, pos = 0 ; i < niters ; i++ ) {
    int now_tiles = (int)tiles(i);
    VectorXd now_rates = rates.segment(pos,now_tiles);
    VectorXd now_xseed = xseed.segment(pos,now_tiles);
    VectorXd now_yseed = yseed.segment(pos,now_tiles);
    MatrixXd now_seeds(now_tiles, 2);
    now_seeds << now_xseed,now_yseed;
    compute_contour_vals(zvals,marks,now_rates,now_seeds,distm);
    Zvals += zvals; pos += now_tiles;
  }
  ///Zvals = Zvals.array() / niters;
  return Zvals.transpose();
}
