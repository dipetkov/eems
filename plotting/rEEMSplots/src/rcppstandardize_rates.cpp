
// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

const double pi_180 = M_PI / 180.0;
const double Earth_radiusX2 = 2.0 * 6378137.0;

// Compute the effective resistance matrix
//
// [[Rcpp::export]]
Eigen::MatrixXd resistance_distance(const Eigen::MatrixXd &M) {
  int d = M.rows();
  Eigen::MatrixXd Hinv = - M; Hinv.diagonal() += M.rowwise().sum(); Hinv.array() += 1.0;
  Eigen::MatrixXd H = Hinv.inverse();
  Eigen::MatrixXd R = - 2.0 * H;
  Eigen::VectorXd u_d = Eigen::VectorXd::Ones(d);
  R.noalias() += H.diagonal() * u_d.transpose();
  R.noalias() += u_d * H.diagonal().transpose();
  return R;
}
// Squared Euclidean distance: Taking the square root is not necessary
// because EEMS uses the distances to find closest points. For example,
//   pairwise_distance(X,Y).col(0).minCoeff( &closest )
// finds the row/point in X that is closest to the first row/point in Y
Eigen::MatrixXd euclidean_dist(const Eigen::MatrixXd &X, const Eigen::MatrixXd &Y) {
  return(  X.rowwise().squaredNorm().eval().replicate(1,Y.rows())
	 + Y.rowwise().squaredNorm().eval().transpose().replicate(X.rows(),1)
	 - 2.0 * X * Y.transpose() );
}
// Great circle distance, using the haversine formula
Eigen::MatrixXd greatcirc_dist(const Eigen::MatrixXd &X, const Eigen::MatrixXd &Y) {
  int nr = X.rows();
  int nc = Y.rows();
  // Convert from degrees to radians
  Eigen::ArrayXXd lon1 = X.col(0).replicate(1, nc).array() * pi_180;
  Eigen::ArrayXXd lat1 = X.col(1).replicate(1, nc).array() * pi_180;
  Eigen::ArrayXXd lon2 = Y.col(0).transpose().replicate(nr, 1).array() * pi_180;
  Eigen::ArrayXXd lat2 = Y.col(1).transpose().replicate(nr, 1).array() * pi_180;
  // The haversine function is hav(theta) = (1 - cos(theta)) / 2
  Eigen::ArrayXXd hav_lon = 0.5 * (1.0 - (lon2 - lon1).cos());
  Eigen::ArrayXXd hav_lat = 0.5 * (1.0 - (lat2 - lat1).cos());
  Eigen::ArrayXXd h = hav_lat + hav_lon * lat1.cos() * lat2.cos();
  // The great circle distance d is given by d = 2 * radius * arcsine( sqrt(h) )
  Eigen::MatrixXd d = (h < 1.0).select(h.sqrt(), 1.0).asin(); // Avoid numerical issues by ensuring h <= 1.0
  return (Earth_radiusX2 * d);
}
// Compute pairwise distances between the rows of X and the rows of Y.
// Choose either Euclidean or great circle distance
Eigen::MatrixXd pairwise_dist(const Eigen::MatrixXd &X, const Eigen::MatrixXd &Y, const std::string &distm) {
  if (!distm.compare("greatcirc")) {
    return (greatcirc_dist(X, Y));
  } else {
    return (euclidean_dist(X, Y));
  }
}
// Compute one contour, by filling in each of the pixels/marks
void compute_contour_vals(Eigen::MatrixXd &zvals, const Eigen::MatrixXd &marks,
			  const Eigen::VectorXd &now_rates, const Eigen::MatrixXd &now_seeds,
			  const std::string &distm) {
  int nxmrks = zvals.cols();
  int nymrks = zvals.rows();
  Eigen::MatrixXd distances = pairwise_dist(marks, now_seeds, distm);
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
// [[Rcpp::export]]
Eigen::MatrixXd rcppstandardize_rates(const Eigen::VectorXd &tiles, const Eigen::VectorXd &rates,
                                      const Eigen::VectorXd &xseed, const Eigen::VectorXd &yseed,
                                      const Eigen::MatrixXd &marks, const Eigen::VectorXd &nmrks,
                                      const std::string &distm) {
  bool use_weighted_mean = true;
  int nxmrks = nmrks(0);
  int nymrks = nmrks(1);
  Eigen::MatrixXd Zvals = Eigen::MatrixXd::Zero(nymrks, nxmrks);
  Eigen::MatrixXd zvals = Eigen::MatrixXd::Zero(nymrks, nxmrks);
  int niters = tiles.size();
  for ( int i = 0, pos = 0 ; i < niters ; i++ ) {
    int now_tiles = (int)tiles(i);
    Eigen::VectorXd now_rates = rates.segment(pos, now_tiles);
    Eigen::VectorXd now_xseed = xseed.segment(pos, now_tiles);
    Eigen::VectorXd now_yseed = yseed.segment(pos, now_tiles);
    Eigen::MatrixXd now_seeds(now_tiles, 2);
    now_seeds << now_xseed, now_yseed;
    if (use_weighted_mean) {
      compute_contour_vals(zvals, marks, now_rates, now_seeds, distm);
      zvals = zvals.array() - zvals.mean();
    } else {
      now_rates = now_rates.array() - now_rates.mean();
      compute_contour_vals(zvals, marks, now_rates, now_seeds, distm);
    }
    Zvals += zvals; pos += now_tiles;
  }
  // Do not divide by niters here but in 'average.eems.contours' instead
  // Zvals = Zvals.array() / niters;
  return Zvals.transpose();
}
// Compute the average contour, by calling compute_contour_vals repeatedly
//
// [[Rcpp::export]]
Eigen::MatrixXd rcppdont_standardize_rates(const Eigen::VectorXd &tiles, const Eigen::VectorXd &rates,
					   const Eigen::VectorXd &xseed, const Eigen::VectorXd &yseed,
					   const Eigen::MatrixXd &marks, const Eigen::VectorXd &nmrks,
					   const std::string &distm) {
  int nxmrks = nmrks(0);
  int nymrks = nmrks(1);
  Eigen::MatrixXd Zvals = Eigen::MatrixXd::Zero(nymrks, nxmrks);
  Eigen::MatrixXd zvals = Eigen::MatrixXd::Zero(nymrks, nxmrks);
  int niters = tiles.size();
  for ( int i = 0, pos = 0 ; i < niters ; i++ ) {
    int now_tiles = (int)tiles(i);
    Eigen::VectorXd now_rates = rates.segment(pos, now_tiles);
    Eigen::VectorXd now_xseed = xseed.segment(pos, now_tiles);
    Eigen::VectorXd now_yseed = yseed.segment(pos, now_tiles);
    Eigen::MatrixXd now_seeds(now_tiles, 2);
    now_seeds << now_xseed, now_yseed;
    compute_contour_vals(zvals, marks, now_rates, now_seeds, distm);
    Zvals += zvals; pos += now_tiles;
  }
  // Do not divide by niters here but in 'average.eems.contours' instead
  // Zvals = Zvals.array() / niters;
  return Zvals.transpose();
}
