
// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

const double pi_180 = M_PI / 180.0;
const double Earth_radiusX2 = 2.0 * 6378137.0;

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
Eigen::VectorXd compute_contours(const Eigen::MatrixXd &marks, const Eigen::VectorXd &now_rates,
				 const Eigen::MatrixXd &now_seeds, const std::string &distm) {
  int nmrks = marks.rows();
  Eigen::VectorXd raster = Eigen::VectorXd::Zero(nmrks);
  Eigen::MatrixXd dists = pairwise_dist(marks, now_seeds, distm);
  for ( int row = 0, closest = 0 ; row < nmrks ; row++ ) {
    dists.row(row).minCoeff( &closest );
    raster(row) = now_rates( closest );
  }
  return raster;
}
// Compute the average contour, by calling compute_contours repeatedly
//
// [[Rcpp::export]]
Rcpp::List tiles2contours_standardize(const Eigen::VectorXd &tiles,
                                      const Eigen::VectorXd &rates,
                                      const Eigen::MatrixXd &seeds,
                                      const Eigen::MatrixXd &marks,
                                      const std::string &distm) {
  bool use_weighted_mean = true;
  int nmrks = marks.rows();
  Eigen::VectorXd zval = Eigen::VectorXd::Zero(nmrks);
  Eigen::VectorXd ones = Eigen::VectorXd::Ones(nmrks);
  // Vector to store the migration rates m (or diversity rates q).
  Eigen::VectorXd raster = Eigen::VectorXd::Zero(nmrks);
  // Vectors to store the number of times (m > 0) and (m < 0).
  Eigen::VectorXd pr_gt0 = Eigen::VectorXd::Zero(nmrks);
  Eigen::VectorXd pr_lt0 = Eigen::VectorXd::Zero(nmrks);
  for ( int iter = 0, pos = 0 ; iter < tiles.size() ; iter++ ) {
    int now_tiles = (int)tiles(iter);
    Eigen::VectorXd now_rates = rates.segment(pos, now_tiles);
    Eigen::MatrixXd now_seeds = seeds.block(pos, 0, now_tiles, 2);
    if (use_weighted_mean) {
      zval = compute_contours(marks, now_rates, now_seeds, distm);
      zval = zval.array() - zval.mean();
    } else {
      now_rates = now_rates.array() - now_rates.mean();
      zval = compute_contours(marks, now_rates, now_seeds, distm);
    }
    pr_gt0 += (zval.array() > 0.0).select(ones, 0.0);
    pr_lt0 += (zval.array() < 0.0).select(ones, 0.0);
    raster += zval;
    pos += now_tiles;
  }
  // Do not divide by the number of iterations here
  // as we might be combing results from several output directories.
  return Rcpp::List::create(Rcpp::Named("zvals") = raster,
                            Rcpp::Named("prgt0") = pr_gt0,
                            Rcpp::Named("prlt0") = pr_lt0);
}
// Compute the average contour, by calling compute_contours repeatedly,
// this time without normalizing the rates to have mean 0. Without normalizing,
// it doesn't make sense to keep track of the number of times (m > 0) and (m < 0).
//
// [[Rcpp::export]]
Rcpp::List tiles2contours(const Eigen::VectorXd &tiles,
                          const Eigen::VectorXd &rates,
                          const Eigen::MatrixXd &seeds,
                          const Eigen::MatrixXd &marks,
                          const std::string &distm) {
  Eigen::VectorXd raster = Eigen::VectorXd::Zero(marks.rows());
  for ( int iter = 0, pos = 0 ; iter < tiles.size() ; iter++ ) {
    int now_tiles = (int)tiles(iter);
    Eigen::VectorXd now_rates = rates.segment(pos, now_tiles);
    Eigen::MatrixXd now_seeds = seeds.block(pos, 0, now_tiles, 2);
    raster += compute_contours(marks, now_rates, now_seeds, distm);
    pos += now_tiles;
  }
  // Do not divide by the number of iterations here
  // as we might be combing results from several output directories.
  return Rcpp::List::create(Rcpp::Named("zvals") = raster);
}
