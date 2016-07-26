
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
    Eigen::MatrixXd Hinv = - M; Hinv.diagonal() += M.rowwise().sum(); Hinv.array() += 1.0;
    Eigen::MatrixXd H = Hinv.inverse();
    Eigen::MatrixXd R = - 2.0 * H;
    Eigen::VectorXd u = Eigen::VectorXd::Ones(M.rows());
    R.noalias() += H.diagonal() * u.transpose();
    R.noalias() += u * H.diagonal().transpose();
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
Eigen::VectorXd compute_contour_vals(const Eigen::MatrixXd &marks,
                                     const Eigen::VectorXd &now_rates, const Eigen::MatrixXd &now_seeds,
                                     const std::string &distm) {
    int nmrks = marks.rows();
    Eigen::VectorXd zvals = Eigen::VectorXd::Zero(nmrks);
    Eigen::MatrixXd dists = pairwise_dist(marks, now_seeds, distm);
    for ( int row = 0, closest = 0 ; row < nmrks ; row++ ) {
        dists.row(row).minCoeff( &closest );
        zvals(row) = now_rates( closest );
    }
    return zvals;
}
// Compute the average contour, by calling compute_contour_vals repeatedly
//
// [[Rcpp::export]]
void rcppstandardize_rates(const Eigen::VectorXd &tiles, const Eigen::VectorXd &rates,
                           const Eigen::MatrixXd &seeds, const Eigen::MatrixXd &marks,
                           const std::string &distm,
                           Eigen::VectorXd &zvals, Eigen::VectorXd &prgt0, Eigen::VectorXd &prlt0) {
    bool use_weighted_mean = true;
    int nmrks = marks.rows();
    Eigen::VectorXd zval = Eigen::VectorXd::Zero(nmrks);
    Eigen::VectorXd ones = Eigen::VectorXd::Ones(nmrks);
    // Vector to store the migration rates m (or diversity rates q).
    zvals = Eigen::VectorXd::Zero(nmrks);
    // Vectors to store the number of times (m > 0) and (m < 0).
    prgt0 = Eigen::VectorXd::Zero(nmrks);
    prlt0 = Eigen::VectorXd::Zero(nmrks);
    for ( int iter = 0, start = 0 ; iter < tiles.size() ; iter++ ) {
        int now_tiles = (int)tiles(iter);
        Eigen::VectorXd now_rates = rates.segment(start, now_tiles);
        Eigen::MatrixXd now_seeds = seeds.block(start, 0, now_tiles, 2);
        if (use_weighted_mean) {
            zval = compute_contour_vals(marks, now_rates, now_seeds, distm);
            zval = zval.array() - zval.mean();
        } else {
            now_rates = now_rates.array() - now_rates.mean();
            zval = compute_contour_vals(marks, now_rates, now_seeds, distm);
        }
        prgt0 += (zval.array() > 0.0).select(ones, 0.0);
        prlt0 += (zval.array() < 0.0).select(ones, 0.0);
        zvals += zval;
        start += now_tiles;
    }
    // Do not divide by the number of iterations here but in 'average.eems.contours' instead
    // as we might be combing results from several output directories.
    // zvals = zvals.array() / tiles.size();
}
// Compute the average contour, by calling compute_contour_vals repeatedly,
// this time without normalizing the rates to have mean 0. Without normalizing,
// it doesn't make sense to keep track of the number of times (m > 0) and (m < 0).
//
// [[Rcpp::export]]
void rcppdont_standardize_rates(const Eigen::VectorXd &tiles, const Eigen::VectorXd &rates,
                                const Eigen::MatrixXd &seeds, const Eigen::MatrixXd &marks,
                                const std::string &distm,
                                Eigen::VectorXd &zvals) {
    zvals = Eigen::VectorXd::Zero(marks.rows());
    for ( int iter = 0, start = 0 ; iter < tiles.size() ; iter++ ) {
        int now_tiles = (int)tiles(iter);
        Eigen::VectorXd now_rates = rates.segment(start, now_tiles);
        Eigen::MatrixXd now_seeds = seeds.block(start, 0, now_tiles, 2);
        zvals += compute_contour_vals(marks, now_rates, now_seeds, distm);
        start += now_tiles;
    }
    // Do not divide by the number of iterations here but in 'average.eems.contours' instead
    // as we might be combing results from several output directories.
    // zvals = zvals.array() / tiles.size();
}
