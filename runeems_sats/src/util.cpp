
#include "util.hpp"

extern string dist_metric;

Params::Params( ) {
  seed = 1;
  diploid = true;
  datapath = "";
  mcmcpath = "";
  prevpath = "";
  gridpath = "";
  nDemes = 0;
  nIndiv = 0;
  nSites = 0;
  numMCMCIter = 0;
  numBurnIter = 0;
  numThinIter = 0;
  mSeedsProposalS2 = 0.01;
  qSeedsProposalS2 = 0.1;
  mEffctProposalS2 = 0.1;
  qEffctProposalS2 = 0.001;
  mrateMuProposalS2 = 0.01;
  qVoronoiPr = 0.25;
  qrateShape_2 = 0.001;
  mrateShape_2 = 0.001;
  sigmaShape_2 = 0.001;
  mrateScale_2 = 1.0;
  qrateScale_2 = 1.0;
  sigmaScale_2 = 1.0;
  negBiProb = 0.67;
  negBiSize = 10;
  mEffctHalfInterval = 2.0;
  qEffctHalfInterval = 0.1;
  mrateMuHalfInterval = 2.4771; // log10(300)
  /*
    There are two functions for testing whether the prior and the likelihood are computed correctly:
    * test_prior(parameters) + test_likelihood(parameters)
    * eval_prior() + eval_likelihood(), which use the current parameter values
    If testing = true, then the function check_ll_computation() will be called after each MCMC iteration,
    and it will check that the current prior, nowpi, is the same as test_prior(current parameters)
    as well as that the currrent likelihood, nowll, is the same as test_likelihood(current parameters)
  */  
  testing = false;
}
Params::~Params( ) { }

void Params::check_input_arguments() {
  check_condition(nDemes > 1, "Check that nDemes > 1 (or type `./runeems_snps --help`)");
  check_condition(nIndiv > 2, "Check that nIndiv > 2 (or type `./runeems_snps --help`)");
  check_condition(nSites > 0, "Check that nSites > 0");
  check_condition(numMCMCIter > 10, "Check that numMCMCIter > 10. But numMCMCIter should be >> 10.");
  check_condition(numBurnIter >= 0, "Check that numBurnIter >= 0");
  check_condition(numThinIter >= 0, "Check that numThinIter >= 0");
  check_condition(numMCMCIter > numBurnIter + numThinIter,
		  "Check that numMCMCIter > numBurnIter + numThinIter");
  check_condition(qVoronoiPr > 0.0 && qVoronoiPr < 0.5,
		  "Check that qVoronoiPr in (0,0.5).\n\tIt is a good idea to update migration more often than diversity.");
  check_condition(mSeedsProposalS2 > 0.0, "Check that mSeedsProposalS2 > 0");
  check_condition(qSeedsProposalS2 > 0.0, "Check that qSeedsProposalS2 > 0");
  check_condition(mEffctProposalS2 > 0.0, "Check that mEffctProposalS2 > 0");
  check_condition(qEffctProposalS2 > 0.0, "Check that qEffctProposalS2 > 0");
  check_condition(mrateMuProposalS2 > 0.0, "Check that mrateMuProposalS2 > 0");
  check_condition(mrateShape_2 > 0.0, "Check that mrateShape > 0");
  check_condition(qrateShape_2 > 0.0, "Check that qrateShape > 0");
  check_condition(sigmaShape_2 > 0.0, "Check that sigmaShape > 0");
  check_condition(mrateScale_2 > 0.0, "Check that mrateScale > 0");
  check_condition(qrateScale_2 > 0.0, "Check that qrateScale > 0");
  check_condition(sigmaScale_2 > 0.0, "Check that sigmaScale > 0");
  check_condition(negBiSize > 0, "Check that negBiSize > 0");
  check_condition(negBiProb > 0.0 && negBiProb < 1.0,
		  "Check that negBiProb in (0,1)");
  check_condition(dist_metric.compare("euclidean") || dist_metric.compare("greatcirc"),
		  "Check that 'distance' is either 'euclidean' or 'greatcirc'");
  /////////////////////////////////////////////////////
  mrateShape_2 /= 2.0;
  qrateScale_2 /= 2.0;
  mrateShape_2 /= 2.0;
  qrateScale_2 /= 2.0;
  sigmaShape_2 /= 2.0;
  sigmaScale_2 /= 2.0;
  /////////////////////////////////////////////////////
  cout << "Using Boost " << BOOST_LIB_VERSION << " (EEMS was tested with version 1_57)" << endl
       << "      Eigen " << EIGEN_WORLD_VERSION << "." << EIGEN_MAJOR_VERSION << "." << EIGEN_MINOR_VERSION
       << " (EEMS was tested with version 3.2.2)" << endl << endl;
  check_condition(numeric_limits<double>::has_infinity, "No representation for infinity");
  check_condition(boost::filesystem::exists(datapath + ".coord"),
		  "Check that " + datapath + ".coord exists");
  check_condition(boost::filesystem::exists(datapath + ".sites"),
		  "Check that " + datapath + ".sites exists");
  check_condition(boost::filesystem::exists(datapath + ".outer"),
		  "Check that " + datapath + ".outer exists");
  if ( !gridpath.empty() ) {
    check_condition(boost::filesystem::exists(gridpath + ".demes"),
		    "Check that " + gridpath + ".demes exists");
    check_condition(boost::filesystem::exists(gridpath + ".edges"),
		    "Check that " + gridpath + ".edges exists");
  }
  boost::filesystem::path mcmcdir(mcmcpath.c_str());
  boost::filesystem::path prevdir(prevpath.c_str());
  if (exists(mcmcdir)) {
    cout << "  Will save EEMS output to " << mcmcpath << "/" << endl;
  } else {
    check_condition(boost::filesystem::create_directory(mcmcdir),
		    "Check that can create directory " + mcmcpath);
  }
  if (exists(prevdir)) {
    cout << "  Will load EEMS output from " << prevpath << "/" << endl;
  } else {
    prevpath.clear();
  }
}
ostream& operator<<(ostream& out, const Params& params) {
  out << "               datapath = " << params.datapath << endl
      << "               mcmcpath = " << params.mcmcpath << endl
      << "               prevpath = " << params.prevpath << endl
      << "               gridpath = " << params.gridpath << endl    
      << "               distance = " << dist_metric << endl
      << "                diploid = " << params.diploid << endl
      << "                 nIndiv = " << params.nIndiv << endl
      << "                 nSites = " << params.nSites << endl
      << "                 nDemes = " << params.nDemes << endl
      << "                   seed = " << params.seed << endl
      << "            numMCMCIter = " << params.numMCMCIter << endl
      << "            numBurnIter = " << params.numBurnIter << endl
      << "            numThinIter = " << params.numThinIter << endl
      << "              negBiSize = " << params.negBiSize << endl
      << fixed << setprecision(6)
      << "              negBiProb = " << params.negBiProb << endl
      << "             qVoronoiPr = " << params.qVoronoiPr << endl
      << "             mrateShape = " << 2.0*params.mrateShape_2 << endl
      << "             qrateShape = " << 2.0*params.qrateShape_2 << endl
      << "             sigmaShape = " << 2.0*params.sigmaShape_2 << endl
      << "             qrateScale = " << 2.0*params.qrateScale_2 << endl
      << "             mrateScale = " << 2.0*params.mrateScale_2 << endl
      << "             sigmaScale = " << 2.0*params.sigmaScale_2 << endl
      << "       mSeedsProposalS2 = " << params.mSeedsProposalS2 << endl
      << "       qSeedsProposalS2 = " << params.qSeedsProposalS2 << endl
      << "       mEffctProposalS2 = " << params.mEffctProposalS2 << endl
      << "       qEffctProposalS2 = " << params.qEffctProposalS2 << endl
      << "      mrateMuProposalS2 = " << params.mrateMuProposalS2 << endl;
  return out;
}
VectorXd split(const string &line) {
  istringstream in(line);
  vector<double> numbers;
  double number;
  while (!in.eof()) { in >> number;
    if (!in.fail()) { numbers.push_back(number); } else { break; }
  }
  if (in.fail()||in.bad()) { return (VectorXd::Zero(0)); }
  // Since we don't know how many numbers there are,
  // first we store the data in a std::vector of doubles,
  // and then typecast it to Eigen::VectorXd
  return (VectorXd::Map(&numbers[0],numbers.size()));
}
bool is_finite(const double x) {
  return ((x != Inf) && (x != -Inf));
}
bool isposdef(const MatrixXd &A) {
  check_condition(A.rows() == A.cols(), "isposdef : matrix A is not square");
  SelfAdjointEigenSolver<MatrixXd> eig(A, EigenvaluesOnly);
  double smallest_eigen_val = eig.eigenvalues().minCoeff();
  return (smallest_eigen_val > 0);
}
bool isdistmat(const MatrixXd &A) {
  check_condition(A.rows() == A.cols(), "isdistmat : matrix A is not square");
  double eps = 1e-12;
  bool A_is_nonnegative = A.minCoeff() >= 0;
  bool diag_is_zeros = (A.diagonal().minCoeff() == 0) && (A.diagonal().maxCoeff() == 0);
  bool A_is_symmetric = abs((A - A.transpose()).maxCoeff()) < eps;
  SelfAdjointEigenSolver<MatrixXd> eig(A, EigenvaluesOnly);
  ArrayXd eigen_vals = eig.eigenvalues().array();
  int num_neg_eigen_vals = (eigen_vals < -eps).count();
  int num_pos_eigen_vals = (eigen_vals > +eps).count();
  int num_zero_eigen_vals = A.rows() - num_neg_eigen_vals - num_pos_eigen_vals;
  return (A_is_nonnegative && diag_is_zeros && A_is_symmetric &&
	  (num_pos_eigen_vals == 1) && (num_zero_eigen_vals == 0));
}
double logdet(const MatrixXd &A) {
  return (A.selfadjointView<Lower>().ldlt().vectorD().array().log().sum());
}
double pseudologdet(const MatrixXd &A, const int rank) {
  SelfAdjointEigenSolver<MatrixXd> x(A);
  return (x.eigenvalues().reverse().array().head(rank).log().sum());
}
double wishpdfln(const MatrixXd &X, const MatrixXd &Sigma, const double df) {
  double ldX = logdet(X);
  double ldS = logdet(Sigma);
  int n = X.rows();
  double ln_p = 0.5 * (- df * ldS - Sigma.selfadjointView<Lower>().llt().solve(X).trace()) + 
    0.5 * ((df - n - 1.0) * ldX - df * n * ln_2) - mvgammaln(0.5 * df, n);
  return ln_p;
}
// Not a general pseudowishpdfln (because L*Diff*L' at a single locus has rank 1
double pseudowishpdfln(const MatrixXd &X, const MatrixXd &Sigma, const int df) {
  int rank = 1;
  double ldX = pseudologdet(X, rank);
  double ldS = logdet(Sigma);
  int n = X.rows();
  int q = (df < n) ? df : n;
  double ln_p = 0.5 * (- df * ldS - Sigma.selfadjointView<Lower>().llt().solve(X).trace()) +
    0.5 * ((df - n - 1.0) * ldX - df * n * ln_2 - df * (n - q) * ln_pi) - mvgammaln(0.5 * df, q);
  return ln_p;
}
double mvgammaln(const double a, const int p) {
  double val = ln_pi * p * (p-1) / 4.0;
  for (int i = 0 ; i < p ; i++ ) {
    val += lgamma(a - 0.5 * i);
  }
  return (val);
}
/*
  Log probability mass function of the negative binomial distribution, with
  parameters size (number of failures) and prob (probability of success)
  * up to a constant of proportionality that depends on size and prob,
  but not on the random variable k (number of successes)
  Also -- log pmf of the zero-truncated negative binomial, with the same parameters,
  as long as k>0, since truncating only changes the constant of proportionality
  Uses gamma(n) = (n-1)!
*/
double dnegbinln(const int k, const int size, const double prob) {
  double ln_p = -Inf;
  if ( k >= 0 && ( size > 0 && prob > 0.0 && prob < 1.0 ) ) {
    ln_p = lgamma(size + k) - lgamma(k + 1) + k * log(prob);
  }
  return ln_p;
}
/*
  Log probability density function of the inverse gamma distribution,
  with parameters shape and scale
  * Up to a constant of proportionality that depends on the shape and the scale,
  but not on the random variable x
*/
double dinvgamln(const double x, const double shape, const double scale) {
  double ln_p = -Inf;
  if ( x > 0.0 && ( shape > 0.0 && scale > 0.0 ) ) {
    ln_p = - (shape + 1.0) * log(x) - scale / x;
  }
  return ln_p;
}
/*
  Truncated normal probability density function with support [ - halfInterval, + halfInterval],
  on the log scale, including the normalizing constant
*/
double dtrnormln(const double x, const double mu, const double sigma2, const double halfInterval) {
  double ln_p = -Inf;
  if ( ( x >= -halfInterval && x <= halfInterval ) && sigma2 > 0.0 ) {
    boost::math::normal pnorm(mu,sqrt(sigma2));
    ln_p = - 0.5 * log(sigma2) - 0.5 * (x-mu) * (x-mu) / sigma2
      - log(cdf(pnorm,halfInterval) - cdf(pnorm,-halfInterval)); // the normalizing constant
  }
  return ln_p;
}
/*
  Squared Euclidean distance: Taking the square root is not necessary
  because EEMS uses the distances to find closest points. For example,
  pairwise_distance(X,Y).col(0).minCoeff( &closest )
  finds the row/point in X that is closest to the first row/point in Y
*/
MatrixXd euclidean_dist(const MatrixXd &X, const MatrixXd &Y) {
  check_condition(X.cols() == Y.cols(), "euclidean_dist : ncol(X) != ncol(Y)");
  return (  X.rowwise().squaredNorm().eval().replicate(1, Y.rows())
	    + Y.rowwise().squaredNorm().eval().transpose().replicate(X.rows(), 1)
	    - 2.0 * X * Y.transpose() );
}
// Great circle distance, up to a constant of proportionality equal to 2*R
// where R is the earth's radius
MatrixXd greatcirc_dist(const MatrixXd &X, const MatrixXd &Y) {
  check_condition(X.cols() == 2 && Y.cols() == 2,
		  "greatcirc_dist : ncol(X) != 2 or ncol(Y) != 2");
  int nr = X.rows();
  int nc = Y.rows();
  // Convert from degrees to radians
  double degree = boost::math::double_constants::degree;
  ArrayXXd lon1 = X.col(0).replicate(1, nc).array() * degree;
  ArrayXXd lat1 = X.col(1).replicate(1, nc).array() * degree;
  ArrayXXd lon2 = Y.col(0).transpose().replicate(nr, 1).array() * degree;
  ArrayXXd lat2 = Y.col(1).transpose().replicate(nr, 1).array() * degree;
  // The haversine function is hav(theta) = (1 - cos(theta)) / 2
  ArrayXXd hav_lon = 0.5 * (1.0 - (lon2 - lon1).cos());
  ArrayXXd hav_lat = 0.5 * (1.0 - (lat2 - lat1).cos());
  ArrayXXd h = hav_lat + hav_lon * lat1.cos() * lat2.cos();
  // The great circle distance d is given by d = 2 * radius * arcsine( sqrt(h) )
  MatrixXd d = (h < 1.0).select(h.sqrt(), 1.0).asin(); // Avoid numerical issues by ensuring h <= 1.0
  return (Earth_radiusX2 * d);
}
// Compute pairwise distances between the rows of X and the rows of Y.
// Choose either Euclidean or great circle distance
MatrixXd pairwise_distance(const MatrixXd &X, const MatrixXd &Y) {
  if (!dist_metric.compare("greatcirc")) {
    return (greatcirc_dist(X,Y));
  } else {
    return (euclidean_dist(X,Y));
  }
}
MatrixXd resistance_distance(const MatrixXd &M) {
  check_condition(M.cols() == M.rows(), "resistance distance : M is not symmetric");
  check_condition(M.minCoeff() >= 0, "resistance distance : M has negative weights");
  MatrixXd Hinv = - M; Hinv.diagonal() += M.rowwise().sum();
  Hinv.array() += 1.0;
  MatrixXd H = Hinv.inverse();
  MatrixXd R = - 2.0 * H;
  VectorXd u = VectorXd::Ones(M.rows());
  R.noalias() += H.diagonal() * u.transpose();
  R.noalias() += u * H.diagonal().transpose();
  return R;
}
// Implements Equation S12 in the Supplementary 
MatrixXd expected_dissimilarities(const MatrixXd &J, const MatrixXd &M, const VectorXd &w) {
  int n = J.rows();
  int o = J.cols();
  VectorXd Jw = J * w.head(o);
  MatrixXd R = resistance_distance(M).topLeftCorner(o,o);
  VectorXd u = VectorXd::Ones(n);
  MatrixXd Delta = J * R * J.transpose();
  Delta.noalias() += 0.5 * Jw * u.transpose();
  Delta.noalias() += 0.5 * u * Jw.transpose();
  Delta.diagonal() -= Jw;
  return Delta;
}
// Read a matrix, with unknown dimensions, from a text file
// Return an empty matrix (0 rows, 0 columns) if there is an error
MatrixXd readMatrixXd(const string &filename) {
  ifstream instrm(filename.c_str(), ios::in);
  if (!instrm.is_open( )) { return (MatrixXd::Zero(0,0)); }
  string line; getline(instrm,line);
  boost::algorithm::trim(line);
  // Split the first line into numbers
  VectorXd row = split(line);
  // This tells us the number of columns
  int cols = row.size();
  if (!cols) { return (MatrixXd::Zero(0,0)); }
  MatrixXd mat(1,cols); mat.row(0) = row;
  while (getline(instrm,line)) {
    boost::algorithm::trim(line);
    row = split(line);
    if (row.size()==cols) {
      // Resize the matrix to add each row (done once at the start)
      insertRow(mat,row);
    } else {
      return (MatrixXd::Zero(0,0));
    }
  }
  instrm.close();
  return mat;
}
// If there are two draws and the first has two tiles and the second -- three tiles,
// then sizes = c(2,3) and array = c(m_{1t_1},m_{1t_2},m_{2t_1},m_{2t_2},m_{2t_3})
void dlmcell(const string &filename, const VectorXd &sizes, const vector<double> &array) {
  check_condition(array.size() == sizes.sum(), "dlmcell : length(array) != sum(sizes).");
  ofstream out(filename.c_str(),ofstream::out);
  check_condition(out.is_open(), "Cannot open " + filename + " for writing");
  vector<double>::const_iterator it = array.begin();
  for ( int i = 0 ; i < sizes.size() ; i++ ) {
    for ( int j = 0 ; j < sizes(i) ; j++ ) {
      out << fixed << setprecision(6) << *it << " "; it++;
      check_condition(!out.bad(), "Failed writing to " + filename);
    }
    out << endl;
  }
  out.close( );
}
void dlmwrite(const string &filename, const MatrixXd & mat) {
  ofstream out(filename.c_str(),ofstream::out);
  out << fixed << setprecision(6) << mat << endl;
  check_condition(!out.bad(), "Failed writing to " + filename);
  out.close( );
}
void removeElem(VectorXd &vec, const int elemToRemove)
{
  int elems = vec.size() - 1;
  if (elemToRemove <= elems) {
    // Copy the last element into the element to remove
    if (elemToRemove < elems) {
      vec(elemToRemove) = vec(elems);
    }
    vec.conservativeResize(elems);
  }
}
void removeRow(MatrixXd &mat, const int rowToRemove)
{
  int rows = mat.rows() - 1;
  int cols = mat.cols();
  if( rowToRemove <= rows ) {
    // Copy the last row into the row to remove
    if ( rowToRemove < rows ) {
      mat.row(rowToRemove) = mat.row(rows);
    }
    mat.conservativeResize(rows,cols);
  }
}
void insertElem(VectorXd &vec, const double &elem)
{
  int elems = vec.size();
  vec.noalias() = (VectorXd(elems+1) << vec, elem).finished();
}
void insertRow(MatrixXd &mat, const VectorXd &row)
{
  check_condition(mat.cols() == row.size(), "insertRow : ncol(mat) != size(row)");
  int rows = mat.rows();
  int cols = mat.cols();
  mat.noalias() = (MatrixXd(rows+1,cols) << mat, row.transpose()).finished();
}
VectorXd slice(const VectorXd &A, const VectorXi &I) {
  check_condition(I.minCoeff() >= 0 && I.maxCoeff() < A.size(),
		  "Check that I is an index vector with elements in the range [0, size - 1].");
  int size = I.size();
  VectorXd B(size);
  for ( int i = 0 ; i < size ; i++ ) {
    B(i) = A(I(i));
  }
  return B;
}
MatrixXd slice(const MatrixXd &A, const VectorXi &R, const VectorXi &C) {
  check_condition(R.minCoeff() >= 0 && R.maxCoeff() < A.rows(),
		  "Check that R is an index vector with elements in the range [0, nrow - 1).");
  check_condition(C.minCoeff() >= 0 && C.maxCoeff() < A.cols(),
		  "Check that C is an index vector with elements in the range [0, ncol - 1].");
  int rows = R.size();
  int cols = C.size();
  MatrixXd B(rows,cols);
  for ( int i = 0 ; i < rows ; i++ ) {
    for ( int j = 0 ; j < cols ; j++ ) { 
      B(i,j) = A(R(i), C(j));
    }
  }
  return B;
}
void check_condition(const bool condition, const string &message) {
  if (condition == false) {
    cerr << "Error: " << message << endl;
    exit(1);
  }
}
