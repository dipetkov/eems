
#include "util.hpp"

Params::Params( ) { }
Params::~Params( ) { }
Params::Params(const string &params_file) {
  po::options_description eems_options("Options");
  eems_options.add_options()
    ("datapath",po::value<string>(&datapath)->default_value(""), "Datapath")
    ("mcmcpath",po::value<string>(&mcmcpath)->default_value(""), "MCMCpath")
    ("diploid", po::value<bool>(&diploid)->default_value(true), "diploid")
    ("nDemes", po::value<int>(&nDemes)->default_value(0), "nDemes")
    ("nIndiv", po::value<int>(&nIndiv)->default_value(0), "nIndiv")
    ("nSites", po::value<int>(&nSites)->default_value(0), "nSites")
    ("numMCMCIter", po::value<int>(&numMCMCIter)->default_value(1), "numMCMCIter")
    ("numBurnIter", po::value<int>(&numBurnIter)->default_value(0), "numBurnIter")
    ("numThinIter", po::value<int>(&numThinIter)->default_value(0), "numThinIter")
    ("mSeedsProposalS2", po::value<double>(&mSeedsProposalS2)->default_value(0.01), "mSeedsProposalS2")
    ("qSeedsProposalS2", po::value<double>(&qSeedsProposalS2)->default_value(0.1), "qSeedsProposalS2")
    ("mEffctProposalS2", po::value<double>(&mEffctProposalS2)->default_value(0.1), "mEffctProposalS2")
    ("qEffctProposalS2", po::value<double>(&qEffctProposalS2)->default_value(0.001), "qEffctProposalS2")
    ("mrateMuProposalS2", po::value<double>(&mrateMuProposalS2)->default_value(0.01), "mrateMuProposalS2")
    ("mrateShape", po::value<double>(&mrateShape)->default_value(0.001), "mrateShape")
    ("qrateShape", po::value<double>(&qrateShape)->default_value(0.001), "qrateShape")
    ("s2locShape", po::value<double>(&s2locShape)->default_value(0.001), "s2locShape")
    ("qrateScale", po::value<double>(&qrateScale)->default_value(1.0), "qrateScale")
    ("mrateScale", po::value<double>(&mrateScale)->default_value(1.0), "mrateScale")
    ("s2locScale", po::value<double>(&s2locScale)->default_value(1.0), "s2locScale")
    ("negBiProb", po::value<double>(&negBiProb)->default_value(0.67), "negBiProb")
    ("negBiSize", po::value<int>(&negBiSize)->default_value(10), "negBiSize") ;
  ////////////////////////////
  // Constants:
  testing = false;
  mEffctHalfInterval = 2.4771; // log10(300)
  qEffctHalfInterval = 0.1;
  mrateMuHalfInterval = 2.0;
  ifstream instrm(params_file.c_str());
  po::variables_map vm;
  try {
    po::store(po::parse_config_file(instrm,eems_options,true),vm);
    po::notify(vm);
  } catch (exception& e) {
    cerr << "[EEMS::Params] Error parsing input parameters in " << params_file << ": " << e.what() << endl; 
  }
  instrm.close();
  ////////////////////////////
  // Required argument: datapath, mcmcpath, nIndiv, nSites, nDemes
  if (!((datapath.length()>0)&&(mcmcpath.length()>0))) {
    cerr << "[EEMS::Params] The params file must specify datapath and mcmcpath" << endl;
    exit(1);
  }
  if (!((nIndiv>0)&&(nSites>0)&&(nDemes>0))) {
    cerr << "[EEMS::Params] The params file must specify nIndiv, nSites, nDemes" << endl;
    exit(1);
  }
  dlmwrite(cerr);
}
void Params::dlmwrite(ostream &out) const {
  out << "Input parameters:" << endl
      << "               datapath = " << datapath << endl
      << "               mcmcpath = " << mcmcpath << endl
      << "                 nIndiv = " << nIndiv << endl
      << "                 nSites = " << nSites << endl
      << "                 nDemes = " << nDemes << endl
      << "            numMCMCIter = " << numMCMCIter << endl
      << "            numBurnIter = " << numBurnIter << endl
      << "            numThinIter = " << numThinIter << endl
      << "              negBiSize = " << negBiSize << endl
      << fixed << setprecision(6)
      << "              negBiProb = " << negBiProb << endl
      << "             mrateShape = " << mrateShape << endl
      << "             qrateShape = " << qrateShape << endl
      << "             s2locShape = " << s2locShape << endl
      << "             qrateScale = " << qrateScale << endl
      << "             mrateScale = " << mrateScale << endl
      << "             s2locScale = " << s2locScale << endl
      << "       mSeedsProposalS2 = " << mSeedsProposalS2 << endl
      << "       qSeedsProposalS2 = " << qSeedsProposalS2 << endl
      << "       mEffctProposalS2 = " << mEffctProposalS2 << endl
      << "       qEffctProposalS2 = " << qEffctProposalS2 << endl
      << "      mrateMuProposalS2 = " << mrateMuProposalS2 << endl << endl;
}
void get_boost_version(ostream& out) {
  out << BOOST_LIB_VERSION;
}
void get_eigen_version(ostream& out) {
  out << EIGEN_WORLD_VERSION << "." << EIGEN_MAJOR_VERSION << "." << EIGEN_MINOR_VERSION;
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
  return (0.5*(-df*ldS - Sigma.selfadjointView<Lower>().llt().solve(X).trace() + 
	       (df-n-1.0)*ldX - df*n*log_2) - mvgammaln(0.5*df,n));
}
// Not a general pseudowishpdfln (because L*Diff*L' at a single locus has rank 1
double pseudowishpdfln(const MatrixXd &X, const MatrixXd &Sigma, const int df) {
  int rank = 1;
  double ldX = pseudologdet(X,rank);
  double ldS = logdet(Sigma);
  int n = X.rows();
  int q = (df<n) ? df : n;
  return (0.5*(-df*ldS - Sigma.selfadjointView<Lower>().llt().solve(X).trace() +
	       (df-n-1.0)*ldX - df*n*log_2 - df*(n-q)*log_pi) - mvgammaln(0.5*df,q));
}
double mvgammaln(const double a, const int p) {
  double val = 0.25*log_pi*p*(p-1);
  for (int i = 0 ; i < p ; i++ ) {
    val += lgamma(a-0.5*i);
  }
  return (val);
}
MatrixXd distEucSq(const MatrixXd &X, const MatrixXd &Y) {
  return (  X.rowwise().squaredNorm().eval().replicate(1,Y.rows())
  	  + Y.rowwise().squaredNorm().eval().transpose().replicate(X.rows(),1)
  	  - 2.0*X*Y.transpose() );
}
MatrixXd resistance_distance(const MatrixXd &M, const int o) {
  MatrixXd Hinv = - M; Hinv.diagonal() += M.rowwise().sum();
  Hinv.array() += 1.0;
  MatrixXd H = Hinv.inverse().topLeftCorner(o,o);
  MatrixXd R = - 2.0*H;
  VectorXd u_o = VectorXd::Ones(o);
  R.noalias() += H.diagonal()*u_o.transpose();
  R.noalias() += u_o*H.diagonal().transpose();
  return (R);
}
MatrixXd expected_dissimilarities(const MatrixXd &J, const MatrixXd &M, const VectorXd &q) {
  int n = J.rows();
  int o = J.cols();
  VectorXd Jq = J*q.head(o);
  VectorXd u_n = VectorXd::Ones(n);
  MatrixXd Delta = J*resistance_distance(M,o)*J.transpose();
  Delta.noalias() += 0.5*Jq*u_n.transpose();
  Delta.noalias() += 0.5*u_n*Jq.transpose();
  Delta.diagonal() -= Jq;
  return (Delta);
}
VectorXd split(const string &line) {
  istringstream in(line);
  vector<double> numbers;
  double number;
  while (!in.eof()) {
    in >> number;
    if (!in.fail()) { numbers.push_back(number); } else { break; }
  }
  if (in.fail()||in.bad()) { return (VectorXd::Zero(0)); }
  // Since we don't know how many numbers there are,
  // first we store the data in a std::vector of doubles,
  // and then typecast it to Eigen::VectorXd
  return (VectorXd::Map(&numbers[0],numbers.size()));
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
      // Resize the matrix to add each row (done once at the start of EEMS, so okay)
      insertRow(mat,row);
    } else {
      return (MatrixXd::Zero(0,0));
    }
  }
  return(mat);
}
double trace_AxB(const MatrixXd &A, const MatrixXd &B) {
  return (A.cwiseProduct(B).sum());
}
int mod(const int a, const int b) { return (a % b); }
int neighbors_in_grid(const int r1, const int c1, int &r2, int &c2, const int pos,
		      const int nx, const int ny) {
  int alpha = r1 * nx + c1;
  int beta = -1; r2 = -1; c2 = -1;
  // mod(alpha  ,nx) > 0 means that alpha is not the first deme in the row
  // mod(alpha+1,nx) > 0 means that alpha is not the last deme in the row
  // r1 < (ny-1) means that alpha is not on the top row
  // r1 > 0 means that alpha is not on the bottom row
  // mod(alpha   ,2*nx) > 0 means that alpha is not the first deme in an odd row
  // mode(alpha+1,2*nx) > 0 means that alpha is not the last deme in an even row
  // mod(r1+1,2) == 1 means that alpha is on an odd row
  if        ( (pos==0) && (mod(alpha  ,nx)>0) ) { 
    r2 = r1; c2 = c1 - 1;
  } else if ( (pos==3) && (mod(alpha+1,nx)>0) ) { 
    r2 = r1; c2 = c1 + 1;
  } else if ( (pos==5) && (r1>0 && mod(alpha,2*nx)>0) ) {
    r2 = r1 - 1; c2 = c1 - mod(r1+1,2);
  } else if ( (pos==4) && (r1>0 && mod(alpha+1,2*nx)>0) ) {       
    r2 = r1 - 1; c2 = c1 + 1 - mod(r1+1,2);
  } else if ( (pos==1) && (r1<(ny-1) && mod(alpha,2*nx)>0) ) {  
    r2 = r1 + 1; c2 = c1 - mod(r1+1,2);
  } else if ( (pos==2) && (r1<(ny-1) && mod(alpha+1,2*nx)>0) ) {  
    r2 = r1 + 1; c2 = c1 + 1 - mod(r1+1,2);
  }
  if ((r2>=0)&&(c2>=0)) { beta = nx*r2 + c2; }
  return (beta);
}
// If there are two draws and the first has two tiles and the second -- three tiles,
// then sizes = c(2,3) and array = c(m_{1t_1},m_{1t_2},m_{2t_1},m_{2t_2},m_{2t_3})
bool dlmcell(const string &filename, const VectorXd &sizes, const vector<double> &array) {
  if (array.size()!=sizes.sum()) { return false; }
  ofstream out(filename.c_str(),ofstream::out);
  if (!out.is_open()) { return false; }
  vector<double>::const_iterator it = array.begin();
  for ( int i = 0 ; i < sizes.size() ; i++ ) {
    for ( int j = 0 ; j < sizes(i) ; j++ ) {
      out << fixed << setprecision(6) << *it << " "; it++;
    }
    out << endl;
  }
  out.close( );
  return true;
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
  int rows = mat.rows();
  int cols = mat.cols();
  mat.noalias() = (MatrixXd(rows+1,cols) << mat, row.transpose()).finished();
}
VectorXd slice(const VectorXd &A, const VectorXi &I) {
  int elems = I.size();
  VectorXd B(elems);
  assert(I.minCoeff() >= 0);
  assert(I.maxCoeff() < A.size());
  for ( int i = 0 ; i < elems ; i++ ) {
    B(i) = A(I(i));
  }
  return (B);
}
MatrixXd slice(const MatrixXd &A, const VectorXi &R, const VectorXi &C) {
  int rows = R.size();
  int cols = C.size();
  MatrixXd B(rows,cols);
  assert(R.minCoeff() >= 0);
  assert(C.minCoeff() >= 0);
  assert(R.maxCoeff() < A.rows());
  assert(C.maxCoeff() < A.cols());
  for ( int i = 0 ; i < rows ; i++ ) {
  for ( int j = 0 ; j< cols ; j++ ) {
    B(i,j) = A(R(i),C(j));
  } }
  return (B);
}
