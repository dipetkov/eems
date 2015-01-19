#include <cmath>
#include <vector>
#include <limits>
#include <iomanip>
#include <fstream>
#include <iostream>
using namespace std;

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
using namespace Eigen;

#include <boost/version.hpp>
#include <boost/limits.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/conversion/bounds.hpp>

#include <boost/config.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/detail/config_file.hpp>
#include <boost/algorithm/string.hpp>
namespace po = boost::program_options;


#ifndef UTILS
#define UTILS

const double Inf = numeric_limits<double>::infinity();
const double pi = boost::math::constants::pi<double>();
const double log_pi = log(boost::math::constants::pi<double>());
const double log_2 = boost::math::constants::ln_two<double>();

typedef Eigen::SparseMatrix<double> SpMat; // Declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> Tri;

class Params {
public:

  Params( );
  ~Params( );
  Params(const string &params_file);
  void dlmwrite(ostream &out) const;

  bool diploid;
  string datapath, mcmcpath;
  double qEffctHalfInterval, mEffctHalfInterval, mrateMuHalfInterval;
  double mSeedsProposalS2, mSeedsProposalS2x, mSeedsProposalS2y;
  double qSeedsProposalS2, qSeedsProposalS2x, qSeedsProposalS2y;
  double qEffctProposalS2, mEffctProposalS2, mrateMuProposalS2;
  double mrateShape, mrateScale;
  double qrateShape, qrateScale;
  double s2locShape, s2locScale;
  double dfProposalS2,negBiProb;
  double df,dflob,dfupb,dfmin,dfmax;
  int numMCMCIter,numBurnIter,numThinIter;
  int nDemes,nIndiv,nSites,negBiSize;
  bool testing;
};

void get_boost_version(ostream& out);
void get_eigen_version(ostream& out);
VectorXd split(const string &line);

int mod(const int a, const int b);
int neighbors_in_grid(const int r1, const int c1, int &r2, int &c2, const int pos,
		      const int nx, const int ny);

bool isdistmat(const MatrixXd &A);
double logdet(const MatrixXd &A);
double mvgammaln(const double a, const int p);
double wishpdfln(const MatrixXd &X, const MatrixXd &Sigma, const double df);
MatrixXd distEucSq(const MatrixXd &X, const MatrixXd &Y);
MatrixXd resistance_distance(const MatrixXd& M, const int o);
MatrixXd expected_dissimilarities(const MatrixXd &J, const MatrixXd& M, const VectorXd& q);

MatrixXd readMatrixXd(const string &filename);
double trace_AxB(const MatrixXd &A, const MatrixXd &B);

bool dlmcell(const string &filename, const VectorXd &sizes, const vector<double> &array);
void removeRow(MatrixXd &matrix, const int rowToRemove);
void removeElem(VectorXd &vector, const int elemToRemove);
void insertRow(MatrixXd &mat, const VectorXd &row);
void insertElem(VectorXd &vec, const double &elem);

#endif
