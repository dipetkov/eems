#pragma once

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
#include <boost/filesystem.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/conversion/bounds.hpp>

#include <boost/config.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/detail/config_file.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/algorithm/string.hpp>
namespace po = boost::program_options;

#ifndef UTIL_H
#define UTIL_H

const double Inf = numeric_limits<double>::infinity();
const double pi = boost::math::constants::pi<double>();
const double log_2 = boost::math::constants::ln_two<double>();
const double log_pi = log(boost::math::constants::pi<double>());
const double pi_180 = pi / 180.0;

typedef Eigen::SparseMatrix<double> SpMat; // Declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> Tri;

class Params {
public:

  Params( );
  ~Params( );
  Params(const string &params_file, const long seed_from_command_line);
  bool check_input_params( ) const;
  
  friend ostream& operator<<(ostream& out, const Params& params);
  
  long seed;
  bool diploid, testing;
  string datapath, mcmcpath, prevpath, gridpath;
  double qEffctHalfInterval, mEffctHalfInterval, mrateMuHalfInterval;
  double mSeedsProposalS2, mSeedsProposalS2x, mSeedsProposalS2y;
  double qSeedsProposalS2, qSeedsProposalS2x, qSeedsProposalS2y;
  double qEffctProposalS2, mEffctProposalS2, mrateMuProposalS2;
  double mrateShape_2, mrateScale_2;
  double qrateShape_2, qrateScale_2;
  double sigmaShape_2, sigmaScale_2;
  double qVoronoiPr, negBiProb;
  int numMCMCIter, numBurnIter, numThinIter;
  int nDemes, nIndiv, nSites, negBiSize;
  string distance;
};

VectorXd split(const string &line);
bool isdistmat(const MatrixXd &A);
double logdet(const MatrixXd &A);
double pseudologdet(const MatrixXd &A, const int rank);
double pseudowishpdfln(const MatrixXd &X, const MatrixXd &Sigma, const int df);
double mvgammaln(const double a, const int p);
double wishpdfln(const MatrixXd &X, const MatrixXd &Sigma, const double df);
MatrixXd pairwise_distance(const MatrixXd &X, const MatrixXd &Y);
MatrixXd resistance_distance(const MatrixXd& M, const int o);
MatrixXd expected_dissimilarities(const MatrixXd &J, const MatrixXd& M, const VectorXd& q);
MatrixXd readMatrixXd(const string &filename);
double trace_AxB(const MatrixXd &A, const MatrixXd &B);

bool dlmcell(const string &filename, const VectorXd &sizes, const vector<double> &array);
void removeRow(MatrixXd &matrix, const int rowToRemove);
void removeElem(VectorXd &vector, const int elemToRemove);
void insertRow(MatrixXd &mat, const VectorXd &row);
void insertElem(VectorXd &vec, const double &elem);

double dnegbinln(const int k, const int size, const double prob);
double dinvgamln(const double x, const double shape, const double scale);
double dmvnormln(const VectorXd &x, const VectorXd &mu, const MatrixXd &sigma);
double dtrnormln(const double x, const double mu, const double sigma2, const double bnd);

VectorXd slice(const VectorXd &A, const VectorXi &I);
MatrixXd slice(const MatrixXd &A, const VectorXi &R, const VectorXi &C);

MatrixXd greatcirc_dist(const MatrixXd &X, const MatrixXd &Y);
MatrixXd euclidean_dist(const MatrixXd &X, const MatrixXd &Y);

#endif
