
#include "mcmc.hpp"

MCMC::MCMC(const Params &params) {
  numMCMCIter = params.numMCMCIter;
  numBurnIter = params.numBurnIter;
  numThinIter = params.numThinIter;
  currIter = 0;
  currStep = 0;
  numTypes = 8;
  iterDone = false;
  finished = false;
  okayMoves = vector<double>(numTypes,0);
  totalMoves = vector<double>(numTypes,0);
}
MCMC::~MCMC( ) { }
void MCMC::start_iteration() {
  currIter++;
  currStep = 0;
  iterDone = false;
}
int MCMC::num_iters_to_save( ) const {
  int a = (numMCMCIter - numBurnIter) / (numThinIter + 1);
  return (a);
}
int MCMC::to_save_iteration( ) const {
  if (currIter>numBurnIter) {
    int a = (currIter - numBurnIter) / (numThinIter + 1);
    int b = (currIter - numBurnIter) % (numThinIter + 1);
    if (b==0) { return (a-1); }
  } 
  return (-1);
}
ostream& operator<<(ostream& out, const MCMC& mcmc) {
  for ( int i = 0 ; i <= 6 ; i++ ) {
    double a = mcmc.okayMoves.at(i);
    double A = mcmc.totalMoves.at(i);
    out << setprecision(2) << "\t(" << (int)a << "/" << (int)A << ") = " << 100.0*(a/A) << "% for ";
    switch (i) {
    case 0:
      out << "type 0 (qEffcts)" << endl;
      break;
    case 1:
      out << "type 1 (qVoronoi)" << endl;
      break;
    case 2:
      out << "type 2 (qBirthDeath)" << endl;
      break;
    case 3:
      out << "type 3 (mEffcts)" << endl;
      break;
    case 4:
      out << "type 4 (mrateMu)" << endl;
      break;
    case 5:
      out << "type 5 (mVoronoi)" << endl;
      break;
    case 6:
      out << "type 6 (mBirthDeath)" << endl;
      break;
    default:
      out << "type 7 (df)" << endl;
    }
  }
  return out;
}
void MCMC::end_iteration( ) {
  if (currIter==numMCMCIter) {
    finished = true;
  }
}
void MCMC::change_update( ) {
  switch (currStep) {
  case 4:
    currStep = 0;
    iterDone = true;
    break;
  default:
    currStep++;
  }
}
void MCMC::add_to_okay_moves(const int type) {
  okayMoves[type]++;
}
void MCMC::add_to_total_moves(const int type) {
  totalMoves[type]++;
}
