
#include <string>
#include <fstream>
#include <sstream>
#include "data.hpp"

using namespace Eigen;

int main(int argc, char * argv[])
{
  
  ////////////////////////////////////////////////////////////////////////////////
  // Parse commandline arguments
  if ( argc!=3 || strcmp(argv[1],"--bfile")!=0 ) {
    std::cerr << "Usage: " << argv[0] << " --bfile PlinkData " << std::endl;
    return 1;
  }
  std::string bfilepath(argv[2]);
  ////////////////////////////////////////////////////////////////////////////////
  // End command line parsing
      
  Data data(bfilepath);
  data.verbose = true;
  data.getsize();
  data.bed2diffs();

  return EXIT_SUCCESS;
}
