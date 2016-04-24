
#include <string>
#include "data.hpp"

static void show_usage( ) {
  std::cerr << "Usage: bed2diffs_v2 --bfile PlinkData\n"
	    << std::endl;
}


int main(int argc, char * argv[])
{
  
  std::string bfilepath;

  ////////////////////////////////////////////////////////////////////////////////
  // Parse commandline arguments
  if (argc<3) { show_usage( ); return EXIT_FAILURE;  }
  for (int i = 1 ; i < argc ; i++ ) {
    std::string arg = argv[i];
    if (arg == "--bfile") {
      if (i+1<argc) { bfilepath.append(argv[++i]); } else { show_usage( ); return EXIT_FAILURE; }
    }
  }
  ////////////////////////////////////////////////////////////////////////////////
  // End command line parsing
   
  std::cout << "Compute the average genetic differences according to: " << std::endl
	    << "  Dij = (1/|Mtot|) sum_{m in Mtot} (z*_{im} - z*_{jm})^2" << std::endl
	    << "  where Mtot is the set of all SNPs and" << std::endl
	    << "  z*_{im} = z_{im} if sample i is called at marker m" << std::endl
	    << "          = zbar_m (the average genotype at m) otherwise" << std::endl
	    << std::endl;
   
  Data data(bfilepath);
  data.getsize();
  data.bed2diffs_v2();
   
  return EXIT_SUCCESS;
}
