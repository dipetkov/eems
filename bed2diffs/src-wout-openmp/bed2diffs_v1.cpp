
#include <string>
#include "data.hpp"

static void show_usage( ) {
  std::cerr << "Usage: bed2diffs_v1 --bfile PlinkData\n"
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
	    << "  Dij = (1/|Mij|) sum_{m in Mij} (z_{im} - z_{jm})^2" << std::endl
	    << "  where Mij is the set of SNPs where both i and j are called" << std::endl
	    << std::endl;
  
  Data data(bfilepath);
  data.getsize();
  data.bed2diffs_v1();
   
  return EXIT_SUCCESS;
}
