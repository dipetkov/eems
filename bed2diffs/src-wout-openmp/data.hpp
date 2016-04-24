
#include <fstream>
#include <iomanip>
#include <iostream>
#include <algorithm>

#include <plinkio/plinkio.h>
#define PLINK_NA 3


class Data {
public:
      
  size_t nIndiv, nSites;
  struct pio_file_t plink_file;
  std::string datapath;
      
  Data(std::string datapath);
  ~Data();

  void getsize();
  void bed2diffs_v1();
  void bed2diffs_v2();

protected:
  
  size_t Index(size_t i,size_t j);

};
