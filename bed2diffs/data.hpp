#pragma once

#include <fstream>
#include <iostream>
#include <algorithm>

#include <plinkio/plinkio.h>
#define PLINK_NA 3


class Data {
public:
      
  size_t nIndiv, nSites;
  struct pio_file_t plink_file;
  std::string datapath;
  int nthreads;
      
  Data(std::string datapath,int nthreads);
  ~Data();

  void getsize();
  void bed2diffs();

protected:
  
  size_t Index(size_t i,size_t j);

};
