#pragma once

#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdexcept>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigen>

#include <plinkio/plinkio.h>

//#define PACK_DENSITY 4
#define PLINK_NA 3

//#define PHENO_BINARY_12 0
//#define PHENO_CONTINUOUS 1
//#define PLINK_PHENO_MISSING -9

// The BED file magic numbers
//#define PLINK_OFFSET 3
//#define COVAR_ACTION_TRAIN_TEST 0
//#define COVAR_ACTION_TRAIN_ONLY 1
//#define COVAR_ACTION_TRAIN_TEST_STR "traintest"
//#define COVAR_ACTION_TRAIN_ONLY_STR "train"

/* 3 is 11 in binary, we need a 2 bit mask for each of the 4 positions */
//#define MASK0 3	  /* 3 << 2 * 0 */
//#define MASK1 12  /* 3 << 2 * 1 */
//#define MASK2 48  /* 3 << 2 * 2 */
//#define MASK3 192 /* 3 << 2 * 3 */

//#define BUFSIZE 100
//#define DATA_MODE_TRAIN 1
//#define DATA_MODE_TEST 2

using namespace Eigen;

class Data {
public:
      
  unsigned int nIndiv, nSites;
  struct pio_file_t plink_file;
  std::string datapath;
  bool verbose;
      
  Data(std::string datapath);
  ~Data();

  void getsize();
  void bed2diffs();

};
