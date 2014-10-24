
#include "data.hpp"

Data::Data(std::string datapath,int nthreads)
{
  nSites = 0;
  nIndiv = 0;
  this->datapath = datapath;
  this->nthreads = nthreads;
  if( pio_open( &plink_file, datapath.c_str() ) != PIO_OK )
    {
      std::cerr << "[Data::getsize] Error opening plink files " << datapath << ".[bed/bim/fam]" << std::endl;
      exit(1);
    }
  if( !pio_one_locus_per_row( &plink_file ) ) 
    {
      std::cerr << "[Data::getsize] This program requires plink files [bed/bim/fam] in SNP-major mode" << std::endl;
      exit(1);
    }
}

Data::~Data( ) { 
  pio_close( &plink_file );
}

void Data::getsize( )
{
  nIndiv = pio_num_samples( &plink_file );
  nSites = pio_num_loci( &plink_file );
  std::cout << "Detected plink dataset " << datapath << ".[bed/bim/fam] "
	    << "with " << nIndiv << " samples and " << nSites << " SNPs" << std::endl;
}

void Data::bed2diffs()
{
  std::string diffsfile = datapath + ".diffs";
  std::string orderfile = datapath + ".order";
  std::ofstream outdiffs(diffsfile.c_str(), std::ios::out);
  std::ofstream outorder(orderfile.c_str(), std::ios::out);

  outdiffs.precision(6);
  outdiffs.setf(std::ios::fixed,std::ios::floatfield);
  
  size_t nPairs = nIndiv*(nIndiv-1)/2;
  size_t nSitesProcessed = 0;

  // Read the genotypes in serial, compute the differences in parallel
  snp_t *snps = (snp_t *) malloc( sizeof(snp_t)*nIndiv );
  size_t *a = (size_t *) malloc( sizeof(size_t)*nPairs );
  size_t *b = (size_t *) malloc( sizeof(size_t)*nPairs );
  double *diffs = (double *) malloc( sizeof(double)*nPairs );
  double *pairs = (double *) malloc( sizeof(double)*nPairs );

  for (size_t i = 0 ; i<(nIndiv-1) ; i++ ) {
    for (size_t j = i+1 ; j<nIndiv ; j++ ) {
      size_t ij = Index(i,j);
      diffs[ij] = 0.0;
      pairs[ij] = 0.0;
      a[ij] = i;
      b[ij] = j;
    }
  }

  while (pio_next_row( &plink_file, snps ) == PIO_OK) {

    nSitesProcessed++;
    
    #pragma omp parallel for
    for (size_t ij = 0 ; ij < nPairs ; ij++ ) {
      size_t zi = snps[a[ij]];
      size_t zj = snps[b[ij]];
      if ((zi!=PLINK_NA)&&(zj!=PLINK_NA)) {
	diffs[ij] += (zi - zj)*(zi - zj);
	pairs[ij] += 1.0;
      }
    }
  }

  std::cout << "Computed average pairwise differences across " << nSitesProcessed << " SNPs" << std::endl;

  if (!outdiffs.is_open())
    {
      std::cerr << "[Data::bed2diffs] Error writing file " << diffsfile << std::endl;
      exit(1);
    }   
  if (!outorder.is_open())
    {
      std::cerr << "[Data::bed2diffs] Error writing file " << orderfile << std::endl;
      exit(1);
    }   

  for (size_t i = 0 ; i<nIndiv ; i++ ) {

    struct pio_sample_t *sample = pio_get_sample( &plink_file, i );
    outorder << sample->fid << " " << sample->iid << std::endl;

    for (size_t j = 0 ; j<nIndiv ; j++ ) {
      if (i==j) {
	outdiffs << " " << 0;
      } else {
	size_t ij = Index(i,j);
	outdiffs << " " << diffs[ij]/pairs[ij];
      }
    }
    outdiffs << std::endl;
  }
  
  outdiffs.close( );
  outorder.close( );
  free( a );
  free( b );
  free( snps );
  free( diffs );
  free( pairs );
}
size_t Data::Index(size_t i, size_t j) {
  size_t ij;
  /*
    For n = 5, return the following indices:
    25  0  1  2  3
     0 25  4  5  6
     1  4 25  7  8
     2  5  7 25  9
     3  6  8  9 25
  */
  if (i==j) {
    ij = nIndiv*nIndiv;
  } else if (i<j) {
    ij = nIndiv*i - i*(i+3)/2 + j - 1;
  } else {
    ij = nIndiv*j - j*(j+3)/2 + i - 1;
  }
  return(ij);
}
