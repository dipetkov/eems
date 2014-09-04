
#include "data.hpp"

Data::Data(std::string datapath)
{
  nSites = 0;
  nIndiv = 0;
  this->datapath = datapath;
}

Data::~Data( ) { }

void Data::getsize( )
{
  std::stringstream errMsg;

  if( pio_open( &plink_file, datapath.c_str() ) != PIO_OK )
    {
      errMsg << "[Data::getsize] Error opening plink files " << datapath << ".[bed/bim/fam]" << std::endl;
      throw std::ios_base::failure(errMsg.str());
    }

  if( !pio_one_locus_per_row( &plink_file ) ) 
    {
      errMsg << "[Data::getsize] This program requires plink files [bed/bim/fam] in SNP-major mode" << std::endl;
      throw std::ios_base::failure(errMsg.str());
    }

  nIndiv = pio_num_samples( &plink_file );
  nSites = pio_num_loci( &plink_file );

  std::cout << "Detected plink files " << datapath << ".[bed/bim/fam] with " << nIndiv << " samples and " << nSites << " SNPs" << std::endl;
}

// Expects PLINK BED in SNP-major format
void Data::bed2diffs()
{
  std::string outfile = datapath + ".diffs";
  std::ofstream out(outfile.c_str(), std::ios::out);
  std::stringstream errMsg;

  RowVectorXd zvec(nIndiv);
  RowVectorXd obsrv(nIndiv);
  RowVectorXd ones = RowVectorXd::Ones(nIndiv);
  MatrixXd Diffs = MatrixXd::Zero(nIndiv,nIndiv);
  MatrixXd Counts = MatrixXd::Zero(nIndiv,nIndiv);
  MatrixXd ztonev(nIndiv,nIndiv);
  MatrixXd diffs(nIndiv,nIndiv);
  MatrixXd diffsSq(nIndiv,nIndiv);
  MatrixXd obsrvBoth(nIndiv,nIndiv);

  snp_t *snp_buffer = (snp_t *) malloc( pio_row_size( &plink_file ) );
  while ( pio_next_row( &plink_file, snp_buffer ) == PIO_OK )
    {

      for (unsigned int j = 0 ; j < nIndiv ; j++) 
	{
	  
	  /* Do not impute missing genotypes:
	     For a pair $i \neq j$, the average pairwise difference is
	     \begin{align}
	       D_{ij} = \frac{1}{|M_{ij}|} \sum_{l\in M_{ij}} (z_{li} - z_{lj})^2
	     \end{align}
	     where $M_{ij}$ is the set of loci where the genotype of both $i$ and $j$ is observed.
	     Thus no implicit imputation is applied by assuming that if $z_{li}$ is missing, then 
	     $z_{li} = \hat{p}_l$ where $\hat{p}_l$ is the observed frequency in the sample.
	  */
	  
	  int s = snp_buffer[ j ];
	  if(s != PLINK_NA) {
	    zvec(j) = s;
	    obsrv(j) = 1;
	  } else {
	    zvec(j) = s;
	    obsrv(j) = 0;
	  }
	}

      /*
         In matrix notation,
	 \begin{align}
	   1'z - z'1 = \Big( z_i - z_j \Big)
	 \end{align}
	 for all pairs $(i,j)$. Then
	 \begin{align}
	   (1'z - z'1) \circ (1'z - z'1) = \Big( (z_i - z_j)^2 \Big)
	 \end{align}
	 where $\circ$ is the elements-wise product.
       */

      ztonev = ones.transpose( ) * zvec;
      diffs = ztonev - ztonev.transpose( );
      diffsSq = diffs.cwiseProduct(diffs);

      obsrvBoth = obsrv.transpose( ) * obsrv;
      diffsSq = diffsSq.cwiseProduct(obsrvBoth);

      Diffs += diffsSq;
      Counts += obsrvBoth;
    }
  
  if (!out.is_open())
    {
      errMsg << "[Data::bed2diffs] Error writing file " << outfile << std::endl;
      throw std::ios_base::failure(errMsg.str());
    }   
  
  Diffs = Diffs.cwiseQuotient(Counts);
  out << Diffs << std::endl;
  
  free( snp_buffer );
  pio_close( &plink_file );
  out.close( );
}
