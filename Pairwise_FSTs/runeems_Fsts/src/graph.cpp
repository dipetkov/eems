
#include "graph.hpp"

Graph::Graph( ) { }
Graph::~Graph( ) { }
void Graph::generate_grid(const string &datapath, const string &gridpath, const Habitat &habitat,
			  const int nDemeDensity, const int nIndiv)
{
  cerr << "[Graph::initialize]" << endl;
  // These two functions can be rewritten to load a pre-constructred population graph
  // First, generate the deme coordinates, DemeCoord, and the list of edges, DemePairs
  // Then, load the sample coordinates and map each sample to the closest deme
  if (gridpath.empty()) {
    cerr << "  Generate population grid and sample assignment"<< endl;
    make_triangular_grid(habitat,nDemeDensity);
    map_indiv_to_deme(datapath,nIndiv);
  } else {
    cerr << "  Load population grid and sample assignment from " << datapath << endl;
    // Read the population grid (demes and edges)
    if (!read_input_grid(datapath,DemeCoord,DemePairs)) {
      cerr << "  Error reading population grid." << endl;
      exit(1);
    }
    // Read the assignment of individuals to demes
    if (!read_indiv_to_deme(datapath,DemeCoord.rows(),indiv2deme)) {
      cerr << "  Error reading sample assignment." << endl;
      exit(1);
    }
  }
  if (!this->is_connected()) {
    cerr << "  The population grid is not connected." << endl; exit(1);
  }
  // Finally, reorder the demes, whether the graph is constructed on the fly or loaded from files
  this->reindex_demes();
  int nDemes = this->get_num_total_demes();
  int oDemes = this->get_num_obsrv_demes();
  int nEdges = this->get_num_edges();
  cerr << "  The population grid has " << nDemes << " demes and " << nEdges << " edges" << endl
       << "  There are " << nIndiv << " samples assigned to " << oDemes << " observed demes" << endl;
  cerr << "[Graph::initialize] Done." << endl << endl;
}
/*
  Construct a regular triangular grid, entirely contained inside the habitat outline
 */
void Graph::make_triangular_grid(const Habitat &habitat, const int nDemeDensity) {
  double xspan = habitat.get_xspan();
  double yspan = habitat.get_yspan();
  double area = habitat.get_area();
  int xDemes = (int)sqrt(nDemeDensity*xspan*xspan/area);
  int yDemes = (int)sqrt(nDemeDensity*yspan*yspan/area);
  double scalex = 1.0;
  double scaley = 1.0;
  // A triangular grid extends half a triangle on the right
  if (xDemes>1) { scalex = xspan/(xDemes-0.5); }
  if (yDemes>1) { scaley = yspan/(yDemes-1.0); }
  // The new index will be -1 for demes to exclude because they fall outside the habitat
  // The rest of the demes will be assigned an unique index, starting with 0
  VectorXi newIndex = -1 * VectorXi::Ones(xDemes*yDemes);
  int nDemes = 0;
  int nEdges = 0;
  // Suggest with a regular triangular xDemes-by-yDemes grid within the rectangle (xmin,ymin)
  // -- (xmax,ymax) but only include demes that fall inside the habitat boundary
  // This first double loop counts the number of demes and edges, and orders the demes
  for ( int r1 = 0, r2 = 0 ; r1 < yDemes ; r1++ ) {
  for ( int c1 = 0, c2 = 0 ; c1 < xDemes ; c1++ ) {
    int alpha = r1 * xDemes + c1;
    double x1 = habitat.get_xmin() + scalex*(c1+0.5*(r1%2));
    double y1 = habitat.get_ymin() + scaley*(r1);
    if (habitat.in_point(x1,y1)) {
      newIndex[alpha] = nDemes++;
      for ( int pos = 0 ; pos < 6 ; pos++ ) {
	int beta = neighbors_in_grid(r1,c1,r2,c2,pos,xDemes,yDemes);
	double x2 = habitat.get_xmin() + scalex*(c2+0.5*(r2%2));
	double y2 = habitat.get_ymin() + scaley*(r2);
	// It is sufficient to enter each edge only once because migration is undirected
	if ((alpha<beta) && habitat.in_point(x2,y2)) { nEdges++; }
      }
    }
  } }
  DemeCoord.resize(nDemes,2);
  DemePairs.resize(nEdges,2);
  // This second double loop assigns coordinates to the demes, and connects neighboring demes
  for ( int r1 = 0, r2 = 0, e = 0 ; r1 < yDemes ; r1++ ) {
  for ( int c1 = 0, c2 = 0        ; c1 < xDemes ; c1++ ) {
    int alpha = r1 * xDemes + c1;
    double x1 = habitat.get_xmin() + scalex*(c1+0.5*(r1%2));
    double y1 = habitat.get_ymin() + scaley*(r1);
    if (habitat.in_point(x1,y1)) {
      DemeCoord(newIndex[alpha],0) = x1;
      DemeCoord(newIndex[alpha],1) = y1;
      for ( int pos = 0 ; pos < 6 ; pos++ ) {
	int beta = neighbors_in_grid(r1,c1,r2,c2,pos,xDemes,yDemes);
	double x2 = habitat.get_xmin() + scalex*(c2+0.5*(r2%2));
	double y2 = habitat.get_ymin() + scaley*(r2);
	if ((alpha<beta) && habitat.in_point(x2,y2)) {
	  DemePairs(e,0) = newIndex[alpha];
	  DemePairs(e,1) = newIndex[beta]; e++;
	}
      }
    }
  } }
}
/*
  Read the sample coordinates and assign each sample to its closest deme
 */
void Graph::map_indiv_to_deme(const string &datapath, const int nIndiv) {
  MatrixXd IndivCoord = readMatrixXd(datapath + ".coord");
  if (IndivCoord.rows()!=nIndiv || IndivCoord.cols()!=2) {
    cerr << "  Error reading sample coordinates from " << datapath + ".coord" << endl
	 << "  Expect a list of " << nIndiv << " points, with two coordinates per line" << endl; exit(1);
  }
  cerr << "  Loaded sample coordinates from " << datapath + ".coord" << endl;
  indiv2deme.resize(nIndiv); // the deme to which each sample is assigned (to be determined)
  for ( int i = 0, alpha = 0 ; i < nIndiv ; i++ ) {
    pairwise_distance(DemeCoord,IndivCoord.row(i)).col(0).minCoeff( &alpha );
    indiv2deme(i) = alpha;
  }
}
/*
  Assign new indices to the demes, so that
   * the observed demes have indices from 0 to oDemes-1, and 
   * the unobserved demes have indices from oDemes to nDemes-1
  Then we can split the nDemes-by-nDemes matrix of pairwise distances between demes
  into two contiguous blocks, which correspond to observed and unobserved demes, resp.
  This makes it easier to exploit the Schur decomposition trick to compute distances
 */
void Graph::reindex_demes( ) {
  int nDemes = DemeCoord.rows();
  int nEdges = DemePairs.rows();
  int oDemes = 0;
  int nIndiv = indiv2deme.size();
  if ( (DemePairs.maxCoeff()>=nDemes) || (DemePairs.minCoeff()<0) ) {
    cerr << "[Graph::reindex_demes] DemePairs input error." << endl; exit(1);
  }
  if ( (indiv2deme.maxCoeff()>=nDemes) || (indiv2deme.minCoeff()<0) ) {
    cerr << "[Graph::reindex_demes] indiv2deme input error." << endl; exit(1);
  }
  VectorXi newIndex = VectorXi::Constant(nDemes,-1);
  for ( int i = 0, alpha = 0 ; i < nIndiv ; i++ ) {
    alpha = indiv2deme(i);
    if (newIndex(alpha) == -1) {
      newIndex(alpha) = oDemes++;
    }
    indiv2deme(i) = newIndex(alpha);
  }
  // Count the number of samples are taken from each observed deme
  DemeSizes = VectorXi::Zero(oDemes);
  for ( int i = 0 ; i < nIndiv ; ++i ) {
    DemeSizes(indiv2deme(i))++;
  }
  cerr << "  There are " << oDemes << " observed demes (out of " << nDemes << " demes)" << endl;
  // So far, have assigned new indices to the observed demes
  // Now assign indices to the unobserved demes, which run from oDemes to nDemes-1
  for ( int i = 0 ; i < nDemes ; i++ ) {
    if (newIndex(i) == -1) { newIndex(i) = oDemes++; }
  }
  // Now indiv2deme and DemeSizes use the new deme indices
  // But DemeCoords and DemePairs still use the old deme indices
  // Have to store the deme coordinates in a temporary matrix
  MatrixXd InputDemes = DemeCoord;
  for ( int i = 0 ; i < DemeCoord.rows( ) ; i++ ) {
    DemeCoord(newIndex(i),0) = InputDemes(i,0);
    DemeCoord(newIndex(i),1) = InputDemes(i,1);
  }
  for ( int i = 0 ; i < DemePairs.rows( ) ; i++ ) {
    int alpha = newIndex(DemePairs(i,0));
    int beta = newIndex(DemePairs(i,1));
    DemePairs(i,0) = alpha;
    DemePairs(i,1) = beta;
  }
}
// The graph is connected if it has exactly one connected component
bool Graph::is_connected( ) const {
  boost::adjacency_list <boost::vecS,boost::vecS,boost::undirectedS> testG;
  for ( int i = 0 ; i < DemePairs.rows() ; i++ ) {
    add_edge(DemePairs(i,0),DemePairs(i,1),testG);
  }
  vector<int> component(num_vertices(testG));
  return (connected_components(testG,&component[0])==1);
}
bool Graph::dlmwrite_grid(const string &mcmcpath) const {
  ofstream out;
  // All the plotting is done in R, and indexing in R starts at 1
  out.open((mcmcpath + "/ipmap.txt").c_str(),ofstream::out);
  if (!out.is_open()) { return false; }
  out << indiv2deme.array() + 1 << endl;
  out.close( );
  out.open((mcmcpath + "/demes.txt").c_str(),ofstream::out);
  if (!out.is_open()) { return false; }
  out << DemeCoord << endl;
  out.close( );
  out.open((mcmcpath + "/edges.txt").c_str(),ofstream::out);
  if (!out.is_open()) { return false; }
  out << DemePairs.array() + 1 << endl;
  out.close( );
  return true;
}
int Graph::get_num_obsrv_demes() const { return (DemeSizes.size()); }
int Graph::get_num_total_demes() const { return (DemeCoord.rows()); }
int Graph::get_num_edges() const { return (DemePairs.rows()); }
int Graph::get_deme_of_indiv(const int i) const { return (indiv2deme(i)); }
void Graph::get_edge(int edge, int &alpha, int &beta) const
{
  if ( (edge >= 0) && (edge < this->get_num_edges() ) ) {
    alpha = DemePairs(edge,0); beta = DemePairs(edge,1);
  } else {
    alpha = beta = -1;
  }
}
MatrixXd Graph::get_the_obsrv_demes() const
{
  return (DemeCoord.topLeftCorner(this->get_num_obsrv_demes(),2));
}
/*
  Compute the pairwise distances of each deme in DemeCoord to each seed in Seeds.
  Once EEMS is initialized, DemeCoord is fixed but Seeds changes as tiles are added and removed.
 */
void Graph::index_closest_to_deme(const MatrixXd &Seeds, VectorXi &Closest) const
{
  int nDemes = this->get_num_total_demes();
  if (Closest.size()!=nDemes) { Closest.resize(nDemes); }
  MatrixXd Dist = pairwise_distance(Seeds,DemeCoord);
  for ( int i = 0 ; i < nDemes ; i++ ) {
    Dist.col(i).minCoeff( &Closest(i) );
  }
}
int Graph::neighbors_in_grid(const int r1, const int c1, int &r2, int &c2, const int pos,
			     const int nx, const int ny) const {
  int alpha = r1 * nx + c1;
  int beta = -1; r2 = -1; c2 = -1;
  // mod(alpha  ,nx) > 0 means that alpha is not the first deme in the row
  // mod(alpha+1,nx) > 0 means that alpha is not the last deme in the row
  // r1 < (ny-1) means that alpha is not on the top row
  // r1 > 0 means that alpha is not on the bottom row
  // mod(alpha   ,2*nx) > 0 means that alpha is not the first deme in an odd row
  // mode(alpha+1,2*nx) > 0 means that alpha is not the last deme in an even row
  // mod(r1+1,2) == 1 means that alpha is on an odd row
  if        ( (pos==0) && ( alpha   %nx>0) ) { 
    r2 = r1; c2 = c1 - 1;
  } else if ( (pos==3) && ((alpha+1)%nx>0) ) { 
    r2 = r1; c2 = c1 + 1;
  } else if ( (pos==5) && (r1>0) && ( alpha   %(2*nx)>0) ) {
    r2 = r1 - 1; c2 = c1 - (r1+1)%2;
  } else if ( (pos==4) && (r1>0) && ((alpha+1)%(2*nx)>0) ) {
    r2 = r1 - 1; c2 = c1 + 1 - (r1+1)%2;
  } else if ( (pos==1) && (r1<(ny-1)) && (alpha    %(2*nx)>0) ) {  
    r2 = r1 + 1; c2 = c1 - (r1+1)%2;
  } else if ( (pos==2) && (r1<(ny-1)) && ((alpha+1)%(2*nx)>0) ) {  
    r2 = r1 + 1; c2 = c1 + 1 - (r1+1)%2;
  }
  if ((r2>=0)&&(c2>=0)) { beta = nx*r2 + c2; }
  return (beta);
}
// Read the population grid (demes and edges)
// The grid to read does not need to be triangular as long as it is connected
bool Graph::read_input_grid(const string &datapath, MatrixXd &DemeCoord, MatrixXi &DemePairs) {
  // Read the deme coordinates
  DemeCoord = readMatrixXd(datapath + ".demes");
  if (DemeCoord.cols()!=2) { return (false); }
  // Read the connected pairs of demes
  MatrixXd tempi = readMatrixXd(datapath + ".edges");
  if (tempi.cols()!=2) { exit(1); }
  DemePairs = tempi.cast<int>();
  DemePairs = DemePairs.array() - 1;
  int nDemes = DemeCoord.rows();
  return (DemePairs.minCoeff()>=0 && DemePairs.maxCoeff()<nDemes);
}
// Read the assignment of individuals to demes
bool Graph::read_indiv_to_deme(const string &datapath, const int nDemes, VectorXi &indiv2deme) {
  MatrixXd tempi = readMatrixXd(datapath + ".ipmap");
  if (tempi.cols()!=1) { exit(1); }
  indiv2deme = tempi.col(0).cast<int>();
  indiv2deme = indiv2deme.array() - 1;
  return (indiv2deme.minCoeff()>=0 && indiv2deme.maxCoeff()<nDemes);
}
