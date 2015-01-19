
#include "graph.hpp"

Graph::Graph( ) { }
Graph::~Graph( ) { }
void Graph::initialize(const string &datapath, const Habitat &habitat,
		       const int nDemesSuggested, const int nIndiv)
{
  cerr << "[Graph::initialize]" << endl;
  double xspan = habitat.get_xspan();
  double yspan = habitat.get_yspan();
  double area = habitat.get_area();
  int xDemes = (int)sqrt(nDemesSuggested*xspan*xspan/area);
  int yDemes = (int)sqrt(nDemesSuggested*yspan*yspan/area);
  double scalex = 1.0;
  double scaley = 1.0;
  // A triangular grid extends half a triangle on the right
  if (xDemes>1) { scalex = xspan/(xDemes-0.5); }
  if (yDemes>1) { scaley = yspan/(yDemes-1.0); }
  // The new index will be -1 for demes to exclude because they fall outside the habitat
  // The rest of the demes will be assigned an unique index, starting with 0
  VectorXi newIndex = -1 * VectorXi::Ones(xDemes*yDemes);
  int nDemes1 = 0; int d = 0;
  int nEdges1 = 0; int e = 0;
  int r1 = -1, r2 = -1, c1 = -1, c2 = -1, r;
  // Suggest a regular triangular grid with xDemes in the x (longitude) direction and
  // yDemes in the y (latitude) direction, but only include demes that fall inside the habitat
  // (which might not be rectangular in shape)
  // This first double loop counts the number of demes and edges, and orders the demes
  for ( r1 = 0 ; r1 < yDemes ; r1++ ) {
  for ( c1 = 0 ; c1 < xDemes ; c1++ ) {
    int alpha = r1 * xDemes + c1;
    double x1 = habitat.get_xmin() + scalex*(c1+0.5*mod(r1,2));
    double y1 = habitat.get_ymin() + scaley*(r1);
    if (habitat.in_point(x1,y1)) {
      newIndex[alpha] = nDemes1++;
      for ( int pos = 0 ; pos < 6 ; pos++ ) {
	int beta = neighbors_in_grid(r1,c1,r2,c2,pos,xDemes,yDemes);
	double x2 = habitat.get_xmin() + scalex*(c2+0.5*mod(r2,2));
	double y2 = habitat.get_ymin() + scaley*(r2);
	// It is sufficient to enter each edge only once
	// in 'Pairs' because migration is undirected
	if ((alpha<beta) && habitat.in_point(x2,y2)) { nEdges1++; }
      }
    }
  } }
  MatrixXd Demes1 = -1 * MatrixXd::Ones(nDemes1,2);
  MatrixXi Edges1 = -1 * MatrixXi::Ones(nDemes1,6);
  MatrixXi Pairs1 = -1 * MatrixXi::Ones(nEdges1,2);
  // This second double loop assigns coordinates to the demes, and connects neighboring demes
  for ( r1 = 0 ; r1 < yDemes ; r1++ ) {
  for ( c1 = 0 ; c1 < xDemes ; c1++ ) {
    int alpha = r1 * xDemes + c1;
    double x1 = habitat.get_xmin() + scalex*(c1+0.5*mod(r1,2));
    double y1 = habitat.get_ymin() + scaley*(r1);
    if (habitat.in_point(x1,y1)) {
      Demes1(newIndex[alpha],0) = x1;
      Demes1(newIndex[alpha],1) = y1; d++;
      for ( int pos = 0 ; pos < 6 ; pos++ ) {
	int beta = neighbors_in_grid(r1,c1,r2,c2,pos,xDemes,yDemes);
	double x2 = habitat.get_xmin() + scalex*(c2+0.5*mod(r2,2));
	double y2 = habitat.get_ymin() + scaley*(r2);
	// It is sufficient to enter each edge only once
	// in 'Pairs' because migration is undirected
	if (habitat.in_point(x2,y2)) {
	  Edges1(newIndex[alpha],pos) = newIndex[beta];
	}
	if ((alpha<beta) && habitat.in_point(x2,y2)) {
	  Pairs1(e,0) = newIndex[alpha];
	  Pairs1(e,1) = newIndex[beta]; e++;
	}
      }
    }
  } }
  nDemes = Demes1.rows();
  nEdges = Pairs1.rows();
  cerr << "  Created a population grid with " << nDemes << " demes and " << nEdges << " edges" << endl;
  /////////////////////////////////////////////
  // Will assign new indices to the demes, so that:
  //    the observed demes have indices from 1 to oDemes and 
  //    the unobserved demes have indices from oDemes+1 to nDemes
  // This will make it easier to exploit the Schur decomposition trick to get the resistances
  // between observed demes
  newIndex = -1 * VectorXi::Ones(nDemes);
  i2alpha = -1 * VectorXi::Ones(nIndiv); // the deme to which each sample is assigned (to be determined)
  Coord = readMatrixXd(datapath + ".coord"); // the coordinates of the sampled individuals
  if (Coord.rows()==0||Coord.cols()!=2||Coord.rows()!=nIndiv) {
    cerr << "  Error reading sample coordinates from " << datapath + ".coord" << endl
	 << "  Expect a list of " << nIndiv << " points, with two coordinates per line" << endl; exit(1);
  }
  cerr << "  Read sample coordinates from " << datapath + ".coord" << endl;
  oDemes = 0;
  for ( int i=0 ; i<nIndiv ; i++ ) {
    // Assign the sample to the closest deme -- should we be using great circle distances instead?
    distEucSq(Demes1,Coord.row(i)).col(0).minCoeff( &r );
    if (newIndex(r) == -1) {
      newIndex(r) = oDemes++;
    }
    i2alpha(i) = newIndex(r);
  }
  // Count the number of samples are taken from each observed deme
  Counts = VectorXi::Zero(oDemes);
  for ( int i = 0 ; i < nIndiv ; ++i ) { Counts(i2alpha(i))++; }
  cerr << "  There are " << oDemes << " observed demes (out of " << nDemes << " demes)" << endl;
  // So far, have assigned new indices to the observed demes. Now assign indices to the unobserved demes.
  r = oDemes; // unobserved demes have indices from oDemes to nDemes-1
  for ( int j = 0 ; j < nDemes ; j++ ) { if (newIndex(j)<0) { newIndex(j) = r++; } }  
  // At this point, i2alpha and Counts use the new deme indices
  // But Demes, Edges and Pairs use the old deme indices
  Demes = -1 * MatrixXd::Ones(nDemes,2);
  Edges = -1 * MatrixXi::Ones(nDemes,6);
  Pairs = -1 * MatrixXi::Ones(nEdges,2);
  reindex_demes(Demes1, Edges1, Pairs1, newIndex);
  /////////////////////////////////////////////
  // Make sure the resulting graph is connected
  BoostGraph testG;
  for ( int i = 0 ; i < nEdges ; i++ ) {
    add_edge((int)Pairs(i,0),(int)Pairs(i,1),testG);
  }
  vector<int> component(num_vertices(testG));
  int num = connected_components(testG,&component[0]);
  if (num!=1) {
    cerr << "  The generated population grid is not connected." << endl; exit(1);
  }
  cerr << "[Graph::initialize] Done." << endl << endl;
}
void Graph::reindex_demes(const MatrixXd &Demes1, const MatrixXi &Edges1, const MatrixXi &Pairs1,
			  const VectorXi &newIndex)  {
  bool err = false;
  // Re-indexing does not change the size of the graph
  if (Demes1.rows()!=nDemes) { err = true; }
  if (Edges1.rows()!=nDemes) { err = true; }
  if (Pairs1.rows()!=nEdges) { err = true; }
  int nDemes = Demes1.rows();
  int nEdges = Pairs1.rows();
  ArrayXi tempIndex = newIndex.array();
  for ( int i = 0 ; i < nDemes ; i++ ) {
    if ((tempIndex==i).sum()!=1) { err = true; }
  }
  if (Edges1.maxCoeff()>=nDemes) { err = true; }
  if (Pairs1.maxCoeff()>=nDemes) { err = true; }
  if (Edges1.minCoeff()<-1) { err = true; }
  if (Pairs1.minCoeff()< 0) { err = true; }
  if (err) { cerr << "[Graph::reindex_demes] Input error." << endl; exit(1); }
  for ( int i = 0 ; i < Demes1.rows( ) ; i++ ) {
    Demes(newIndex(i),0) = Demes1(i,0);
    Demes(newIndex(i),1) = Demes1(i,1);
    for ( int pos = 0 ; pos < 6 ; pos++ ) {
      if (Edges1(i,pos) < 0) {
	Edges(newIndex(i),pos) = -1;
      } else {
	Edges(newIndex(i),pos) = newIndex(Edges1(i,pos));
      }
    }
  }
  for ( int i = 0 ; i < Pairs1.rows( ) ; i++ ) {
    Pairs(i,0) = newIndex(Pairs1(i,0));
    Pairs(i,1) = newIndex(Pairs1(i,1));
  }
}
bool Graph::dlmwrite(const string &mcmcpath) const {
  ofstream out;
  // All the plotting is done in R, and indexing in R starts at 1
  out.open((mcmcpath + "/edges.txt").c_str(),ofstream::out);
  if (!out.is_open()) { return false; }
  out << Edges.array() + 1 << endl;
  out.close( );
  out.open((mcmcpath + "/ipmap.txt").c_str(),ofstream::out);
  if (!out.is_open()) { return false; }
  out << i2alpha.array() + 1 << endl;
  out.close( );
  out.open((mcmcpath + "/demes.txt").c_str(),ofstream::out);
  if (!out.is_open()) { return false; }
  out << Demes << endl;
  out.close( );
  return true;
}
MatrixXd Graph::get_the_obsrv_demes() const { return (Demes.topLeftCorner(oDemes,2)); }
int Graph::get_num_obsrv_demes() const { return (oDemes); }
int Graph::get_num_total_demes() const { return (nDemes); }
int Graph::get_num_edges() const { return (nEdges); }
int Graph::get_deme_of_indiv(const int i) const { return (i2alpha(i)); }
void Graph::get_edge(int edge, int &alpha, int &beta) const
{
  alpha = Pairs(edge,0); beta = Pairs(edge,1);
}
void Graph::pdist2(const MatrixXd &Seeds, VectorXi &closest) const
{
  if (closest.size()!=nDemes) { closest.resize(nDemes); }
  MatrixXd Dist = distEucSq(Seeds,Demes);
  for ( int i = 0 ; i < nDemes ; i++ ) {
    Dist.col(i).minCoeff( &closest(i) );
  }
}
