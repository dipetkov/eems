#pragma once

#include "util.hpp"
#include "habitat.hpp"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

typedef boost::adjacency_list <boost::vecS,boost::vecS,boost::undirectedS> BoostGraph;

// Triangular graph, possibly irregularly shaped. The shape is determined by the habitat outline.

class Graph {
public:

  Graph( );
  ~Graph( );

  void initialize(const string &datapath, const Habitat &habitat,
		  const int nDemesSuggested, const int nIndiv);
  bool dlmwrite(const string &mcmcpath) const;

  MatrixXd get_the_obsrv_demes() const;
  int get_num_edges() const;
  int get_num_obsrv_demes() const;
  int get_num_total_demes() const;
  int get_deme_of_indiv(const int i) const;
  void get_edge(int edge, int &alpha, int &beta) const;
  void pdist2(const MatrixXd &X, VectorXi &closest) const;

private:

  //Demes: list of vertices (one vertex per row, with x- and y-coordinates)
  //Edges: list of neighbors (one vertex per row, at most six neighbors)
  //Pairs: pairs of connected vertices (one edge per row)
  int nDemes, oDemes, nEdges;
  MatrixXd Demes, Coord;
  MatrixXi Edges, Pairs;
  VectorXi Counts, i2alpha;
  void reindex_demes(const MatrixXd &Demes1, const MatrixXi &Edges1, const MatrixXi &Pairs1,
  		     const VectorXi &newIndex);
  
};
