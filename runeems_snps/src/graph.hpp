#pragma once

#include "util.hpp"
#include "habitat.hpp"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

#ifndef GRAPH_H
#define GRAPH_H

/* 
   Triangular graph, possibly irregularly shaped. (The shape is determined by the habitat outline.)
   Or alternatively, the function Graph::generate_grid() can be modified to load a pre-constructed
   population grid, with samples already assigned to vertices.
 */

class Graph {
public:

  Graph( );
  ~Graph( );

  void generate_grid(const string &datapath, const string &gridpath, const Habitat &habitat,
		     const int nDemeDensity, const int nIndiv);
  bool dlmwrite_grid(const string &mcmcpath) const;
  
  bool is_connected() const;
  int get_num_edges() const;
  int get_num_obsrv_demes() const;
  int get_num_total_demes() const;
  int get_deme_of_indiv(const int i) const;
  void get_edge(int edge, int &alpha, int &beta) const;
  void index_closest_to_deme(const MatrixXd &X, VectorXi &Closest) const;
  MatrixXd get_the_obsrv_demes() const;

private:

  // DemeCoord: list of vertices (one vertex per row, with two coordinates)
  // DemePairs: pairs of connected vertices (one edge per row)
  MatrixXd DemeCoord;
  VectorXi DemeSizes, indiv2deme;
  MatrixXi DemePairs;
  
  int neighbors_in_grid(const int r1, const int c1, int &r2, int &c2, const int pos,
			const int nx, const int ny) const;
  void make_triangular_grid(const Habitat &habitat, const int nDemeDensity);
  void map_indiv_to_deme(const string &datapath, const int nIndiv);
  void reindex_demes();

  // The grid to read does not need to be triangular as long as it is connected
  bool read_input_grid(const string &datapath, MatrixXd &DemeCoord, MatrixXi &DemePairs);
  bool read_indiv_to_deme(const string &datapath, const int nDemes, VectorXi &indiv2deme);
  
};

#endif
