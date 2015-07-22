#pragma once

#include "util.hpp"

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/ring.hpp>
#include <boost/geometry/algorithms/covered_by.hpp>

#ifndef HABITAT_H
#define HABITAT_H

typedef boost::geometry::model::point<double,2,boost::geometry::cs::cartesian> Point;
typedef boost::geometry::model::ring<Point> Ring;

/*
  The habitat is represented as a ring (a polygon without holes).
  The vertices of the ring are specified by the user in the file `datapath.outer`.
  For example, suppose that `datapath.outer` contains the following six lines:
      0 0
      3 3
      0 6
      12 6
      12 0
      0 0
  Then the habitat is given by (0,0) - (3,3) - (0,6) - (12,6) - (12,0) - (0,0).
  The vertices of the population graph (the demes) will fall inside the habitat by construction.
  
  However, EEMS will not check that the sampling locations fall inside the habitat -- instead,
  each sample will be assigned to the closest deme. Therefore, the user should specify a habitat
  that is sufficiently large to cover all sampling locations.
*/

class Habitat {
public:

  Habitat( );
  ~Habitat( );

  void generate_outer(const string &datapath);
  bool dlmwrite_outer(const string &mcmcpath) const;
  bool in_point(const double x, const double y) const;

  double get_area( ) const;
  double get_xmin( ) const;
  double get_xmax( ) const;
  double get_ymin( ) const;
  double get_ymax( ) const;
  double get_xspan( ) const; 
  double get_yspan( ) const;
 
private:

  Ring domain;
  double xmin, xmax, xspan;
  double ymin, ymax, yspan;
  
};

#endif
