
#include "habitat.hpp"

Habitat::Habitat( ) { }
Habitat::~Habitat( ) { }
void Habitat::generate_outer(const string &datapath, const string &mcmcpath)
{
  cout << "[Habitat::initialize]" << endl;
  MatrixXd Outer = readMatrixXd(datapath + ".outer");
  // A valid ring is formed by at least 3 points
  check_condition(Outer.rows() >= 3 && Outer.cols() == 2,
		  "Check that " + datapath + ".outer is a list of locations, two coordinates per line.");
  cout << "  Loaded habitat points from " << datapath + ".outer" << endl;
  xmin = Outer.col(0).minCoeff();
  xmax = Outer.col(0).maxCoeff();
  ymin = Outer.col(1).minCoeff();
  ymax = Outer.col(1).maxCoeff();
  xspan = xmax - xmin;
  yspan = ymax - ymin;
  vector<Point> points;
  for ( int i=0 ; i<Outer.rows() ; i++ ) {
    points.push_back(Point(Outer(i,0),Outer(i,1)));
  }
  boost::geometry::append(domain,points);
  // Correction is limited to fixing the orientation and a closing point (if necessary)
  boost::geometry::correct(domain);
  cout << "  Input habitat: " << endl << "    " << boost::geometry::wkt<Ring>(domain) << endl;
  dlmwrite_outer(mcmcpath);
  check_condition(boost::geometry::is_valid(domain) && boost::geometry::area(domain) > 0,
		  "Specify a valid ring habitat (a simple closed polygon).");
  cout << "[Habitat::initialize] Done." << endl << endl; 
}
double Habitat::get_area( ) const { return (boost::geometry::area(domain)); }
double Habitat::get_xmin( ) const { return (xmin); }
double Habitat::get_xmax( ) const { return (xmax); }
double Habitat::get_ymin( ) const { return (ymin); }
double Habitat::get_ymax( ) const { return (ymax); }
double Habitat::get_xspan( ) const { return (xspan); }
double Habitat::get_yspan( ) const { return (yspan); }
bool Habitat::in_point(const double x, const double y) const
{
  return (boost::geometry::covered_by(Point(x,y),domain));
}
void Habitat::dlmwrite_outer(const string &mcmcpath) const {
  ofstream out((mcmcpath + "/outer.txt").c_str(),ofstream::out);
  check_condition(out.is_open(), "Cannot open " + mcmcpath + "/outer.txt for writing");
  vector<Point>::const_iterator it;
  for ( it = domain.begin(); it != domain.end() ; it++ ) {
    out << boost::geometry::get<0>(*it) << " "
	<< boost::geometry::get<1>(*it) << endl;
  }
  out.close();
}
