
#include "habitat.hpp"

Habitat::Habitat( ) { }
Habitat::~Habitat( ) { }
void Habitat::generate_outer(const string &datapath)
{
  cerr << "[Habitat::initialize]" << endl;
  MatrixXd Outer = readMatrixXd(datapath + ".outer");
  if (Outer.rows()==0||Outer.cols()!=2) {
    cerr << "  Error reading habitat points from " << datapath + ".outer" << endl
	 << "  Expect a list of points, with two coordinates per line" << endl; exit(1);     
  }
  cerr << "  Loaded habitat points from " << datapath + ".outer" << endl;
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
  cerr << "  Input habitat: " << endl << "    " << boost::geometry::wkt<Ring>(domain) << endl;
  if (!boost::geometry::is_valid(domain)) {
    boost::geometry::correct(domain);
    cerr << "  The habitat is not a valid ring (a simple closed polygon)" << endl;
    cerr << "  Corrected habitat: " << endl << "    " << boost::geometry::wkt<Ring>(domain) << endl;
  }
  double area = this->get_area( );
  if (area<=0)
    { cerr << "  The habitat has area " << area << ". Specify a valid ring habitat" << endl; exit(1); }
  cerr << "[Habitat::initialize] Done." << endl << endl; 
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
bool Habitat::dlmwrite_outer(const string &mcmcpath) const {
  ofstream out((mcmcpath + "/outer.txt").c_str(),ofstream::out);
  if (!out.is_open()) { return false; }
  vector<Point>::const_iterator it;
  for ( it = domain.begin(); it != domain.end() ; it++ ) {
    out << boost::geometry::get<0>(*it) << " "
	<< boost::geometry::get<1>(*it) << endl;
  }
  out.close();
  return true;
}
