
#include "habitat.hpp"

Habitat::Habitat( ) { }
Habitat::~Habitat( ) { }
void Habitat::initialize(const string &datapath)
{
  cerr << "[Habitat::initialize]" << endl;
  string infile = datapath + ".outer";
  ifstream instrm(infile.c_str(), ios::in);
  if(!instrm.is_open( ))
    { cerr << "  Error opening habitat points " << infile << endl; exit(1); }
  vector<Point> points; points.clear(); string line;
  xmin = 0.0; xmax = 0.0;
  ymin = 0.0; ymax = 0.0;
  while (getline(instrm,line))
    {
      vector<double> pt = split(line);
      if (pt.size()!=2)	{
	cerr << "  Error reading habitat points " << infile << endl
	     << "  Expect a list of points, with two coordinates per line" << endl;  exit(1);
      }
      if (points.size()==1) {
	xmin = pt[0]; xmax = pt[0];
	ymin = pt[1]; ymax = pt[1];
      } else {
	if (pt[0]<xmin) { xmin = pt[0]; }
	if (pt[0]>xmax) { xmax = pt[0]; }
	if (pt[1]<ymin) { ymin = pt[1]; }
	if (pt[1]>ymax) { ymax = pt[1]; }
      }
      points.push_back(Point(pt[0],pt[1]));
      pt.clear( );
    }
  instrm.close( );
  cerr << "  Read habitat points from " << infile << endl;
  xspan = xmax - xmin;
  yspan = ymax - ymin;
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
bool Habitat::dlmwrite(const string &mcmcpath) const {
  ofstream out;
  out.open((mcmcpath + "/outer.txt").c_str(),ofstream::out);
  if (!out.is_open()) { return false; }
  vector<Point>::const_iterator it;
  for ( it = domain.begin(); it != domain.end() ; it++ ) {
    out << boost::geometry::get<0>(*it) << " " << boost::geometry::get<1>(*it) << endl;
  }
  out.close();
  return true;
}
