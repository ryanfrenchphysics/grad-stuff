#ifndef _GRID_H
#define _GRID_H

#include <vector>

using vecf = std::vector<double>;
using vecmatf = std::vector<vecf>;

class Grid
{
public:
  Grid(int dim);
  ~Grid();

  void set_density(double density_val);
  void set_separations(double x_sep, double y_sep, double z_sep);
  void rectangle(double xmin, double xmax, double ymin, double ymax);

  void gen_grid();


  void plot_grid();
  void print_points();

private:
  int dimensions;
  bool shape_determined = false;
  // Use shape type constants
  int shape_type = 0;
  double delx = 1.0;
  double dely = 1.0;
  double delz = 1.0;
  double x_min = 0.0;
  double x_max = 0.0;
  double y_min = 0.0;
  double y_max = 0.0;
  double z_min = 0.0;
  double z_max = 0.0;
  double rect_xrange = 0.0;
  double rect_yrange = 0.0;
  double rect_zrange = 0.0;

  class Point;
  std::vector<Point> * point_vec;

  Point create_point(bool is_boundary, int bc_type, double x, double y, double z);

  bool is_boundary(auto xpt, auto ypt, auto zpt);

};


#endif // _GRID_H
