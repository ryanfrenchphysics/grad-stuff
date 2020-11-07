#define WITHOUT_NUMPY

#include "Grid.hpp"
#include "SCUniversal.hpp"
#include "matplotlibcpp.hpp"
#include <cmath>
#include <vector>
#include <iostream>

namespace plt = matplotlibcpp;

/*
 *
 */
class Grid::Point
{
public:
  // Constructor
  Point(double x, double y, double z)
  {
    ptvec = {x, y, z};
  }
  // Destructor
  ~Point()
  {

  }

  auto get_pt()
  {
    return ptvec;
  }

  void is_boundary(bool tf)
  {
    boundary = tf;
  }

  void setBC(int bctype)
  {
    bc_type = bctype;
  }

  int getBC()
  {
    return bc_type;
  }

private:
  int bc_type;
  bool boundary;
  std::vector<double> ptvec = {0.0, 0.0, 0.0};

};




/*
 *
 */
Grid::Grid(int dim) : dimensions(dim)
{

}


/*
 *
 */
Grid::~Grid()
{
  delete point_vec;
}


/*
 *
 */
void Grid::set_density(double density_val)
{
  if (shape_type == GRID_SHAPE_RECT) {
    delx = rect_xrange / 100; //* exp(-density_val);
    dely = rect_yrange / 100; //* exp(-density_val);
  }
}


/*
 *
 */
void Grid::set_separations(double x_sep, double y_sep, double z_sep)
{
  delx = x_sep;
  dely = y_sep;
  delz = z_sep;
}


/*
 *
 */
 void Grid::rectangle(double xmin, double xmax, double ymin, double ymax)
{
  if (shape_determined == false) {
    shape_determined = true;
    shape_type = GRID_SHAPE_RECT;
    if (xmin > xmax || ymin > ymax) {
      logSCerror(GRID_RECT_LIM_ERROR);
      return;
    }
  } else {
    logSCerror(GRID_REDEF_ERROR);
    return;
  }
  x_min = xmin;
  x_max = xmax;
  y_min = ymin;
  y_max = ymax;
  rect_xrange = xmax - xmin;
  rect_yrange = ymax - ymin;

}


/*
 * Temporarily: All boundaries are specular
 */
void Grid::gen_grid()
{
  if (shape_determined == false) {
    logSCerror(GRID_SHAPE_NOT_DEF);
  }
  long num_xpoints = (long) (rect_xrange) / (delx);
  long num_ypoints = (long) (rect_yrange) / (dely);
  long num_points = (long) num_xpoints * num_ypoints;

  Grid::Point nullpt = create_point(false, 1, 0.0, 0.0, 0.0);
  point_vec = new std::vector<Point>(num_points, nullpt);
  std::cout << "vec size: " << point_vec->size() << "\n";
  std::cout << "xpts : " << num_xpoints << "\n";
  std::cout << "ypts: " << num_ypoints << "\n";


  auto x = 0;
  auto x_val = 0.0;
  auto y_val = 0.0;
  auto z_val = 0.0;
  for (auto i=0; i < num_xpoints; i++) {
    x_val = x_min + (i * delx);
    for (auto j=0; j < num_ypoints; j++) {
      std::cout << "x: " << x <<"\n";
      y_val = y_min + (j * dely);
      //auto z_val = z_min + (i * delz);

      auto bound = is_boundary(x_val, y_val, z_val);
      point_vec->at(x) = create_point(bound, SPECULAR_BOUNDARY, x_val, y_val, z_val);
      x++;
    }
  }
  std::cout << "Finished grid";

  for (unsigned long i=0; i < point_vec->size(); i++) {
    for (int j=0; j<2; j++) {
      std::cout << point_vec->at(i).get_pt()[j] << " ";
    }
    std::cout << "\n";
  }
}


/*
 *
 */
void Grid::plot_grid()
{
  //std::vector<double> xvals(point_vec->size(), 0.0);
  //std::vector<double> yvals(point_vec->size(), 0.0);
  // for (unsigned long i=0; i < point_vec->size(); i++) {
  //   xvals.at(i) = (point_vec->at(i).get_pt()[0]);
  //   yvals.at(i) = (point_vec->at(i).get_pt()[1]);
  // }
  //
  // for (unsigned long i=0; i < xvals.size(); i++) {
  //   std::cout << "(" << i << "): " << xvals.at(i) << " " << yvals.at(i) << "\n";
  // }
  std::cout<<"here";
  double x[3] = {0, 1, 2};
  double y[3] = {0, 1, 2};

  std::vector<double> xs(3, 0.0);
  std::vector<double> ys(3, 0.0);

  for (int i=0; i<3; i++) {
    xs.at(i) = x[i];
    ys.at(i) = y[i];
  }

  plt::plot({1, 2, 3});
  plt::show();
}


/*
 *
 */
void Grid::print_points()
{
  for (unsigned long i=0; i < point_vec->size(); i++) {
    for (auto j=0; j < 3; j++) {
      std::cout << point_vec->at(i).get_pt()[j] << " ";
    }
    std::cout << std::endl;
  }
  return;
}


/*
 *
 */
Grid::Point Grid::create_point(bool is_boundary, int bc_type, double x, double y, double z)
{
  Point p(x, y, z);
  p.is_boundary(is_boundary);
  p.setBC(bc_type);

  // for (int i=0; i<3; i++) {
  //   std::cout << p.get_pt()[i] << " ";
  // }
  // std::cout << "\n";
  return p;
}


/*
 *
 */
bool Grid::is_boundary(auto xpt, auto ypt, auto zpt)
{
  if (
    (abs(xpt - x_min) <= GRID_BOUNDARY_TOL ||
    abs(xpt - x_max) <= GRID_BOUNDARY_TOL) &&
    (abs(ypt - y_min) <= GRID_BOUNDARY_TOL ||
    abs(ypt - y_max) <= GRID_BOUNDARY_TOL) &&
    (abs(zpt - z_min) <= GRID_BOUNDARY_TOL ||
    abs(zpt - z_max) <= GRID_BOUNDARY_TOL)
  ) {
    return true;
  }
  return false;
}
