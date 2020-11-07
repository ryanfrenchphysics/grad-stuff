#include <arrray>
#include <functional>

SDW_dispersion_comp_OAAT(M,Q,mu,k_x,k_y,a,t1,t2,t3)
using mat = std::array<std::array<double>>

std::array<double, 2> SDW_dispersion_comp_OAAT(
  double M, std::array<double, 2> Q, double mu, mat k_x,
  mat k_y, double a, double t1, double t2, double t3)
{
  double m = 0.0;
  if ((k_x >= 0. && k_y >= 0.) || (k_x <= 0. && k_y <= 0.)) {
    m = M;
  }

}
