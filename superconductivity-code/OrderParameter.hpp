#ifndef _ORDER_PARAMETER_HPP
#define _ORDER_PARAMETER_HPP

//#include <cmath>    // Maybe?
#include <vector>
//#include <math>

using vec = std::vector<float>;


class OrderParameter
{
public:
  OrderParameter(int num_op_components) :
    op_mat(num_op_components + 1, vec(1, 0.0));
  ~OrderParameter();
  std::vector<vec> get_mat();
  void push_vals(vec const &in_vec);

// Example Δhat(pf) = [Δ sin(2*ϕ) * iσ₂]
// Δ(R), Y = [sin(2ϕ), 0, 0, 0] (singlet basis fxn)
protected:

private:
  int num_op_components;  // Number of order parameter components

  // op_array = some C++ array of num_op_components + 1
  std::vector<vec> op_mat;

  // -> op_array = new complex<double>[num_op_components]
};

/*
Example use:
OrderParameter op(2);

valarr = array(complex, len=2)
valarr = op.get_array();

*/



#endif //_ORDER_PARAMETER_HPP
