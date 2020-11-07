#include "OrderParameter.hpp"


/*
 *  OrderParameter Constructor:
 *    args: # of order parameter components
 *    returns: N/A
 *
 *  Reserves space for the # of order parameter components
 */
OrderParameter::OrderParameter(int num_op_components) :
  op_mat(num_op_components + 1, vec(1, 0.0))
{
  this->num_op_components = num_op_components;
}


/*
 *  OrderParameter Destructor:
 *    args: N/A
 *    returns: N/A
 */
OrderParameter::~OrderParameter()
{

}


/*
 *  get_mat:
 *    args: N/A
 *    returns: order parameter matrix
 *
 *  Returns the OP vector
 */
vec<vec<float>> OrderParameter::get_mat()
{
  return op_mat;
}


/*
 *  push_vals:
 *    args: vector of float values
 *    returns: N/A
 *
 *  Push a vector to the existing matrix (updated values of the OP for each
 *  component)
 */
void OrderParameter::push_vals(vec const &in_vec)
{
  this->op_mat.push_back(in_vec);
  return;
}
