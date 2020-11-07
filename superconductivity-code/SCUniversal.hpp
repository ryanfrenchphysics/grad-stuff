#ifndef _SCERRORS_HPP
#define _SCERRORS_HPP
#include <string>
#include <iostream>

/*
 *  Error codes
 */
// Grid shape not defined
const std::string GRID_SHAPE_NOT_DEF = "Grid shape is not defined.";
// Grid redefinition error
const std::string GRID_REDEF_ERROR = "Grid was already defined, cannot redefine.";

// Grid rectangular limits error (min > max)
const std::string GRID_RECT_LIM_ERROR = "Grid minimum limit greater than maximum limit.";


/*
 *  Constant defintions
 */
const auto GRID_DENSITY_FINE = 6.0;
const auto GRID_DENSITY_MEDIUM = 4.0;
const auto GRID_DENSITY_ROUGH = 2.0;

const auto GRID_SHAPE_RECT = 1;

// How close to boundary is considered boundary?
const auto GRID_BOUNDARY_TOL = 0.01;

// Boundary types
const auto SPECULAR_BOUNDARY = 1;



// namespace scerror
// {

void logSCerror(auto error_type)
{
  std::cerr << "Error: " << error_type << std::endl;
}

//} // end namespace scerror

#endif // _SCERRORS_HPP
