#pragma once
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/block_vector_base.h>
#include <deal.II/numerics/vector_tools_common.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/vector_tools_common.h>
#include <fstream>
#include <iostream>
#include <ostream>
#include <vector>
#include <cmath>
#include "d_ij.h"

using namespace dealii; 

double continuous_flux(const double u); 

Vector<double> numerical_flux(const Vector<double> &U_old); 


