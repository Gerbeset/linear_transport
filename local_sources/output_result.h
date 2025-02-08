#pragma once 

#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <iostream>
#include <ostream>

using namespace dealii; 

template<int dim>
void output_result(const DoFHandler<dim> &in_dof_handler, const Vector<double> &result, const unsigned int timestep_number, const double time); 



