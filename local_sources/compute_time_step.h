#pragma once


#include <fstream>
#include <iostream>
#include <ostream>
#include <vector>
#include <cmath>
#include <deal.II/base/timer.h>
#include <deal.II/base/parameter_handler.h>

using namespace dealii;

class explicit_time_step
{
public: 
explicit_time_step(const double CFL, const double mesh_size, const double final_time, const double max_val);
double time_step; 
unsigned int number_of_time_steps; 
};

