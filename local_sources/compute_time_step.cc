#include "compute_time_step.h"
using namespace dealii;

explicit_time_step::explicit_time_step(const double CFL, const double m_i, const double final_time, const double max_val)
{
  time_step = CFL*m_i/(4.*max_val); 

  double temp_val = final_time/time_step; 
  
  number_of_time_steps = static_cast<int>(temp_val); 
  time_step = final_time/number_of_time_steps; 
}