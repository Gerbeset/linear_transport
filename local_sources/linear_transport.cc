#include <__functional/unary_negate.h>
#include <array>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/patterns.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <fstream>
#include <iostream>
#include <ostream>
#include <string>
#include <vector>
#include <cmath>
#include <deal.II/base/parameter_handler.h>
#include "flux.h"
#include "generate_mesh.h"
#include "compute_mesh_size.h"
#include "compute_time_step.h"
#include "interpolate_exact_at_time.h"
#include "error_analysis.h"
#include "vectors.h"
#include "error_analysis.h"
#include "explicit.h"


using namespace dealii;


int main()
{
  // Space dimension of the problem
  const unsigned int dim = 1; 

  // Parameter handling
  ParameterHandler prm; 
  prm.declare_entry("number_of_refinements", "1", Patterns::Integer(0));
  prm.declare_entry("CFL", "1.0", Patterns::Double(0.0000000001)); 
  prm.declare_entry("final_time", "0.01", Patterns::Double(0.0));

  prm.parse_input("parameter.prm", "", true); 

  const unsigned int refinement = prm.get_integer("number_of_refinements"); 
  const double final_time = prm.get_double("final_time"); 
  const double CFL = prm.get_double("CFL"); 

  // End of parameters 

  FE_Q<dim> fe(1); 

  Triangulation<dim> triangulation;

  generate_mesh(triangulation , refinement, 0, 1);

  DoFHandler<dim> dof_handler(triangulation); 
  
  dof_handler.distribute_dofs(fe); 

  std::cout<<"Number of dofs:"<<dof_handler.n_dofs()<<std::endl; 

  vectors_for_problem vectors; 
  vectors.initialize(dof_handler.n_dofs()); 

  interpolate_exact_at_time(dof_handler, vectors.initial_condition, 0.0); 

  DataOut<1> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(vectors.initial_condition, "initial_condition");
  data_out.build_patches();
  std::ofstream output("initial_condition.vtu");
  data_out.write_vtu(output);

  double mesh_size = compute_mesh_size(triangulation); 

  std::cout<<"mesh size = "<<mesh_size<<std::endl<<std::endl; 

  double max_initial_value; 

 for(unsigned int i = 0; i < vectors.initial_condition.size(); ++i)
  {
    max_initial_value = std::max(max_initial_value, abs(vectors.initial_condition[i])); 
  }

  std::cout<<"max_initial_value = " << max_initial_value<<std::endl; 

  explicit_time_step time_step_data(CFL, mesh_size, final_time, max_initial_value); 
  double time_step = time_step_data.time_step; 
  std::cout<<"time_step = "<<time_step<<std::endl;

  unsigned int number_of_time_steps = time_step_data.number_of_time_steps; 
  std::cout<<"number of time steps = "<<number_of_time_steps<<std::endl; 

  explicit_step explicit_(mesh_size); 
  
  double time =0.0; 
  unsigned int time_step_number = 0; 

  vectors.old_solution = vectors.initial_condition; 

  while(time_step_number<number_of_time_steps)
  {
    time += time_step; 
    ++ time_step_number; 

    if(time_step_number % 100 == 0.0)
    {
      std::cout<<"Current time step: "<< time_step_number<<std::endl; 
    }
    
    vectors.solution = explicit_.step(vectors.old_solution, time_step); 
    vectors.old_solution = vectors.solution; 
  }

  std::cout<<"Final time = "<<time<<std::endl; 

  vectors.final_solution = vectors.old_solution; 

  DataOut<1> data_out_1;
  data_out_1.attach_dof_handler(dof_handler);
  data_out_1.add_data_vector(vectors.final_solution, "final_solution");
  data_out_1.build_patches();
  std::ofstream output_1("final_solution.vtu");
  data_out_1.write_vtu(output_1);

  // std::cout<<"Final solution = "<<vectors.final_solution<<std::endl; 

  compute_error(triangulation, vectors.final_solution, time, dof_handler, fe); 
  }


