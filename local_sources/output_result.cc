#include "output_result.h"

template<int dim>
void output_result(const DoFHandler<dim> &in_dof_handler, const Vector<double> &result, const unsigned int timestep_number, const double time)
{
  DataOut<dim> data_out; 
  data_out.attach_dof_handler(in_dof_handler);
  data_out.add_data_vector(result, "u");
  data_out.build_patches();
  data_out.set_flags(DataOutBase::VtkFlags(time, timestep_number));
  const std::string filename =
    "solution-" + Utilities::int_to_string(timestep_number, 3) + ".vtu";
  std::ofstream output(filename);
  data_out.write_vtu(output);
}

template void output_result<1>(const DoFHandler<1> &in_dof_handler, const Vector<double> &result, const unsigned int timestep_number, const double time);
template void output_result<2>(const DoFHandler<2> &in_dof_handler, const Vector<double> &result, const unsigned int timestep_number, const double time);
template void output_result<3>(const DoFHandler<3> &in_dof_handler, const Vector<double> &result, const unsigned int timestep_number, const double time);

