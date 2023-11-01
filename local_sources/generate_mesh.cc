#include "generate_mesh.h"

using namespace dealii; 

template<int dim> 
void generate_mesh(Triangulation<dim> &tria, const unsigned int number_of_refinements, const double left_endpoint, const double right_endpoint) 
{   
    GridGenerator::hyper_cube(tria, left_endpoint, right_endpoint, true); // Mesh is on a square 
    tria.begin_active()->face(0)->set_boundary_id(0); // Label first boundary 
    tria.begin_active()->face(1)->set_boundary_id(1); // Label second boundary 
    tria.refine_global(number_of_refinements); // Refine

}

template void generate_mesh<1>(Triangulation<1> &tria, const unsigned int number_of_refinements, const double left_endpoint, const double right_endpoint); 
template void generate_mesh<2>(Triangulation<2> &tria, const unsigned int number_of_refinements,const double left_endpoint, const double right_endpoint); 
template void generate_mesh<3>(Triangulation<3> &tria, const unsigned int number_of_refinements,const double left_endpoint, const double right_endpoint); 
