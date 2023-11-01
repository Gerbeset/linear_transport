#include "explicit.h"
#include "flux.h"
#include <deal.II/lac/vector.h>

explicit_step::explicit_step(const double mesh_size) : m_i(mesh_size)
{}

Vector<double> explicit_step::step(const Vector<double> &U_old, const double tau)
{
    Vector<double> U_new(U_old.size());
    Vector<double> F_of_U(U_old.size());  

    //std::cout<<"In explicit: tau = "<<tau<<std::endl<<std::endl;

    F_of_U = numerical_flux(U_old); 
    //std::cout<<"Numercial flux = "<<F_of_U<<std::endl<<std::endl; 
    U_new = U_old;
    //std::cout<<"U_new before flux = "<<U_new<<std::endl<<std::endl; 
    //std::cout<<"-tau/m_i = "<<-tau/m_i<<std::endl<<std::endl; 
    U_new.add(-tau/m_i, F_of_U);  

    //std::cout<<"New U_new = "<<U_new<<std::endl<<std::endl; 

    return U_new; 
}


