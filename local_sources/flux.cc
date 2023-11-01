#include "flux.h"
#include "d_ij.h"
#include <cmath>

using namespace dealii; 

double continuous_flux(const double u)
{
        return u; 
     //return 0.5*u*u; 
}

Vector<double> numerical_flux(const Vector<double> &U)
{
    double c_iim1 = -0.5; 
    double c_iip1 = 0.5; 

    Vector<double> U_new(U.size()); 

    U_new[0] = (continuous_flux(U[U.size()-2])*c_iim1 + continuous_flux(U[1])*c_iip1) - d_ij(U[U.size()-2], U[0])*(U[U.size() -2] - U[0]) - d_ij(U[0], U[1])*(U[1] - U[0]); 

    for(unsigned int i = 1; i<U.size()-1; ++i)
    {
        U_new[i] = (continuous_flux(U[i-1])*c_iim1 + continuous_flux(U[i+1])*c_iip1) - d_ij(U[i],U[i-1])*(U[i-1] - U[i]) - d_ij(U[i],U[i+1])*(U[i+1] - U[i]); 
    }
    U_new[U.size()-1] = U_new[0]; 

   // std::cout<<"numerical flux = "<<U_new<<std::endl; 

    return U_new;
}
