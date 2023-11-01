#include "exact_solution.h"
#include <deal.II/base/numbers.h>
#include <cmath>


template <int dim>
ExactSolution<dim>::ExactSolution()
{}

template <int dim>
double ExactSolution<dim>::value(const Point<dim> &p, const unsigned int /*component*/) const
{
    const double time = this->get_time();

    #if 0
    // Trapezoid //
    if(p[0]>=0 && p[0]<=1 + time)
        {
            return (1./(1.+time))*p[0]; 
        }
    
    else if(p[0] > (1. + time) && p[0]<(2 +time))
        {
            return 1;
        }

    else if(p[0]>= (2. + time) && p[0] <= 3)
        {
            return -(1./(1.-time))*p[0] + 3*(1./(1.-time)); 
        }

    else {
                return 0;
         }
#endif 
 

    // Block //
   #if 0
    if(p[0] >= 0.25 && p[0]<= 0.75)
    {
        return 1;
    }
    else {
    return 0; 
    }
    #endif 
   
    // Cosine wave //
    return  cos(2 * numbers::PI * p[0] + 2*numbers::PI*time);


}

// Explicit instantiation for the desired dimensions
template class ExactSolution<1>; 
template class ExactSolution<2>;
template class ExactSolution<3>;

