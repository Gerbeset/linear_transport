#include "vectors.h"

vectors_for_problem::vectors_for_problem()
{
}

void vectors_for_problem::initialize(const unsigned int &n)
{
    solution.reinit(n);
    old_solution.reinit(n); 
    rhs_data.reinit(n);
    initial_condition.reinit(n); 
    final_solution.reinit(n); 
}