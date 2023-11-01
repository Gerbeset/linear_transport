#pragma once
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/block_vector_base.h>
#include <deal.II/numerics/vector_tools_common.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/vector_tools_common.h>
#include <deal.II/lac/vector.h>
#include <fstream>
#include <iostream>
#include <ostream>
#include <vector>
#include <cmath>
#include "d_ij.h"
#include "flux.h"

class explicit_step
{
    public: 
    explicit_step(const double mesh_size); 
    Vector<double> step(const Vector<double> &U_old, const double tau);

    private: 
    double m_i; 
};

