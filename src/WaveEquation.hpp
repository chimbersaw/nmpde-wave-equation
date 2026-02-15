#ifndef WAVE_EQUATION_HPP
#define WAVE_EQUATION_HPP

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h> 

#include <fstream>
#include <iostream>
#include <cmath>

using namespace dealii;

// Classe per la condizione iniziale: una gaussiana centrata nel dominio
template <int dim>
class InitialCondition : public Function<dim> {
public:
  InitialCondition() : Function<dim>() {}
  virtual double value(const Point<dim> &p, const unsigned int component = 0) const override {
    (void)component;
    const double x = p[0] - 0.5;
    const double y = p[1] - 0.5;
    const double sigma = 0.1;
    return std::exp(-(x * x + y * y) / (sigma * sigma));
  }
};

template <int dim>
class WaveEquation {
public:
  WaveEquation(const std::string &mesh_file_name_, 
               const unsigned int &r_,
               const double &T_, 
               const double &delta_t_);

  void run();

private:
  void setup();
  void assemble_system();
  void solve_timestep();
  void output_results() const;

  const std::string mesh_file_name;
  const unsigned int r;
  const double T;
  const double delta_t;

  double time;
  unsigned int timestep_number;

  const unsigned int mpi_size;
  const unsigned int mpi_rank;
  parallel::fullydistributed::Triangulation<dim> mesh;
  DoFHandler<dim> dof_handler;
  std::unique_ptr<FiniteElement<dim>> fe;
  std::unique_ptr<Quadrature<dim>> quadrature;

  TrilinosWrappers::SparseMatrix system_matrix;
  TrilinosWrappers::MPI::Vector system_rhs;
  
  TrilinosWrappers::MPI::Vector solution;       
  TrilinosWrappers::MPI::Vector old_solution;   
  TrilinosWrappers::MPI::Vector older_solution; 
  TrilinosWrappers::MPI::Vector solution_owned;

  ConditionalOStream pcout;
};

#endif