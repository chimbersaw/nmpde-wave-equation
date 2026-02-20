#pragma once

#include "config.hpp"

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/numerics/data_out.h>

#include <map>
#include <memory>
#include <string>

struct ConvergenceResult
{
  double       t_final = 0.0;
  double       h = 0.0;
  unsigned int ndofs = 0;
  double       l2_error = 0.0;
  double       h1_error = 0.0;
};

class WaveSolver
{
public:
  using VectorType = dealii::TrilinosWrappers::MPI::Vector;
  using MatrixType = dealii::TrilinosWrappers::SparseMatrix;

  WaveSolver(const WaveProblemConfig &config, MPI_Comm mpi_communicator);

  void
  run();

  void
  set_output_enabled(bool enabled);

  ConvergenceResult
  compute_convergence_result() const;

  const WaveProblemConfig &
  get_config() const;

  const MatrixType &
  get_mass_matrix() const;

  const MatrixType &
  get_stiffness_matrix() const;

  VectorType &
  get_solution();

  VectorType &
  get_velocity();

  VectorType &
  get_acceleration();

  VectorType
  create_owned_vector() const;

  MatrixType
  create_matrix_like_system() const;

  void
  assemble_force(const double time, VectorType &force) const;

  void
  solve_spd_system(const MatrixType &A, VectorType &x, const VectorType &b) const;

  void
  solve_displacement_system(const MatrixType &A,
                            VectorType       &u,
                            const VectorType &rhs,
                            const double      time) const;

  void
  enforce_displacement_bc(VectorType &u, double time) const;

  void
  enforce_velocity_bc(VectorType &v,
                      double      previous_time,
                      double      current_time,
                      double      dt) const;

  void
  output_results(unsigned int step, double time) const;

  bool
  should_output(unsigned int step) const;

private:
  void
  setup_triangulation_and_dofs();

  void
  initialize_state();

  void
  assemble_system_matrices();

  void
  build_constraints(
    double time,
    dealii::AffineConstraints<double> &constraints,
    std::map<dealii::types::global_dof_index, double> &boundary_values) const;

  double
  compute_l2_error(double time) const;

  double
  compute_h1_error(double time) const;

  double
  estimate_mesh_h() const;

  MPI_Comm mpi_communicator;

  WaveProblemConfig config;

  dealii::ConditionalOStream pcout;

  dealii::parallel::fullydistributed::Triangulation<2> triangulation;
  dealii::FE_SimplexP<2>                               fe;
  dealii::DoFHandler<2>                                dof_handler;
  dealii::QGaussSimplex<2>                             quadrature;

  dealii::IndexSet locally_owned_dofs;
  dealii::IndexSet locally_relevant_dofs;

  MatrixType mass_matrix;
  MatrixType stiffness_matrix;

  VectorType solution;
  VectorType velocity;
  VectorType acceleration;

  std::shared_ptr<dealii::Function<2>> u0_function;
  std::shared_ptr<dealii::Function<2>> u1_function;
  std::shared_ptr<dealii::Function<2>> forcing_function;
  std::shared_ptr<dealii::Function<2>> boundary_function;
  std::shared_ptr<dealii::Function<2>> exact_solution_function;

  bool output_enabled = true;
};
