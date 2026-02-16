#ifndef WAVE_EQUATION_HPP
#define WAVE_EQUATION_HPP

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/vector.h>

#include <map>
#include <string>

namespace WaveEquation
{
  template <int dim>
  class WaveEquationProblem
  {
  public:
    WaveEquationProblem(unsigned int       global_refinements,
                        double             final_time,
                        double             cfl_number,
                        const std::string &mesh_path);

    void
    run();

  private:
    class ExactSolution : public dealii::Function<dim>
    {
    public:
      explicit ExactSolution(double time = 0.0);

      ~ExactSolution() override = default;

      double
      value(const dealii::Point<dim> &point, unsigned int component = 0) const override;
    };

    class InitialVelocity : public dealii::Function<dim>
    {
    public:
      InitialVelocity();

      ~InitialVelocity() override = default;

      double
      value(const dealii::Point<dim> &point, unsigned int component = 0) const override;
    };

    class RightHandSide : public dealii::Function<dim>
    {
    public:
      explicit RightHandSide(double time = 0.0);

      ~RightHandSide() override = default;

      double
      value(const dealii::Point<dim> &point, unsigned int component = 0) const override;
    };

    void
    make_grid();

    void
    setup_system();

    void
    assemble_matrices();

    void
    compute_time_step();

    void
    initialize_solution();

    void
    assemble_forcing_term(double time, dealii::Vector<double> &forcing) const;

    std::map<dealii::types::global_dof_index, double>
    get_boundary_values(double time) const;

    void
    solve_mass_system(const std::map<dealii::types::global_dof_index, double> &boundary_values,
                      dealii::Vector<double>                                  &destination,
                      const dealii::Vector<double>                            &rhs) const;

    void
    output_results(unsigned int timestep) const;

    void
    compute_error() const;

    const unsigned int global_refinements;
    const double       final_time;
    const double       cfl_number;
    const std::string  mesh_path;
    const unsigned int minimum_time_steps;
    const std::string  output_directory;

    double       time;
    double       time_step;
    unsigned int timestep_number;
    unsigned int output_every;

    dealii::Triangulation<dim> triangulation;
    dealii::FE_Q<dim>          fe;
    dealii::DoFHandler<dim>    dof_handler;

    dealii::AffineConstraints<double> constraints;

    dealii::SparsityPattern      sparsity_pattern;
    dealii::SparseMatrix<double> mass_matrix;
    dealii::SparseMatrix<double> stiffness_matrix;

    dealii::Vector<double> solution_old;
    dealii::Vector<double> solution;
    dealii::Vector<double> solution_new;
    dealii::Vector<double> system_rhs;
    dealii::Vector<double> forcing_vector;
  };
} // namespace WaveEquation

#endif // WAVE_EQUATION_HPP
