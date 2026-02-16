#include "WaveEquation.hpp"

#include <deal.II/base/exceptions.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_direct.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>

namespace WaveEquation
{
  template <int dim>
  WaveEquationProblem<dim>::ExactSolution::ExactSolution(const double time)
    : dealii::Function<dim>(1, time)
  {}

  template <int dim>
  double
  WaveEquationProblem<dim>::ExactSolution::value(const dealii::Point<dim> &point, const unsigned int component) const
  {
    (void)component;
    Assert(dim == 2, dealii::ExcNotImplemented());

    double pi    = dealii::numbers::PI;
    double omega = std::sqrt(2.0) * pi;
    return std::sin(pi * point[0]) * std::sin(pi * point[1]) * std::cos(omega * this->get_time());
  }

  template <int dim>
  WaveEquationProblem<dim>::InitialVelocity::InitialVelocity()
    : dealii::Function<dim>(1, 0.0)
  {}

  template <int dim>
  double
  WaveEquationProblem<dim>::InitialVelocity::value(const dealii::Point<dim> &point, const unsigned int component) const
  {
    (void)point;
    (void)component;
    return 0.0;
  }

  template <int dim>
  WaveEquationProblem<dim>::RightHandSide::RightHandSide(const double time)
    : dealii::Function<dim>(1, time)
  {}

  template <int dim>
  double
  WaveEquationProblem<dim>::RightHandSide::value(const dealii::Point<dim> &point, const unsigned int component) const
  {
    (void)point;
    (void)component;
    return 0.0;
  }

  template <int dim>
  WaveEquationProblem<dim>::WaveEquationProblem(const unsigned int global_refinements,
                                                const double       final_time,
                                                const double       cfl_number,
                                                const std::string &mesh_path)
    : global_refinements(global_refinements)
    , final_time(final_time)
    , cfl_number(cfl_number)
    , mesh_path(mesh_path)
    , minimum_time_steps(500)
    , output_directory(std::filesystem::exists("src") ? "solution" : "../solution")
    , time(0.0)
    , time_step(0.0)
    , timestep_number(0)
    , output_every(1)
    , fe(1)
    , dof_handler(triangulation)
  {}

  template <int dim>
  void
  WaveEquationProblem<dim>::make_grid()
  {
    std::filesystem::path resolved_mesh_path;
    if (!mesh_path.empty())
      {
        const std::filesystem::path user_mesh_path(mesh_path);
        if (std::filesystem::exists(user_mesh_path))
          {
            resolved_mesh_path = user_mesh_path;
          }
        else
          {
            const std::filesystem::path fallback_path = std::filesystem::path("..") / user_mesh_path;
            if (std::filesystem::exists(fallback_path))
              {
                resolved_mesh_path = fallback_path;
              }
          }
      }

    if (!resolved_mesh_path.empty())
      {
        dealii::GridIn<dim> grid_in;
        grid_in.attach_triangulation(triangulation);

        std::ifstream input(resolved_mesh_path);
        AssertThrow(input.is_open(), dealii::ExcMessage("Could not open mesh file: " + resolved_mesh_path.string()));

        grid_in.read_msh(input);
        std::cout << "Loaded mesh from " << resolved_mesh_path.string() << '\n';
      }
    else
      {
        dealii::GridGenerator::hyper_cube(triangulation, 0.0, 1.0);
        triangulation.refine_global(global_refinements);
        std::cout << "Mesh file not found, generated unit square mesh with " << triangulation.n_active_cells()
                  << " cells.\n";
      }
  }

  template <int dim>
  void
  WaveEquationProblem<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe);

    constraints.clear();
    dealii::DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    constraints.close();

    dealii::DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    dealii::DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
    sparsity_pattern.copy_from(dsp);

    mass_matrix.reinit(sparsity_pattern);
    stiffness_matrix.reinit(sparsity_pattern);

    solution_old.reinit(dof_handler.n_dofs());
    solution.reinit(dof_handler.n_dofs());
    solution_new.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
    forcing_vector.reinit(dof_handler.n_dofs());

    std::cout << "DoFs: " << dof_handler.n_dofs() << '\n';
  }

  template <int dim>
  void
  WaveEquationProblem<dim>::assemble_matrices()
  {
    dealii::QGauss<dim>   quadrature_formula(fe.degree + 1);
    const auto            flags = dealii::update_values | dealii::update_gradients | dealii::update_JxW_values;
    dealii::FEValues<dim> fe_values(fe, quadrature_formula, flags);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    dealii::FullMatrix<double>                   cell_mass(dofs_per_cell, dofs_per_cell);
    dealii::FullMatrix<double>                   cell_stiffness(dofs_per_cell, dofs_per_cell);
    std::vector<dealii::types::global_dof_index> local_dof_indices(dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        fe_values.reinit(cell);
        cell_mass      = 0.0;
        cell_stiffness = 0.0;

        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                const double                 phi_i      = fe_values.shape_value(i, q);
                const dealii::Tensor<1, dim> grad_phi_i = fe_values.shape_grad(i, q);
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                  {
                    const double                 phi_j      = fe_values.shape_value(j, q);
                    const dealii::Tensor<1, dim> grad_phi_j = fe_values.shape_grad(j, q);

                    cell_mass(i, j) += phi_i * phi_j * fe_values.JxW(q);
                    cell_stiffness(i, j) += grad_phi_i * grad_phi_j * fe_values.JxW(q);
                  }
              }
          }

        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(cell_mass, local_dof_indices, mass_matrix);
        constraints.distribute_local_to_global(cell_stiffness, local_dof_indices, stiffness_matrix);
      }
  }

  template <int dim>
  void
  WaveEquationProblem<dim>::compute_time_step()
  {
    const double h_min  = dealii::GridTools::minimal_cell_diameter(triangulation);
    const double safety = std::max(1.0, std::sqrt(static_cast<double>(dim)));
    const double raw_dt = cfl_number * h_min / safety;

    const unsigned int cfl_steps = std::max(1U, static_cast<unsigned int>(std::ceil(final_time / raw_dt)));
    const unsigned int n_steps   = std::max(cfl_steps, minimum_time_steps);
    time_step                    = final_time / n_steps;

    std::cout << "Computed dt = " << time_step << " with " << n_steps << " time steps.\n";
  }

  template <int dim>
  void
  WaveEquationProblem<dim>::assemble_forcing_term(const double current_time, dealii::Vector<double> &forcing) const
  {
    forcing = 0.0;

    RightHandSide         rhs_function(current_time);
    dealii::QGauss<dim>   quadrature_formula(fe.degree + 1);
    dealii::FEValues<dim> fe_values(
      fe, quadrature_formula, dealii::update_values | dealii::update_quadrature_points | dealii::update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    dealii::Vector<double>                       cell_rhs(dofs_per_cell);
    std::vector<dealii::types::global_dof_index> local_dof_indices(dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        fe_values.reinit(cell);
        cell_rhs = 0.0;
        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            const double f_value = rhs_function.value(fe_values.quadrature_point(q));
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                cell_rhs(i) += f_value * fe_values.shape_value(i, q) * fe_values.JxW(q);
              }
          }

        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(cell_rhs, local_dof_indices, forcing);
      }
  }

  template <int dim>
  std::map<dealii::types::global_dof_index, double>
  WaveEquationProblem<dim>::get_boundary_values(const double current_time) const
  {
    ExactSolution                                     exact_solution(current_time);
    std::map<dealii::types::global_dof_index, double> boundary_values;

    for (const dealii::types::boundary_id boundary_id : triangulation.get_boundary_ids())
      {
        dealii::VectorTools::interpolate_boundary_values(dof_handler, boundary_id, exact_solution, boundary_values);
      }

    return boundary_values;
  }

  template <int dim>
  void
  WaveEquationProblem<dim>::solve_mass_system(const std::map<dealii::types::global_dof_index, double> &boundary_values,
                                              dealii::Vector<double>                                  &destination,
                                              const dealii::Vector<double>                            &rhs) const
  {
    dealii::SparseMatrix<double> system_matrix;
    system_matrix.reinit(sparsity_pattern);
    system_matrix.copy_from(mass_matrix);

    dealii::Vector<double> rhs_with_bc(rhs);
    destination = 0.0;
    dealii::MatrixTools::apply_boundary_values(boundary_values, system_matrix, destination, rhs_with_bc);
    dealii::SparseDirectUMFPACK direct_solver;
    direct_solver.initialize(system_matrix);
    direct_solver.vmult(destination, rhs_with_bc);

    constraints.distribute(destination);
  }

  template <int dim>
  void
  WaveEquationProblem<dim>::initialize_solution()
  {
    ExactSolution   initial_displacement(0.0);
    InitialVelocity initial_velocity;

    dealii::Vector<double> velocity(dof_handler.n_dofs());
    dealii::VectorTools::interpolate(dof_handler, initial_displacement, solution);
    dealii::VectorTools::interpolate(dof_handler, initial_velocity, velocity);
    constraints.distribute(solution);
    constraints.distribute(velocity);

    output_results(0);

    dealii::Vector<double> mass_u0(dof_handler.n_dofs());
    dealii::Vector<double> mass_v0(dof_handler.n_dofs());
    dealii::Vector<double> stiffness_u0(dof_handler.n_dofs());

    mass_matrix.vmult(mass_u0, solution);
    mass_matrix.vmult(mass_v0, velocity);
    stiffness_matrix.vmult(stiffness_u0, solution);

    assemble_forcing_term(0.0, forcing_vector);

    dealii::Vector<double> first_step_rhs(dof_handler.n_dofs());
    first_step_rhs = mass_u0;
    first_step_rhs.add(time_step, mass_v0);
    first_step_rhs.add(0.5 * time_step * time_step, forcing_vector);
    first_step_rhs.add(-0.5 * time_step * time_step, stiffness_u0);

    solution_old = solution;
    solve_mass_system(get_boundary_values(time_step), solution, first_step_rhs);
    time            = time_step;
    timestep_number = 1;

    if (timestep_number % output_every == 0 || time >= final_time)
      {
        output_results(timestep_number);
      }
  }

  template <int dim>
  void
  WaveEquationProblem<dim>::output_results(const unsigned int timestep) const
  {
    std::filesystem::create_directories(output_directory);

    dealii::DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "u");
    data_out.build_patches();

    std::ofstream output(output_directory + "/solution-" + dealii::Utilities::int_to_string(timestep, 4) + ".vtu");
    data_out.write_vtu(output);
  }

  template <int dim>
  void
  WaveEquationProblem<dim>::compute_error() const
  {
    dealii::Vector<double> difference_per_cell(triangulation.n_active_cells());
    ExactSolution          exact_solution(time);

    dealii::VectorTools::integrate_difference(dof_handler,
                                              solution,
                                              exact_solution,
                                              difference_per_cell,
                                              dealii::QGauss<dim>(fe.degree + 2),
                                              dealii::VectorTools::L2_norm);

    std::cout << "Final time: " << time << '\n';
    std::cout << "L2 error: " << difference_per_cell.l2_norm() << '\n';
  }

  template <int dim>
  void
  WaveEquationProblem<dim>::run()
  {
    make_grid();
    setup_system();
    assemble_matrices();
    compute_time_step();
    initialize_solution();

    dealii::Vector<double> mass_u_n(dof_handler.n_dofs());
    dealii::Vector<double> mass_u_nm1(dof_handler.n_dofs());
    dealii::Vector<double> stiffness_u_n(dof_handler.n_dofs());

    while (time + 1e-14 < final_time)
      {
        assemble_forcing_term(time, forcing_vector);

        mass_matrix.vmult(mass_u_n, solution);
        mass_matrix.vmult(mass_u_nm1, solution_old);
        stiffness_matrix.vmult(stiffness_u_n, solution);

        system_rhs = 0.0;
        system_rhs.add(2.0, mass_u_n);
        system_rhs.add(-1.0, mass_u_nm1);
        system_rhs.add(-time_step * time_step, stiffness_u_n);
        system_rhs.add(time_step * time_step, forcing_vector);

        const double next_time = std::min(time + time_step, final_time);
        solve_mass_system(get_boundary_values(next_time), solution_new, system_rhs);

        solution_old = solution;
        solution     = solution_new;
        time         = next_time;
        ++timestep_number;

        if (timestep_number % output_every == 0 || time + 1e-14 >= final_time)
          {
            output_results(timestep_number);
          }
      }

    compute_error();
  }

  template class WaveEquationProblem<2>;
} // namespace WaveEquation
