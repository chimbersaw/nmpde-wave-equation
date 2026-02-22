#include "wave_solver.hpp"

#include <deal.II/base/mpi.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_description.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <cmath>
#include <filesystem>
#include <limits>
#include <set>

#include "function_factory.hpp"
#include "newmark_integrator.hpp"
#include "theta_integrator.hpp"

namespace
{
  bool
  is_absorbing_boundary_condition(const std::string &bc_name)
  {
    return (bc_name == "absorbing");
  }
} // namespace

WaveSolver::WaveSolver(const WaveProblemConfig &config, MPI_Comm mpi_communicator)
  : mpi_communicator(mpi_communicator)
  , config(config)
  , pcout(std::cout, dealii::Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
  , triangulation(mpi_communicator)
  , fe(static_cast<unsigned int>(config.fe_degree))
  , dof_handler(triangulation)
  , quadrature(static_cast<unsigned int>(config.fe_degree + 2))
{
  u0_function             = make_named_function(config.scenario_u0, config.wave_speed);
  u1_function             = make_named_function(config.scenario_u1, config.wave_speed);
  forcing_function        = make_named_function(config.scenario_f, config.wave_speed);
  sigma_function          = make_named_function(config.scenario_sigma, config.wave_speed);
  boundary_function       = make_named_function(config.scenario_bc, config.wave_speed);
  exact_solution_function = make_named_function(config.convergence_reference_case, config.wave_speed);

  setup_triangulation_and_dofs();
  assemble_system_matrices();
  initialize_state();
}

void
WaveSolver::set_output_enabled(const bool enabled)
{
  output_enabled = enabled;
}

const WaveProblemConfig &
WaveSolver::get_config() const
{
  return config;
}

const WaveSolver::MatrixType &
WaveSolver::get_mass_matrix() const
{
  return mass_matrix;
}

const WaveSolver::MatrixType &
WaveSolver::get_stiffness_matrix() const
{
  return stiffness_matrix;
}

const WaveSolver::MatrixType &
WaveSolver::get_damping_matrix() const
{
  return damping_matrix;
}

WaveSolver::VectorType &
WaveSolver::get_solution()
{
  return solution;
}

WaveSolver::VectorType &
WaveSolver::get_velocity()
{
  return velocity;
}

WaveSolver::VectorType &
WaveSolver::get_acceleration()
{
  return acceleration;
}

WaveSolver::VectorType
WaveSolver::create_owned_vector() const
{
  VectorType vec;
  vec.reinit(locally_owned_dofs, mpi_communicator);
  return vec;
}

WaveSolver::MatrixType
WaveSolver::create_matrix_like_system() const
{
  dealii::DynamicSparsityPattern dsp(locally_relevant_dofs);
  dealii::DoFTools::make_sparsity_pattern(dof_handler, dsp);
  dealii::SparsityTools::distribute_sparsity_pattern(dsp, locally_owned_dofs, mpi_communicator, locally_relevant_dofs);

  MatrixType matrix;
  matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp, mpi_communicator);
  return matrix;
}

void
WaveSolver::setup_triangulation_and_dofs()
{
  dealii::Triangulation<2> serial_triangulation;

  dealii::GridIn<2> grid_in;
  grid_in.attach_triangulation(serial_triangulation);
  grid_in.read_msh(config.mesh_file);

  dealii::GridTools::partition_triangulation(dealii::Utilities::MPI::n_mpi_processes(mpi_communicator),
                                             serial_triangulation);

  const auto description =
    dealii::TriangulationDescription::Utilities::create_description_from_triangulation(serial_triangulation,
                                                                                       mpi_communicator);

  triangulation.create_triangulation(description);

  dof_handler.distribute_dofs(fe);

  locally_owned_dofs    = dof_handler.locally_owned_dofs();
  locally_relevant_dofs = dealii::DoFTools::extract_locally_relevant_dofs(dof_handler);

  dealii::DynamicSparsityPattern dsp(locally_relevant_dofs);
  dealii::DoFTools::make_sparsity_pattern(dof_handler, dsp);
  dealii::SparsityTools::distribute_sparsity_pattern(dsp, locally_owned_dofs, mpi_communicator, locally_relevant_dofs);

  mass_matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp, mpi_communicator);
  stiffness_matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp, mpi_communicator);
  damping_matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp, mpi_communicator);

  solution.reinit(locally_owned_dofs, mpi_communicator);
  velocity.reinit(locally_owned_dofs, mpi_communicator);
  acceleration.reinit(locally_owned_dofs, mpi_communicator);

  pcout << "Triangulation and DoFs ready: n_dofs=" << dof_handler.n_dofs() << "\n";
}

void
WaveSolver::assemble_system_matrices()
{
  dealii::FEValues<2>     fe_values(fe,
                                quadrature,
                                dealii::update_values | dealii::update_gradients | dealii::update_JxW_values);
  const dealii::QGauss<1> face_quadrature(static_cast<unsigned int>(config.fe_degree + 2));
  dealii::FEFaceValues<2> fe_face_values(fe, face_quadrature, dealii::update_values | dealii::update_JxW_values);
  const bool              absorbing_bc = is_absorbing_boundary_condition(config.scenario_bc);

  const unsigned int dofs_per_cell   = fe.n_dofs_per_cell();
  const unsigned int n_q_points      = quadrature.size();
  const unsigned int n_face_q_points = face_quadrature.size();

  dealii::FullMatrix<double>                   cell_mass(dofs_per_cell, dofs_per_cell);
  dealii::FullMatrix<double>                   cell_stiffness(dofs_per_cell, dofs_per_cell);
  dealii::FullMatrix<double>                   cell_damping(dofs_per_cell, dofs_per_cell);
  std::vector<dealii::types::global_dof_index> local_dof_indices(dofs_per_cell);
  std::vector<double>                          sigma_values(n_q_points, 0.0);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      cell_mass      = 0.0;
      cell_stiffness = 0.0;
      cell_damping   = 0.0;
      fe_values.reinit(cell);
      sigma_function->value_list(fe_values.get_quadrature_points(), sigma_values);

      for (unsigned int q = 0; q < n_q_points; ++q)
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            {
              cell_mass(i, j) += fe_values.shape_value(i, q) * fe_values.shape_value(j, q) * fe_values.JxW(q);

              cell_stiffness(i, j) += fe_values.shape_grad(i, q) * fe_values.shape_grad(j, q) * fe_values.JxW(q);

              cell_damping(i, j) +=
                sigma_values[q] * fe_values.shape_value(i, q) * fe_values.shape_value(j, q) * fe_values.JxW(q);
            }

      if (absorbing_bc)
        {
          for (const auto face_no : cell->face_indices())
            if (cell->face(face_no)->at_boundary())
              {
                fe_face_values.reinit(cell, face_no);
                for (unsigned int q = 0; q < n_face_q_points; ++q)
                  for (unsigned int i = 0; i < dofs_per_cell; ++i)
                    for (unsigned int j = 0; j < dofs_per_cell; ++j)
                      cell_damping(i, j) += config.wave_speed * fe_face_values.shape_value(i, q) *
                                            fe_face_values.shape_value(j, q) * fe_face_values.JxW(q);
              }
        }

      cell->get_dof_indices(local_dof_indices);
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
          {
            mass_matrix.add(local_dof_indices[i], local_dof_indices[j], cell_mass(i, j));
            stiffness_matrix.add(local_dof_indices[i], local_dof_indices[j], cell_stiffness(i, j));
            damping_matrix.add(local_dof_indices[i], local_dof_indices[j], cell_damping(i, j));
          }
    }

  mass_matrix.compress(dealii::VectorOperation::add);
  stiffness_matrix.compress(dealii::VectorOperation::add);
  damping_matrix.compress(dealii::VectorOperation::add);
}

void
WaveSolver::initialize_state()
{
  dealii::VectorTools::interpolate(dof_handler, *u0_function, solution);
  dealii::VectorTools::interpolate(dof_handler, *u1_function, velocity);

  enforce_displacement_bc(solution, 0.0);
  enforce_velocity_bc(velocity, 0.0, config.dt, config.dt);

  VectorType f0 = create_owned_vector();
  assemble_force(0.0, f0);

  VectorType ku = create_owned_vector();
  stiffness_matrix.vmult(ku, solution);

  VectorType rhs = create_owned_vector();
  rhs            = f0;
  VectorType cv  = create_owned_vector();
  damping_matrix.vmult(cv, velocity);
  rhs.add(-1.0, cv);
  rhs.add(-config.wave_speed * config.wave_speed, ku);

  solve_spd_system(mass_matrix, acceleration, rhs);
  enforce_acceleration_bc(acceleration, -config.dt, 0.0, config.dt, config.dt);
}

void
WaveSolver::assemble_force(const double time, VectorType &force) const
{
  force = 0.0;

  forcing_function->set_time(time);

  dealii::FEValues<2> fe_values(fe,
                                quadrature,
                                dealii::update_values | dealii::update_quadrature_points | dealii::update_JxW_values);

  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
  const unsigned int n_q_points    = quadrature.size();

  dealii::Vector<double>                       cell_rhs(dofs_per_cell);
  std::vector<double>                          forcing_values(n_q_points, 0.0);
  std::vector<dealii::types::global_dof_index> local_dof_indices(dofs_per_cell);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      fe_values.reinit(cell);
      cell_rhs = 0.0;

      forcing_function->value_list(fe_values.get_quadrature_points(), forcing_values);

      for (unsigned int q = 0; q < n_q_points; ++q)
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          cell_rhs(i) += forcing_values[q] * fe_values.shape_value(i, q) * fe_values.JxW(q);

      cell->get_dof_indices(local_dof_indices);
      std::vector<double> local_values(dofs_per_cell, 0.0);
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        local_values[i] = cell_rhs(i);
      force.add(local_dof_indices, local_values);
    }

  force.compress(dealii::VectorOperation::add);
}

void
WaveSolver::build_constraints(const double                                       time,
                              dealii::AffineConstraints<double>                 &constraints,
                              std::map<dealii::types::global_dof_index, double> &boundary_values) const
{
  constraints.clear();
  boundary_values.clear();

  dealii::DoFTools::make_hanging_node_constraints(dof_handler, constraints);

  if (is_absorbing_boundary_condition(config.scenario_bc))
    {
      constraints.close();
      return;
    }

  auto bc_function = make_named_function(config.scenario_bc, config.wave_speed);
  bc_function->set_time(time);

  std::set<dealii::types::boundary_id> boundary_ids;
  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      for (const auto face_no : cell->face_indices())
        if (cell->face(face_no)->at_boundary())
          boundary_ids.insert(cell->face(face_no)->boundary_id());

  for (const auto id : boundary_ids)
    dealii::VectorTools::interpolate_boundary_values(dof_handler, id, *bc_function, boundary_values);

  for (const auto &entry : boundary_values)
    {
      constraints.add_line(entry.first);
      constraints.set_inhomogeneity(entry.first, entry.second);
    }

  constraints.close();
}

void
WaveSolver::solve_spd_system(const MatrixType &A, VectorType &x, const VectorType &b) const
{
  dealii::SolverControl              solver_control(2000, 1e-12 * b.l2_norm() + 1e-14);
  dealii::TrilinosWrappers::SolverCG solver(solver_control);

  dealii::TrilinosWrappers::PreconditionAMG                 preconditioner;
  dealii::TrilinosWrappers::PreconditionAMG::AdditionalData amg_data;
  amg_data.elliptic              = true;
  amg_data.higher_order_elements = (config.fe_degree > 1);
  preconditioner.initialize(A, amg_data);

  x = 0.0;
  solver.solve(A, x, b, preconditioner);
}

void
WaveSolver::solve_displacement_system(const MatrixType &A,
                                      VectorType       &u,
                                      const VectorType &rhs,
                                      const double      time) const
{
  MatrixType system_matrix = create_matrix_like_system();
  system_matrix.copy_from(A);

  VectorType system_rhs = create_owned_vector();
  system_rhs            = rhs;

  dealii::AffineConstraints<double>                 constraints;
  std::map<dealii::types::global_dof_index, double> boundary_values;
  build_constraints(time, constraints, boundary_values);

  dealii::MatrixTools::apply_boundary_values(boundary_values, system_matrix, u, system_rhs, false);

  solve_spd_system(system_matrix, u, system_rhs);

  constraints.distribute(u);
}

void
WaveSolver::enforce_displacement_bc(VectorType &u, const double time) const
{
  dealii::AffineConstraints<double>                 constraints;
  std::map<dealii::types::global_dof_index, double> boundary_values;
  build_constraints(time, constraints, boundary_values);

  for (const auto &entry : boundary_values)
    if (u.in_local_range(entry.first))
      u[entry.first] = entry.second;

  u.compress(dealii::VectorOperation::insert);
  constraints.distribute(u);
}

void
WaveSolver::enforce_velocity_bc(VectorType  &v,
                                const double previous_time,
                                const double current_time,
                                const double dt) const
{
  if (is_absorbing_boundary_condition(config.scenario_bc))
    return;

  std::map<dealii::types::global_dof_index, double> boundary_old;
  std::map<dealii::types::global_dof_index, double> boundary_new;

  {
    auto bc_old = make_named_function(config.scenario_bc, config.wave_speed);
    bc_old->set_time(previous_time);

    std::set<dealii::types::boundary_id> boundary_ids;
    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        for (const auto face_no : cell->face_indices())
          if (cell->face(face_no)->at_boundary())
            boundary_ids.insert(cell->face(face_no)->boundary_id());

    for (const auto id : boundary_ids)
      dealii::VectorTools::interpolate_boundary_values(dof_handler, id, *bc_old, boundary_old);
  }

  {
    auto bc_new = make_named_function(config.scenario_bc, config.wave_speed);
    bc_new->set_time(current_time);

    std::set<dealii::types::boundary_id> boundary_ids;
    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        for (const auto face_no : cell->face_indices())
          if (cell->face(face_no)->at_boundary())
            boundary_ids.insert(cell->face(face_no)->boundary_id());

    for (const auto id : boundary_ids)
      dealii::VectorTools::interpolate_boundary_values(dof_handler, id, *bc_new, boundary_new);
  }

  for (const auto &entry : boundary_new)
    {
      const auto   it_old      = boundary_old.find(entry.first);
      const double old_value   = (it_old == boundary_old.end()) ? 0.0 : it_old->second;
      const double new_value   = entry.second;
      const double bc_velocity = (new_value - old_value) / dt;

      if (v.in_local_range(entry.first))
        v[entry.first] = bc_velocity;
    }

  v.compress(dealii::VectorOperation::insert);
}

void
WaveSolver::enforce_acceleration_bc(VectorType  &a,
                                    const double previous_time,
                                    const double current_time,
                                    const double next_time,
                                    const double dt) const
{
  if (is_absorbing_boundary_condition(config.scenario_bc))
    return;

  std::map<dealii::types::global_dof_index, double> boundary_prev;
  std::map<dealii::types::global_dof_index, double> boundary_curr;
  std::map<dealii::types::global_dof_index, double> boundary_next;

  std::set<dealii::types::boundary_id> boundary_ids;
  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      for (const auto face_no : cell->face_indices())
        if (cell->face(face_no)->at_boundary())
          boundary_ids.insert(cell->face(face_no)->boundary_id());

  {
    auto bc_prev = make_named_function(config.scenario_bc, config.wave_speed);
    bc_prev->set_time(previous_time);
    for (const auto id : boundary_ids)
      dealii::VectorTools::interpolate_boundary_values(dof_handler, id, *bc_prev, boundary_prev);
  }

  {
    auto bc_curr = make_named_function(config.scenario_bc, config.wave_speed);
    bc_curr->set_time(current_time);
    for (const auto id : boundary_ids)
      dealii::VectorTools::interpolate_boundary_values(dof_handler, id, *bc_curr, boundary_curr);
  }

  {
    auto bc_next = make_named_function(config.scenario_bc, config.wave_speed);
    bc_next->set_time(next_time);
    for (const auto id : boundary_ids)
      dealii::VectorTools::interpolate_boundary_values(dof_handler, id, *bc_next, boundary_next);
  }

  for (const auto &entry : boundary_curr)
    {
      const auto   it_prev = boundary_prev.find(entry.first);
      const auto   it_next = boundary_next.find(entry.first);
      const double g_prev  = (it_prev == boundary_prev.end()) ? 0.0 : it_prev->second;
      const double g_curr  = entry.second;
      const double g_next  = (it_next == boundary_next.end()) ? 0.0 : it_next->second;
      const double a_bc    = (g_next - 2.0 * g_curr + g_prev) / (dt * dt);

      if (a.in_local_range(entry.first))
        a[entry.first] = a_bc;
    }

  a.compress(dealii::VectorOperation::insert);
}

bool
WaveSolver::should_output(const unsigned int step) const
{
  if (!output_enabled)
    return false;

  return (step == 0 || step == static_cast<unsigned int>(config.n_steps) ||
          step % static_cast<unsigned int>(config.output_interval) == 0);
}

void
WaveSolver::output_results(const unsigned int step, const double) const
{
  if (!output_enabled)
    return;

  std::filesystem::create_directories(config.output_dir);
  const std::string out_dir =
    (config.output_dir.empty() || config.output_dir.back() == '/') ? config.output_dir : config.output_dir + "/";

  VectorType u_relevant;
  u_relevant.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  u_relevant = solution;

  VectorType v_relevant;
  v_relevant.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  v_relevant = velocity;

  dealii::DataOut<2> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(u_relevant, "u");
  data_out.add_data_vector(v_relevant, "v");
  data_out.build_patches();

  data_out.write_vtu_with_pvtu_record(out_dir, "solution", step, mpi_communicator, 4);
}

void
WaveSolver::run()
{
  if (should_output(0))
    output_results(0, 0.0);

  if (config.method == WaveProblemConfig::Method::theta)
    {
      ThetaTimeIntegrator integrator(*this, config.theta);
      integrator.run();
    }
  else
    {
      NewmarkTimeIntegrator integrator(*this, config.beta, config.gamma);
      integrator.run();
    }
}

double
WaveSolver::compute_l2_error(const double time) const
{
  auto exact = make_named_function(config.convergence_reference_case, config.wave_speed);
  exact->set_time(time);

  VectorType u_relevant;
  u_relevant.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  u_relevant = solution;
  u_relevant.compress(dealii::VectorOperation::insert);

  dealii::Vector<float> difference_per_cell(triangulation.n_active_cells());
  difference_per_cell = 0.0;
  dealii::VectorTools::integrate_difference(dof_handler,
                                            u_relevant,
                                            *exact,
                                            difference_per_cell,
                                            dealii::QGaussSimplex<2>(static_cast<unsigned int>(config.fe_degree + 3)),
                                            dealii::VectorTools::L2_norm);

  double local_sq = 0.0;
  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        const double value = difference_per_cell(cell->active_cell_index());
        local_sq += value * value;
      }

  const double global_sq       = dealii::Utilities::MPI::sum(local_sq, mpi_communicator);
  const double vector_tools_l2 = std::sqrt(global_sq);
  if (std::isfinite(vector_tools_l2))
    return vector_tools_l2;

  dealii::FEValues<2> fe_values(fe,
                                quadrature,
                                dealii::update_values | dealii::update_gradients | dealii::update_quadrature_points |
                                  dealii::update_JxW_values);
  const unsigned int  dofs_per_cell = fe.n_dofs_per_cell();
  const unsigned int  n_q_points    = quadrature.size();

  std::vector<dealii::types::global_dof_index> local_dof_indices(dofs_per_cell);
  std::vector<double>                          exact_values(n_q_points);

  double local_sq_manual = 0.0;
  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        fe_values.reinit(cell);
        cell->get_dof_indices(local_dof_indices);
        exact->value_list(fe_values.get_quadrature_points(), exact_values);

        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            double uh_q = 0.0;
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              uh_q += u_relevant(local_dof_indices[i]) * fe_values.shape_value(i, q);

            const double diff = uh_q - exact_values[q];
            local_sq_manual += diff * diff * fe_values.JxW(q);
          }
      }

  const double global_sq_manual = dealii::Utilities::MPI::sum(local_sq_manual, mpi_communicator);
  return std::sqrt(global_sq_manual);
}

double
WaveSolver::compute_h1_error(const double time) const
{
  auto exact = make_named_function(config.convergence_reference_case, config.wave_speed);
  exact->set_time(time);

  VectorType u_relevant;
  u_relevant.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  u_relevant = solution;
  u_relevant.compress(dealii::VectorOperation::insert);

  dealii::Vector<float> difference_per_cell(triangulation.n_active_cells());
  difference_per_cell = 0.0;
  dealii::VectorTools::integrate_difference(dof_handler,
                                            u_relevant,
                                            *exact,
                                            difference_per_cell,
                                            dealii::QGaussSimplex<2>(static_cast<unsigned int>(config.fe_degree + 3)),
                                            dealii::VectorTools::H1_norm);

  double local_sq = 0.0;
  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        const double value = difference_per_cell(cell->active_cell_index());
        local_sq += value * value;
      }

  const double global_sq       = dealii::Utilities::MPI::sum(local_sq, mpi_communicator);
  const double vector_tools_h1 = std::sqrt(global_sq);
  if (std::isfinite(vector_tools_h1))
    return vector_tools_h1;

  dealii::FEValues<2> fe_values(fe,
                                quadrature,
                                dealii::update_values | dealii::update_gradients | dealii::update_quadrature_points |
                                  dealii::update_JxW_values);
  const unsigned int  dofs_per_cell = fe.n_dofs_per_cell();
  const unsigned int  n_q_points    = quadrature.size();

  std::vector<dealii::types::global_dof_index> local_dof_indices(dofs_per_cell);
  std::vector<double>                          exact_values(n_q_points);
  std::vector<dealii::Tensor<1, 2>>            exact_gradients(n_q_points);

  double local_sq_manual = 0.0;
  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        fe_values.reinit(cell);
        cell->get_dof_indices(local_dof_indices);
        exact->value_list(fe_values.get_quadrature_points(), exact_values);
        exact->gradient_list(fe_values.get_quadrature_points(), exact_gradients);

        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            double               uh_q = 0.0;
            dealii::Tensor<1, 2> grad_uh_q;
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                const double u_i = u_relevant(local_dof_indices[i]);
                uh_q += u_i * fe_values.shape_value(i, q);
                grad_uh_q += u_i * fe_values.shape_grad(i, q);
              }

            const double value_diff = uh_q - exact_values[q];
            const auto   grad_diff  = grad_uh_q - exact_gradients[q];
            local_sq_manual += (value_diff * value_diff + grad_diff.norm_square()) * fe_values.JxW(q);
          }
      }

  const double global_sq_manual = dealii::Utilities::MPI::sum(local_sq_manual, mpi_communicator);
  return std::sqrt(global_sq_manual);
}

double
WaveSolver::estimate_mesh_h() const
{
  double local_h = 0.0;
  for (const auto &cell : triangulation.active_cell_iterators())
    if (cell->is_locally_owned())
      local_h = std::max(local_h, cell->diameter());

  return dealii::Utilities::MPI::max(local_h, mpi_communicator);
}

ConvergenceResult
WaveSolver::compute_convergence_result() const
{
  ConvergenceResult result;
  result.t_final  = config.dt * static_cast<double>(config.n_steps);
  result.h        = estimate_mesh_h();
  result.ndofs    = static_cast<unsigned int>(dof_handler.n_dofs());
  result.l2_error = compute_l2_error(result.t_final);
  result.h1_error = compute_h1_error(result.t_final);
  return result;
}
