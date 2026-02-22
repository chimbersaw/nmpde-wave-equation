#include "../theta_integrator.hpp"

#include "../wave_solver.hpp"

#include <cmath>

ThetaTimeIntegrator::ThetaTimeIntegrator(WaveSolver &solver, const double theta)
  : solver(solver)
  , theta(theta)
{}

void
ThetaTimeIntegrator::run()
{
  auto &u = solver.get_solution();
  auto &v = solver.get_velocity();
  auto &a = solver.get_acceleration();

  const auto &config = solver.get_config();
  const auto &M = solver.get_mass_matrix();
  const auto &K = solver.get_stiffness_matrix();

  const double dt = config.dt;
  const double c2 = config.wave_speed * config.wave_speed;

  auto force_n = solver.create_owned_vector();
  auto force_np1 = solver.create_owned_vector();
  solver.assemble_force(0.0, force_n);

  WaveSolver::MatrixType effective_matrix = solver.create_matrix_like_system();
  if (theta > 0.0)
    {
      effective_matrix.copy_from(M);
      effective_matrix *= 1.0 / (dt * dt);
      effective_matrix.add(c2 * theta * theta, K);
    }

  for (int step = 0; step < config.n_steps; ++step)
    {
      const double t = step * dt;
      const double t_next = (step + 1) * dt;

      solver.assemble_force(t_next, force_np1);

      auto u_new = solver.create_owned_vector();
      auto v_new = solver.create_owned_vector();

      if (std::abs(theta) < 1e-14)
        {
          auto ku = solver.create_owned_vector();
          K.vmult(ku, u);

          auto rhs = solver.create_owned_vector();
          M.vmult(rhs, v);
          rhs.add(dt, force_n);
          rhs.add(-dt * c2, ku);

          solver.solve_spd_system(M, v_new, rhs);

          u_new = u;
          u_new.add(dt, v);
          solver.enforce_displacement_bc(u_new, t_next);
          solver.enforce_velocity_bc(v_new, t, t_next, dt);
        }
      else
        {
          auto rhs = solver.create_owned_vector();
          auto temp = solver.create_owned_vector();

          temp = u;
          temp *= 1.0 / (dt * dt);
          temp.add(1.0 / dt, v);
          M.vmult(rhs, temp);

          auto ku = solver.create_owned_vector();
          K.vmult(ku, u);
          rhs.add(-c2 * theta * (1.0 - theta), ku);

          rhs.add(theta * (1.0 - theta), force_n);
          rhs.add(theta * theta, force_np1);

          u_new = u;
          solver.solve_displacement_system(effective_matrix, u_new, rhs, t_next);

          v_new = u_new;
          v_new.add(-1.0, u);
          v_new *= 1.0 / (theta * dt);
          v_new.add(-((1.0 - theta) / theta), v);
          solver.enforce_velocity_bc(v_new, t, t_next, dt);
        }

      auto rhs_acc = solver.create_owned_vector();
      auto ku_new = solver.create_owned_vector();
      K.vmult(ku_new, u_new);
      rhs_acc = force_np1;
      rhs_acc.add(-c2, ku_new);
      solver.solve_spd_system(M, a, rhs_acc);
      solver.enforce_acceleration_bc(a, t, t_next, t_next + dt, dt);

      u = u_new;
      v = v_new;

      if (solver.should_output(static_cast<unsigned int>(step + 1)))
        solver.output_results(static_cast<unsigned int>(step + 1), t_next);

      force_n = force_np1;
    }
}
