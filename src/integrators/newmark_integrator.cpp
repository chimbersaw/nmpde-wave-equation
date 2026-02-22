#include "../newmark_integrator.hpp"

#include <cmath>

#include "../wave_solver.hpp"

NewmarkTimeIntegrator::NewmarkTimeIntegrator(WaveSolver &solver, const double beta, const double gamma)
  : solver(solver)
  , beta(beta)
  , gamma(gamma)
{}

void
NewmarkTimeIntegrator::run()
{
  auto &u = solver.get_solution();
  auto &v = solver.get_velocity();
  auto &a = solver.get_acceleration();

  const auto &config = solver.get_config();
  const auto &M      = solver.get_mass_matrix();
  const auto &K      = solver.get_stiffness_matrix();
  const auto &C      = solver.get_damping_matrix();

  const double dt = config.dt;
  const double c2 = config.wave_speed * config.wave_speed;

  WaveSolver::MatrixType displacement_matrix = solver.create_matrix_like_system();
  if (beta > 0.0)
    {
      displacement_matrix.copy_from(M);
      displacement_matrix *= 1.0 / (beta * dt * dt);
      displacement_matrix.add(gamma / (beta * dt), C);
      displacement_matrix.add(c2, K);
    }

  WaveSolver::MatrixType acceleration_matrix = solver.create_matrix_like_system();
  if (std::abs(beta) < 1e-14)
    {
      acceleration_matrix.copy_from(M);
      acceleration_matrix.add(gamma * dt, C);
    }

  for (int step = 0; step < config.n_steps; ++step)
    {
      const double t      = step * dt;
      const double t_next = (step + 1) * dt;

      auto force_next = solver.create_owned_vector();
      solver.assemble_force(t_next, force_next);

      auto u_new = solver.create_owned_vector();
      auto v_new = solver.create_owned_vector();
      auto a_new = solver.create_owned_vector();

      if (beta > 0.0)
        {
          auto u_pred = solver.create_owned_vector();
          u_pred      = u;
          u_pred.add(dt, v);
          u_pred.add(dt * dt * (0.5 - beta), a);
          auto v_pred = solver.create_owned_vector();
          v_pred      = v;
          v_pred.add(dt * (1.0 - gamma), a);

          auto rhs      = solver.create_owned_vector();
          auto mu_pred  = solver.create_owned_vector();
          auto cv_input = solver.create_owned_vector();
          auto cv_term  = solver.create_owned_vector();

          M.vmult(mu_pred, u_pred);
          cv_input = u_pred;
          cv_input *= gamma / (beta * dt);
          cv_input.add(-1.0, v_pred);
          C.vmult(cv_term, cv_input);

          rhs = force_next;
          rhs.add(1.0 / (beta * dt * dt), mu_pred);
          rhs.add(1.0, cv_term);

          u_new = u;
          solver.solve_displacement_system(displacement_matrix, u_new, rhs, t_next);

          a_new = u_new;
          a_new.add(-1.0, u_pred);
          a_new *= 1.0 / (beta * dt * dt);
          solver.enforce_acceleration_bc(a_new, t, t_next, t_next + dt, dt);

          v_new = v;
          v_new.add(dt * (1.0 - gamma), a);
          v_new.add(dt * gamma, a_new);
          solver.enforce_velocity_bc(v_new, t, t_next, dt);
        }
      else
        {
          u_new = u;
          u_new.add(dt, v);
          u_new.add(0.5 * dt * dt, a);
          solver.enforce_displacement_bc(u_new, t_next);

          auto rhs    = solver.create_owned_vector();
          auto ku     = solver.create_owned_vector();
          auto v_pred = solver.create_owned_vector();
          auto cv     = solver.create_owned_vector();
          K.vmult(ku, u_new);
          v_pred = v;
          v_pred.add(dt * (1.0 - gamma), a);
          C.vmult(cv, v_pred);
          rhs = force_next;
          rhs.add(-1.0, cv);
          rhs.add(-c2, ku);

          solver.solve_spd_system(acceleration_matrix, a_new, rhs);
          solver.enforce_acceleration_bc(a_new, t, t_next, t_next + dt, dt);

          v_new = v;
          v_new.add(dt * (1.0 - gamma), a);
          v_new.add(dt * gamma, a_new);
          solver.enforce_velocity_bc(v_new, t, t_next, dt);
        }

      u = u_new;
      v = v_new;
      a = a_new;

      if (solver.should_output(static_cast<unsigned int>(step + 1)))
        solver.output_results(static_cast<unsigned int>(step + 1), t_next);
    }
}
