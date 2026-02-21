#include "../newmark_integrator.hpp"

#include "../wave_solver.hpp"

#include <cmath>

NewmarkTimeIntegrator::NewmarkTimeIntegrator(WaveSolver &solver,
                                             const double beta,
                                             const double gamma)
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
  const auto &M = solver.get_mass_matrix();
  const auto &K = solver.get_stiffness_matrix();

  const double dt = config.dt;
  const double c2 = config.wave_speed * config.wave_speed;

  WaveSolver::MatrixType displacement_matrix = solver.create_matrix_like_system();
  if (beta > 0.0)
    {
      displacement_matrix.copy_from(M);
      displacement_matrix *= 1.0 / (beta * dt * dt);
      displacement_matrix.add(c2, K);
    }

  for (int step = 0; step < config.n_steps; ++step)
    {
      const double t = step * dt;
      const double t_next = (step + 1) * dt;

      auto force_next = solver.create_owned_vector();
      solver.assemble_force(t_next, force_next);

      auto u_new = solver.create_owned_vector();
      auto v_new = solver.create_owned_vector();
      auto a_new = solver.create_owned_vector();

      if (beta > 0.0)
        {
          auto u_pred = solver.create_owned_vector();
          u_pred = u;
          u_pred.add(dt, v);
          u_pred.add(dt * dt * (0.5 - beta), a);

          auto rhs = solver.create_owned_vector();
          auto mu_pred = solver.create_owned_vector();

          M.vmult(mu_pred, u_pred);
          rhs = force_next;
          rhs.add(1.0 / (beta * dt * dt), mu_pred);

          u_new = u;
          solver.solve_displacement_system(displacement_matrix, u_new, rhs, t_next);

          a_new = u_new;
          a_new.add(-1.0, u_pred);
          a_new *= 1.0 / (beta * dt * dt);

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

          auto rhs = solver.create_owned_vector();
          auto ku = solver.create_owned_vector();
          K.vmult(ku, u_new);
          rhs = force_next;
          rhs.add(-c2, ku);

          solver.solve_spd_system(M, a_new, rhs);

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
