#include "convergence_runner.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>

namespace
{
  double
  safe_order(const double e_coarse,
             const double e_fine,
             const double x_coarse,
             const double x_fine)
  {
    if (e_coarse <= 0.0 || e_fine <= 0.0 || x_coarse <= 0.0 || x_fine <= 0.0 ||
        std::abs(x_coarse - x_fine) < 1e-14)
      return 0.0;

    return std::log(e_coarse / e_fine) / std::log(x_coarse / x_fine);
  }
} // namespace

ConvergenceRunner::ConvergenceRunner(const WaveProblemConfig &config,
                                     MPI_Comm                 mpi_communicator)
  : base_config(config)
  , mpi_communicator(mpi_communicator)
{}

void
ConvergenceRunner::run()
{
  if (base_config.mode == WaveProblemConfig::Mode::convergence_space ||
      base_config.mode == WaveProblemConfig::Mode::convergence_both)
    {
      const auto entries = run_space_study();
      write_csv(base_config.convergence_csv_space, entries);
    }

  if (base_config.mode == WaveProblemConfig::Mode::convergence_time ||
      base_config.mode == WaveProblemConfig::Mode::convergence_both)
    {
      const auto entries = run_time_study();
      write_csv(base_config.convergence_csv_time, entries);
    }
}

std::vector<ConvergenceRunner::Entry>
ConvergenceRunner::run_space_study() const
{
  std::vector<Entry> entries;
  entries.reserve(base_config.convergence_mesh_files.size());

  for (const auto &mesh_file : base_config.convergence_mesh_files)
    {
      WaveProblemConfig cfg = base_config;
      cfg.mode = WaveProblemConfig::Mode::solve;
      cfg.mesh_file = mesh_file;
      cfg.output_interval = cfg.n_steps + 1;

      WaveSolver solver(cfg, mpi_communicator);
      solver.set_output_enabled(false);
      solver.run();

      const auto result = solver.compute_convergence_result();

      Entry entry;
      entry.study = "space";
      entry.mesh_file = mesh_file;
      entry.dt = cfg.dt;
      entry.n_steps = cfg.n_steps;
      entry.t_final = result.t_final;
      entry.h = result.h;
      entry.ndofs = result.ndofs;
      entry.l2_error = result.l2_error;
      entry.h1_error = result.h1_error;
      entries.push_back(entry);
    }

  std::sort(entries.begin(), entries.end(), [](const Entry &a, const Entry &b) {
    return a.h > b.h;
  });

  for (std::size_t i = 1; i < entries.size(); ++i)
    {
      entries[i].observed_order_l2 = safe_order(entries[i - 1].l2_error,
                                                entries[i].l2_error,
                                                entries[i - 1].h,
                                                entries[i].h);
      entries[i].observed_order_h1 = safe_order(entries[i - 1].h1_error,
                                                entries[i].h1_error,
                                                entries[i - 1].h,
                                                entries[i].h);
    }

  return entries;
}

std::vector<ConvergenceRunner::Entry>
ConvergenceRunner::run_time_study() const
{
  std::vector<Entry> entries;
  entries.reserve(base_config.convergence_dt_values.size());

  const double t_final = base_config.dt * static_cast<double>(base_config.n_steps);

  for (const auto dt_value : base_config.convergence_dt_values)
    {
      WaveProblemConfig cfg = base_config;
      cfg.mode = WaveProblemConfig::Mode::solve;
      cfg.dt = dt_value;
      cfg.n_steps = std::max(1, static_cast<int>(std::round(t_final / dt_value)));
      cfg.output_interval = cfg.n_steps + 1;

      WaveSolver solver(cfg, mpi_communicator);
      solver.set_output_enabled(false);
      solver.run();

      const auto result = solver.compute_convergence_result();

      Entry entry;
      entry.study = "time";
      entry.mesh_file = cfg.mesh_file;
      entry.dt = cfg.dt;
      entry.n_steps = cfg.n_steps;
      entry.t_final = result.t_final;
      entry.h = result.h;
      entry.ndofs = result.ndofs;
      entry.l2_error = result.l2_error;
      entry.h1_error = result.h1_error;
      entries.push_back(entry);
    }

  std::sort(entries.begin(), entries.end(), [](const Entry &a, const Entry &b) {
    return a.dt > b.dt;
  });

  for (std::size_t i = 1; i < entries.size(); ++i)
    {
      entries[i].observed_order_l2 = safe_order(entries[i - 1].l2_error,
                                                entries[i].l2_error,
                                                entries[i - 1].dt,
                                                entries[i].dt);
      entries[i].observed_order_h1 = safe_order(entries[i - 1].h1_error,
                                                entries[i].h1_error,
                                                entries[i - 1].dt,
                                                entries[i].dt);
    }

  return entries;
}

void
ConvergenceRunner::write_csv(const std::string &path,
                             const std::vector<Entry> &entries) const
{
  std::filesystem::path out_path(path);
  if (out_path.has_parent_path())
    std::filesystem::create_directories(out_path.parent_path());

  std::ofstream out(path);
  out << "study,method,fe_degree,theta,beta,gamma,mesh_file,dt,n_steps,t_final,h,ndofs,l2_error,h1_error,observed_order_l2,observed_order_h1\n";

  const std::string method_name =
    (base_config.method == WaveProblemConfig::Method::theta) ? "theta" : "newmark";

  out << std::setprecision(14);
  for (const auto &e : entries)
    {
      out << e.study << ',' << method_name << ',' << base_config.fe_degree << ','
          << base_config.theta << ',' << base_config.beta << ',' << base_config.gamma
          << ',' << e.mesh_file << ',' << e.dt << ',' << e.n_steps << ',' << e.t_final
          << ',' << e.h << ',' << e.ndofs << ',' << e.l2_error << ',' << e.h1_error
          << ',' << e.observed_order_l2 << ',' << e.observed_order_h1 << '\n';
    }
}
