#pragma once

#include <string>
#include <vector>

struct WaveProblemConfig
{
  enum class Mode
  {
    solve,
    convergence_space,
    convergence_time,
    convergence_both
  };

  enum class Method
  {
    theta,
    newmark
  };

  Mode   mode   = Mode::solve;
  Method method = Method::theta;

  std::string mesh_file;
  int         fe_degree       = 1;
  double      wave_speed      = 1.0;
  double      dt              = 1e-3;
  int         n_steps         = 100;
  int         output_interval = 10;
  std::string output_dir      = "solution";

  std::string scenario_u0    = "standing_wave";
  std::string scenario_u1    = "standing_wave_velocity";
  std::string scenario_f     = "zero";
  std::string scenario_sigma = "zero";
  std::string scenario_bc    = "zero_dirichlet";

  double theta = 0.5;
  double beta  = 0.25;
  double gamma = 0.5;

  std::vector<std::string> convergence_mesh_files;
  std::vector<double>      convergence_dt_values;
  std::string              convergence_reference_case = "standing_wave_exact";
  std::string              convergence_csv_space      = "results/convergence_space.csv";
  std::string              convergence_csv_time       = "results/convergence_time.csv";
};

WaveProblemConfig
parse_config_file(const std::string &path);

std::string
trim_copy(const std::string &input);
