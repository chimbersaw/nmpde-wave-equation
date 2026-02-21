#pragma once

#include "config.hpp"
#include "wave_solver.hpp"

#include <string>
#include <vector>

class ConvergenceRunner
{
public:
  ConvergenceRunner(const WaveProblemConfig &config, MPI_Comm mpi_communicator);

  void
  run();

private:
  struct Entry
  {
    std::string study;
    std::string mesh_file;
    double      dt = 0.0;
    int         n_steps = 0;
    double      t_final = 0.0;
    double      h = 0.0;
    unsigned int ndofs = 0;
    double      l2_error = 0.0;
    double      h1_error = 0.0;
    double      observed_order_l2 = 0.0;
    double      observed_order_h1 = 0.0;
  };

  std::vector<Entry>
  run_space_study() const;

  std::vector<Entry>
  run_time_study() const;

  void
  write_csv(const std::string &path, const std::vector<Entry> &entries) const;

  WaveProblemConfig base_config;
  MPI_Comm          mpi_communicator;
};
