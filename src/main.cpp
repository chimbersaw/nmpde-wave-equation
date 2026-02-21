#include "config.hpp"
#include "convergence_runner.hpp"
#include "wave_solver.hpp"

#include <deal.II/base/mpi.h>

#include <iostream>
#include <string>

namespace
{
  void
  print_help(const char *argv0)
  {
    std::cout << "Usage: " << argv0 << " --config <path/to/config.cfg>\n";
  }
} // namespace

int
main(int argc, char *argv[])
{
  try
    {
      dealii::Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);
      const MPI_Comm                           mpi_communicator(MPI_COMM_WORLD);

      std::string config_path;

      for (int i = 1; i < argc; ++i)
        {
          const std::string arg(argv[i]);
          if (arg == "--help")
            {
              print_help(argv[0]);
              return 0;
            }
          if (arg == "--config" && i + 1 < argc)
            {
              config_path = argv[++i];
              continue;
            }

          std::cerr << "Unknown argument: " << arg << '\n';
          print_help(argv[0]);
          return 1;
        }

      if (config_path.empty())
        {
          std::cerr << "Missing --config argument.\n";
          print_help(argv[0]);
          return 1;
        }

      const auto config = parse_config_file(config_path);

      if (config.mode == WaveProblemConfig::Mode::solve)
        {
          WaveSolver solver(config, mpi_communicator);
          solver.run();
        }
      else
        {
          ConvergenceRunner runner(config, mpi_communicator);
          runner.run();
        }

      return 0;
    }
  catch (const std::exception &e)
    {
      std::cerr << "Error: " << e.what() << '\n';
      return 1;
    }
}
