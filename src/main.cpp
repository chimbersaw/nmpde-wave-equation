#include <deal.II/base/mpi.h>

#include <exception>
#include <iostream>

#include "Options.hpp"
#include "WaveEquation.hpp"

int
main(int argc, char *argv[])
{
  try
    {
      const CliOptions cli = parse_cli_options(argc, argv);
      if (cli.show_help)
        {
          print_help();
          return 0;
        }

      const WaveEquation::Options options = load_options_from_config(cli.config_path);

      int    dealii_argc     = 1;
      char  *dealii_argv[]   = {argv[0], nullptr};
      char **dealii_argv_ptr = dealii_argv;

      dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(dealii_argc, dealii_argv_ptr, 1);

      WaveEquation::WaveEquationProblem problem(options);
      problem.run();
    }
  catch (const std::exception &exception)
    {
      std::cerr << exception.what() << '\n';
      return 1;
    }

  return 0;
}
