#include <deal.II/base/mpi.h>

#include <exception>
#include <iostream>

#include "WaveEquation.hpp"

int
main(int argc, char *argv[])
{
  try
    {
      dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      const std::string                    mesh_path = (argc > 1) ? argv[1] : "mesh/unit_square.msh";
      WaveEquation::WaveEquationProblem<2> problem(5, 5.0, 0.25, mesh_path);

      problem.run();
    }
  catch (const std::exception &exception)
    {
      std::cerr << "Exception: " << exception.what() << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << "Unknown exception." << std::endl;
      return 1;
    }

  return 0;
}
