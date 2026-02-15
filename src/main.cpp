#include "WaveEquation.hpp"

int main(int argc, char *argv[])
{
  // Inizializzazione MPI 
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);
  
  // Dimensione spaziale 2D per l'equazione delle onde
  const unsigned int dim = 2;

  // Parametri della simulazione:
  // 1. Percorso del file mesh (.msh)
  // 2. Grado polinomiale degli elementi finiti (1 = lineari)
  // 3. Tempo finale della simulazione (T)
  // 4. Passo temporale (delta_t)
  WaveEquation<dim> wave_problem("../mesh/mesh-square.msh", 1, 2.0, 0.01);

  // Lancio della simulazione
  wave_problem.run();

  return 0;
}


















