#include "WaveEquation.hpp"

template <int dim>
WaveEquation<dim>::WaveEquation(const std::string &mesh_file_name_, 
                                 const unsigned int &r_,
                                 const double &T_, 
                                 const double &delta_t_)
  : mesh_file_name(mesh_file_name_), r(r_), T(T_), delta_t(delta_t_),
    time(0.0), timestep_number(0),
    mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)),
    mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)),
    mesh(MPI_COMM_WORLD),
    pcout(std::cout, mpi_rank == 0)
{}

template <int dim>
void WaveEquation<dim>::setup() {
  pcout << "===============================================" << std::endl;
  pcout << "Initializing the mesh" << std::endl;

  GridIn<dim> grid_in;
  grid_in.attach_triangulation(mesh);
  std::ifstream mesh_file(mesh_file_name);
  grid_in.read_msh(mesh_file);

  // Raffinamento globale per avere una risoluzione decente
  mesh.refine_global(4); 

  pcout << "Number of elements = " << mesh.n_global_active_cells() << std::endl;

  fe = std::make_unique<FE_SimplexP<dim>>(r);
  quadrature = std::make_unique<QGaussSimplex<dim>>(r + 1);
  dof_handler.reinit(mesh);
  dof_handler.distribute_dofs(*fe);

  pcout << "Number of DoFs     = " << dof_handler.n_dofs() << std::endl;

  const IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();
  const IndexSet locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);

  TrilinosWrappers::SparsityPattern sparsity(locally_owned_dofs, MPI_COMM_WORLD);
  DoFTools::make_sparsity_pattern(dof_handler, sparsity);
  sparsity.compress();

  system_matrix.reinit(sparsity);
  system_rhs.reinit(locally_owned_dofs, MPI_COMM_WORLD);
  
  solution.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
  old_solution.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
  older_solution.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
  solution_owned.reinit(locally_owned_dofs, MPI_COMM_WORLD);
}

template <int dim>
void WaveEquation<dim>::assemble_system() {
  system_matrix = 0;
  system_rhs = 0;

  FEValues<dim> fe_values(*fe, *quadrature, 
                          update_values | update_gradients | update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q = quadrature->size();
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double> cell_rhs(dofs_per_cell);
  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  std::vector<double> old_values(n_q);
  std::vector<double> older_values(n_q);
  std::vector<Tensor<1, dim>> old_grads(n_q);

  for (const auto &cell : dof_handler.active_cell_iterators()) {
    if (!cell->is_locally_owned()) continue;

    fe_values.reinit(cell);
    cell_matrix = 0;
    cell_rhs = 0;

    fe_values.get_function_values(old_solution, old_values);
    fe_values.get_function_values(older_solution, older_values);
    fe_values.get_function_gradients(old_solution, old_grads);

    for (unsigned int q = 0; q < n_q; ++q) {
      for (unsigned int i = 0; i < dofs_per_cell; ++i) {
        for (unsigned int j = 0; j < dofs_per_cell; ++j) {
          cell_matrix(i, j) += fe_values.shape_value(i, q) * fe_values.shape_value(j, q) * fe_values.JxW(q);
        }
        // Schema centrato: M*u_next = M*(2*u_n - u_old) - dt^2 * K*u_n
        cell_rhs(i) += (2.0 * old_values[q] - older_values[q]) * fe_values.shape_value(i, q) * fe_values.JxW(q);
        cell_rhs(i) -= (delta_t * delta_t) * (old_grads[q] * fe_values.shape_grad(i, q)) * fe_values.JxW(q);
      }
    }
    cell->get_dof_indices(dof_indices);
    system_matrix.add(dof_indices, cell_matrix);
    system_rhs.add(dof_indices, cell_rhs);
  }
  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);

  std::map<types::global_dof_index, double> boundary_values;
  VectorTools::interpolate_boundary_values(dof_handler, 0, Functions::ZeroFunction<dim>(), boundary_values);
  MatrixTools::apply_boundary_values(boundary_values, system_matrix, solution_owned, system_rhs);
}

template <int dim>
void WaveEquation<dim>::solve_timestep() {
  SolverControl solver_control(5000, 1e-12 * system_rhs.l2_norm());
  SolverCG<TrilinosWrappers::MPI::Vector> solver(solver_control);
  TrilinosWrappers::PreconditionSSOR preconditioner;
  preconditioner.initialize(system_matrix);

  solver.solve(system_matrix, solution_owned, system_rhs, preconditioner);
  solution = solution_owned;
}

template <int dim>
void WaveEquation<dim>::output_results() const {
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "u");
  data_out.build_patches();

  std::string filename = "solution-" + std::to_string(timestep_number) + ".vtu";
  std::ofstream output(filename);
  data_out.write_vtu(output);
}

/*template <int dim>
void WaveEquation<dim>::run() {
  setup();

  // Applichiamo la gaussiana iniziale
  VectorTools::interpolate(dof_handler, InitialCondition<dim>(), old_solution);
  older_solution = old_solution; 

  output_results();

  while (time < T - 0.5 * delta_t) {
    time += delta_t;
    timestep_number++;
    pcout << "Step " << timestep_number << " at t=" << time << " : ";

    assemble_system();
    solve_timestep();

    older_solution = old_solution;
    old_solution = solution;

    if (timestep_number % 5 == 0) {
        output_results();
        pcout << "Output saved.";
    }
    pcout << std::endl;
  }
}*/

template <int dim>
void WaveEquation<dim>::run() {
  setup();

  // INIZIALIZZAZIONE CORRETTA
  // Creiamo l'oggetto della condizione iniziale
  InitialCondition<dim> initial_condition_function;

  // Proiettiamo la gaussiana sui due passi temporali iniziali (u_n e u_{n-1})
  // Iniziamo con velocità zero, quindi u_n = u_{n-1}
  VectorTools::interpolate(dof_handler, initial_condition_function, old_solution);
  older_solution = old_solution; 
  
  // Molto importante: aggiorniamo anche solution_owned per l'output iniziale
  solution_owned = old_solution;

  // Salviamo lo stato iniziale (Step 0)
  output_results();

  while (time < T - 0.5 * delta_t) {
    time += delta_t;
    timestep_number++;
    
    // Assemblaggio e risoluzione per u_{n+1}
    assemble_system();
    solve_timestep();

    // Aggiornamento dei passi temporali
    older_solution = old_solution;
    old_solution = solution;

    // Salvataggio periodico
    if (timestep_number % 5 == 0) {
        output_results();
        pcout << "Step " << timestep_number << " at t=" << time << std::endl;
    }
  }
}

template class WaveEquation<2>;