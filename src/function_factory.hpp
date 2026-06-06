#pragma once

#include <deal.II/base/function.h>

#include <memory>
#include <string>

#include "wave_solver.hpp"

std::shared_ptr<dealii::Function<WaveSolver::dim>>
make_named_function(const std::string &name, double wave_speed);
