#pragma once

#include <deal.II/base/function.h>

#include <memory>
#include <string>

std::shared_ptr<dealii::Function<2>>
make_named_function(const std::string &name, double wave_speed);

bool
is_exact_solution_name(const std::string &name);
