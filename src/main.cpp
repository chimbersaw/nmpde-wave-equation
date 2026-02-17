#include <deal.II/base/mpi.h>

#include <cctype>
#include <exception>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include "WaveEquation.hpp"

namespace
{
  struct Options
  {
    WaveEquation::SetupId setup_id             = WaveEquation::SetupId::StandingWave;
    unsigned int          refinements          = 5;
    double                final_time           = 5.0;
    double                cfl_number           = 0.25;
    unsigned int          minimum_time_steps   = 500;
    unsigned int          output_every_n_steps = 1;
    double                source_amplitude     = 20.0;
    double                source_frequency     = 0.25;
    std::string           mesh_path;
  };

  std::string
  trim(const std::string &text)
  {
    std::size_t begin = 0;
    while (begin < text.size() && std::isspace(static_cast<unsigned char>(text[begin])))
      {
        ++begin;
      }

    std::size_t end = text.size();
    while (end > begin && std::isspace(static_cast<unsigned char>(text[end - 1])))
      {
        --end;
      }

    return text.substr(begin, end - begin);
  }

  WaveEquation::SetupId
  parse_setup_id(const std::string &name)
  {
    if (name == "standing")
      {
        return WaveEquation::SetupId::StandingWave;
      }
    if (name == "radial")
      {
        return WaveEquation::SetupId::RadialPulse;
      }
    if (name == "radial_absorbing")
      {
        return WaveEquation::SetupId::RadialAbsorbing;
      }
    throw std::runtime_error("Unknown setup '" + name + "'. Use 'standing', 'radial', or 'radial_absorbing'.");
  }

  Options
  load_options_from_config(const std::string &config_path)
  {
    std::ifstream input(config_path);
    if (!input.is_open())
      {
        throw std::runtime_error("Could not open config file: " + config_path);
      }

    Options      options;
    std::string  line;
    unsigned int line_number = 0;
    while (std::getline(input, line))
      {
        ++line_number;
        const auto comment_pos = line.find('#');
        if (comment_pos != std::string::npos)
          {
            line = line.substr(0, comment_pos);
          }

        line = trim(line);
        if (line.empty())
          {
            continue;
          }

        const auto equal_pos = line.find('=');
        if (equal_pos == std::string::npos)
          {
            throw std::runtime_error("Invalid config entry at line " + std::to_string(line_number) + ": " + line);
          }

        const std::string key   = trim(line.substr(0, equal_pos));
        const std::string value = trim(line.substr(equal_pos + 1));

        if (key == "setup")
          {
            options.setup_id = parse_setup_id(value);
          }
        else if (key == "refine")
          {
            options.refinements = static_cast<unsigned int>(std::stoul(value));
          }
        else if (key == "final_time")
          {
            options.final_time = std::stod(value);
          }
        else if (key == "cfl")
          {
            options.cfl_number = std::stod(value);
          }
        else if (key == "minimum_time_steps")
          {
            options.minimum_time_steps = static_cast<unsigned int>(std::stoul(value));
          }
        else if (key == "output_every")
          {
            options.output_every_n_steps = static_cast<unsigned int>(std::stoul(value));
          }
        else if (key == "mesh_path")
          {
            options.mesh_path = value;
          }
        else if (key == "source_amplitude")
          {
            options.source_amplitude = std::stod(value);
          }
        else if (key == "source_frequency")
          {
            options.source_frequency = std::stod(value);
          }
        else
          {
            throw std::runtime_error("Unknown config key at line " + std::to_string(line_number) + ": " + key);
          }
      }

    if (options.final_time <= 0.0)
      {
        throw std::runtime_error("Config error: final_time must be > 0");
      }
    if (options.cfl_number <= 0.0)
      {
        throw std::runtime_error("Config error: cfl must be > 0");
      }
    if (options.minimum_time_steps == 0)
      {
        throw std::runtime_error("Config error: minimum_time_steps must be >= 1");
      }
    if (options.output_every_n_steps == 0)
      {
        throw std::runtime_error("Config error: output_every must be >= 1");
      }
    if (options.source_frequency <= 0.0)
      {
        throw std::runtime_error("Config error: source_frequency must be > 0");
      }
    if (options.source_amplitude < 0.0)
      {
        throw std::runtime_error("Config error: source_amplitude must be >= 0");
      }

    return options;
  }

  void
  print_help()
  {
    std::cout << "Usage: wave-equation --config PATH\n\n"
              << "Options:\n"
              << "  --config PATH             Load simulation parameters from file\n"
              << "  --help                    Show this help\n";
  }
} // namespace

int
main(int argc, char *argv[])
{
  try
    {
      std::string config_path = "configs/standing.cfg";
      for (int i = 1; i < argc; ++i)
        {
          const std::string arg           = argv[i];
          const auto        require_value = [&](const std::string &flag) {
            if (i + 1 >= argc)
              {
                throw std::runtime_error("Missing value for " + flag);
              }
            return std::string(argv[++i]);
          };

          if (arg == "--help")
            {
              print_help();
              return 0;
            }
          if (arg == "--setup")
            {
              throw std::runtime_error("Use config file key 'setup' instead of CLI --setup.");
            }
          if (arg == "--config")
            {
              config_path = require_value(arg);
              continue;
            }
          throw std::runtime_error("Unknown option: " + arg);
        }

      const Options options = load_options_from_config(config_path);

      // Avoid deal.II treating custom flags as ParameterHandler options.
      int                                      dealii_argc     = 1;
      char                                    *dealii_argv[]   = {argv[0], nullptr};
      char                                   **dealii_argv_ptr = dealii_argv;
      dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(dealii_argc, dealii_argv_ptr, 1);

      WaveEquation::WaveEquationProblem<2> problem(options.refinements,
                                                   options.final_time,
                                                   options.cfl_number,
                                                   options.mesh_path,
                                                   options.setup_id,
                                                   options.minimum_time_steps,
                                                   options.output_every_n_steps,
                                                   options.source_amplitude,
                                                   options.source_frequency);

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
