#include "Options.hpp"

#include <cctype>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace
{
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
    throw std::runtime_error("Unknown setup: " + name);
  }
} // namespace

CliOptions
parse_cli_options(const int argc, char *argv[])
{
  CliOptions cli;
  for (int i = 1; i < argc; ++i)
    {
      const std::string arg = argv[i];
      if (arg == "--help")
        {
          cli.show_help = true;
          continue;
        }
      if (arg == "--config" && i + 1 < argc)
        {
          cli.config_path = argv[++i];
          continue;
        }
      throw std::runtime_error("Usage: wave-equation --config PATH");
    }

  return cli;
}

WaveEquation::Options
load_options_from_config(const std::string &config_path)
{
  std::ifstream input(config_path);
  if (!input)
    {
      throw std::runtime_error("Could not open config file: " + config_path);
    }

  WaveEquation::Options options;
  std::string           line;
  unsigned int          line_number = 0;
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
          throw std::runtime_error("Invalid config line " + std::to_string(line_number));
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
      else if (key == "source_amplitude")
        {
          options.source_amplitude = std::stod(value);
        }
      else if (key == "source_frequency")
        {
          options.source_frequency = std::stod(value);
        }
      else if (key == "mesh_path")
        {
          options.mesh_path = value;
        }
      else
        {
          throw std::runtime_error("Unknown config key: " + key);
        }
    }

  if (options.final_time <= 0.0 || options.cfl_number <= 0.0 || options.minimum_time_steps == 0 ||
      options.output_every_n_steps == 0 || options.source_frequency <= 0.0 || options.source_amplitude < 0.0)
    {
      throw std::runtime_error("Invalid config values");
    }

  return options;
}

void
print_help()
{
  std::cout << "Usage: wave-equation --config PATH\n";
}
