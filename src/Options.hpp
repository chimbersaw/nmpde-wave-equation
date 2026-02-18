#ifndef OPTIONS_HPP
#define OPTIONS_HPP

#include <string>

#include "WaveEquation.hpp"

struct CliOptions
{
  std::string config_path = "configs/standing.cfg";
  bool        show_help   = false;
};

CliOptions
parse_cli_options(int argc, char *argv[]);

WaveEquation::Options
load_options_from_config(const std::string &config_path);

void
print_help();

#endif
