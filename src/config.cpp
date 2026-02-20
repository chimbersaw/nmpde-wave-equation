#include "config.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <filesystem>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace
{
  std::vector<std::string>
  split_csv_strings(const std::string &value)
  {
    std::vector<std::string> tokens;
    std::stringstream        ss(value);
    std::string              token;

    while (std::getline(ss, token, ','))
      {
        token = trim_copy(token);
        if (!token.empty())
          tokens.push_back(token);
      }

    return tokens;
  }

  std::vector<double>
  split_csv_doubles(const std::string &value)
  {
    std::vector<double> out;
    for (const auto &entry : split_csv_strings(value))
      out.push_back(std::stod(entry));
    return out;
  }

  WaveProblemConfig::Mode
  parse_mode(const std::string &value)
  {
    if (value == "solve")
      return WaveProblemConfig::Mode::solve;
    if (value == "convergence_space")
      return WaveProblemConfig::Mode::convergence_space;
    if (value == "convergence_time")
      return WaveProblemConfig::Mode::convergence_time;
    if (value == "convergence_both")
      return WaveProblemConfig::Mode::convergence_both;

    throw std::runtime_error("Invalid mode: " + value);
  }

  WaveProblemConfig::Method
  parse_method(const std::string &value)
  {
    if (value == "theta")
      return WaveProblemConfig::Method::theta;
    if (value == "newmark")
      return WaveProblemConfig::Method::newmark;

    throw std::runtime_error("Invalid method: " + value);
  }

  void
  validate_config(const WaveProblemConfig &cfg)
  {
    if (cfg.mesh_file.empty())
      throw std::runtime_error("mesh_file must be specified.");
    if (cfg.fe_degree < 1)
      throw std::runtime_error("fe_degree must be >= 1.");
    if (cfg.wave_speed <= 0.0)
      throw std::runtime_error("wave_speed must be > 0.");
    if (cfg.dt <= 0.0)
      throw std::runtime_error("dt must be > 0.");
    if (cfg.n_steps <= 0)
      throw std::runtime_error("n_steps must be > 0.");
    if (cfg.output_interval <= 0)
      throw std::runtime_error("output_interval must be > 0.");

    if (cfg.method == WaveProblemConfig::Method::theta)
      {
        if (cfg.theta < 0.0 || cfg.theta > 1.0)
          throw std::runtime_error("theta must be in [0, 1].");
      }

    if (cfg.method == WaveProblemConfig::Method::newmark)
      {
        if (cfg.beta < 0.0)
          throw std::runtime_error("beta must be >= 0.");
      }

    if (cfg.mode == WaveProblemConfig::Mode::convergence_space ||
        cfg.mode == WaveProblemConfig::Mode::convergence_both)
      {
        if (cfg.convergence_mesh_files.empty())
          throw std::runtime_error(
            "convergence_mesh_files must be provided for spatial convergence.");
      }

    if (cfg.mode == WaveProblemConfig::Mode::convergence_time ||
        cfg.mode == WaveProblemConfig::Mode::convergence_both)
      {
        if (cfg.convergence_dt_values.empty())
          throw std::runtime_error(
            "convergence_dt_values must be provided for temporal convergence.");
      }
  }
} // namespace

std::string
trim_copy(const std::string &input)
{
  auto begin = std::find_if_not(
    input.begin(), input.end(), [](unsigned char c) { return std::isspace(c); });
  auto end = std::find_if_not(input.rbegin(),
                              input.rend(),
                              [](unsigned char c) { return std::isspace(c); })
               .base();

  if (begin >= end)
    return "";

  return std::string(begin, end);
}

WaveProblemConfig
parse_config_file(const std::string &path)
{
  std::ifstream in(path);
  if (!in)
    throw std::runtime_error("Cannot open config file: " + path);

  WaveProblemConfig cfg;
  const std::filesystem::path config_path(path);
  const std::filesystem::path config_dir = config_path.parent_path();

  const std::set<std::string> known_keys = {
    "mode",
    "method",
    "mesh_file",
    "fe_degree",
    "wave_speed",
    "dt",
    "n_steps",
    "output_interval",
    "output_dir",
    "scenario_u0",
    "scenario_u1",
    "scenario_f",
    "scenario_bc",
    "theta",
    "beta",
    "gamma",
    "convergence_mesh_files",
    "convergence_dt_values",
    "convergence_reference_case",
    "convergence_csv_space",
    "convergence_csv_time"};

  std::string line;
  int         line_no = 0;
  while (std::getline(in, line))
    {
      ++line_no;
      line = trim_copy(line);
      if (line.empty() || line[0] == '#')
        continue;

      const auto pos = line.find('=');
      if (pos == std::string::npos)
        throw std::runtime_error("Malformed config line " +
                                 std::to_string(line_no) + ": " + line);

      const std::string key = trim_copy(line.substr(0, pos));
      const std::string value = trim_copy(line.substr(pos + 1));

      if (known_keys.find(key) == known_keys.end())
        throw std::runtime_error("Unknown config key '" + key +
                                 "' on line " + std::to_string(line_no));

      if (key == "mode")
        cfg.mode = parse_mode(value);
      else if (key == "method")
        cfg.method = parse_method(value);
      else if (key == "mesh_file")
        cfg.mesh_file = value;
      else if (key == "fe_degree")
        cfg.fe_degree = std::stoi(value);
      else if (key == "wave_speed")
        cfg.wave_speed = std::stod(value);
      else if (key == "dt")
        cfg.dt = std::stod(value);
      else if (key == "n_steps")
        cfg.n_steps = std::stoi(value);
      else if (key == "output_interval")
        cfg.output_interval = std::stoi(value);
      else if (key == "output_dir")
        cfg.output_dir = value;
      else if (key == "scenario_u0")
        cfg.scenario_u0 = value;
      else if (key == "scenario_u1")
        cfg.scenario_u1 = value;
      else if (key == "scenario_f")
        cfg.scenario_f = value;
      else if (key == "scenario_bc")
        cfg.scenario_bc = value;
      else if (key == "theta")
        cfg.theta = std::stod(value);
      else if (key == "beta")
        cfg.beta = std::stod(value);
      else if (key == "gamma")
        cfg.gamma = std::stod(value);
      else if (key == "convergence_mesh_files")
        cfg.convergence_mesh_files = split_csv_strings(value);
      else if (key == "convergence_dt_values")
        cfg.convergence_dt_values = split_csv_doubles(value);
      else if (key == "convergence_reference_case")
        cfg.convergence_reference_case = value;
      else if (key == "convergence_csv_space")
        cfg.convergence_csv_space = value;
      else if (key == "convergence_csv_time")
        cfg.convergence_csv_time = value;
    }

  auto resolve_relative = [&config_dir](const std::string &p) -> std::string {
    if (p.empty())
      return p;
    const std::filesystem::path candidate(p);
    if (candidate.is_absolute() || config_dir.empty())
      return p;
    return (config_dir / candidate).lexically_normal().string();
  };

  cfg.mesh_file = resolve_relative(cfg.mesh_file);
  cfg.output_dir = resolve_relative(cfg.output_dir);
  cfg.convergence_csv_space = resolve_relative(cfg.convergence_csv_space);
  cfg.convergence_csv_time = resolve_relative(cfg.convergence_csv_time);
  for (auto &mesh : cfg.convergence_mesh_files)
    mesh = resolve_relative(mesh);

  validate_config(cfg);
  return cfg;
}
