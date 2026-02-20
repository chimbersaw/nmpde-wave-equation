#include "function_factory.hpp"

#include <deal.II/base/numbers.h>

#include <cmath>
#include <stdexcept>

namespace
{
  using dealii::Function;
  using dealii::Point;
  using dealii::Tensor;

  class ZeroScalarFunction : public Function<2>
  {
  public:
    double
    value(const Point<2> &, const unsigned int) const override
    {
      return 0.0;
    }

    Tensor<1, 2>
    gradient(const Point<2> &, const unsigned int) const override
    {
      return {};
    }
  };

  class StandingWaveFunction : public Function<2>
  {
  public:
    explicit StandingWaveFunction(const double wave_speed)
      : Function<2>()
      , omega(wave_speed * dealii::numbers::PI * std::sqrt(2.0))
    {}

    double
    value(const Point<2> &p, const unsigned int) const override
    {
      return std::sin(dealii::numbers::PI * p[0]) *
             std::sin(dealii::numbers::PI * p[1]) *
             std::cos(omega * this->get_time());
    }

    Tensor<1, 2>
    gradient(const Point<2> &p, const unsigned int) const override
    {
      Tensor<1, 2> g;
      const double factor = std::cos(omega * this->get_time());
      g[0] = dealii::numbers::PI * std::cos(dealii::numbers::PI * p[0]) *
             std::sin(dealii::numbers::PI * p[1]) * factor;
      g[1] = dealii::numbers::PI * std::sin(dealii::numbers::PI * p[0]) *
             std::cos(dealii::numbers::PI * p[1]) * factor;
      return g;
    }

  private:
    const double omega;
  };

  class StandingWaveVelocityFunction : public Function<2>
  {
  public:
    explicit StandingWaveVelocityFunction(const double wave_speed)
      : Function<2>()
      , omega(wave_speed * dealii::numbers::PI * std::sqrt(2.0))
    {}

    double
    value(const Point<2> &p, const unsigned int) const override
    {
      return -omega * std::sin(dealii::numbers::PI * p[0]) *
             std::sin(dealii::numbers::PI * p[1]) *
             std::sin(omega * this->get_time());
    }

  private:
    const double omega;
  };

  class StandingWaveForcingFunction : public Function<2>
  {
  public:
    explicit StandingWaveForcingFunction(const double wave_speed)
      : Function<2>()
      , wave_speed(wave_speed)
      , omega(wave_speed * dealii::numbers::PI * std::sqrt(2.0))
    {}

    double
    value(const Point<2> &p, const unsigned int) const override
    {
      const double spatial = std::sin(dealii::numbers::PI * p[0]) *
                             std::sin(dealii::numbers::PI * p[1]);
      const double temporal = std::cos(omega * this->get_time());
      const double utt = -omega * omega * spatial * temporal;
      const double laplace = -2.0 * dealii::numbers::PI * dealii::numbers::PI *
                             spatial * temporal;
      return utt - wave_speed * wave_speed * laplace;
    }

  private:
    const double wave_speed;
    const double omega;
  };

  class GaussianPulseFunction : public Function<2>
  {
  public:
    double
    value(const Point<2> &p, const unsigned int) const override
    {
      const double dx = p[0] - 0.5;
      const double dy = p[1] - 0.5;
      return std::exp(-80.0 * (dx * dx + dy * dy));
    }
  };
} // namespace

std::shared_ptr<dealii::Function<2>>
make_named_function(const std::string &name, const double wave_speed)
{
  if (name == "zero" || name == "zero_dirichlet")
    return std::make_shared<ZeroScalarFunction>();
  if (name == "standing_wave" || name == "standing_wave_exact" ||
      name == "standing_wave_boundary")
    return std::make_shared<StandingWaveFunction>(wave_speed);
  if (name == "standing_wave_velocity")
    return std::make_shared<StandingWaveVelocityFunction>(wave_speed);
  if (name == "standing_wave_forcing")
    return std::make_shared<StandingWaveForcingFunction>(wave_speed);
  if (name == "gaussian_pulse")
    return std::make_shared<GaussianPulseFunction>();

  throw std::runtime_error("Unknown scenario function name: " + name);
}

bool
is_exact_solution_name(const std::string &name)
{
  return (name == "standing_wave_exact" || name == "standing_wave");
}
