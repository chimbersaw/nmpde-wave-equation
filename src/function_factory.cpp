#include "function_factory.hpp"

#include <deal.II/base/numbers.h>

#include <algorithm>
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
      return std::sin(dealii::numbers::PI * p[0]) * std::sin(dealii::numbers::PI * p[1]) *
             std::cos(omega * this->get_time());
    }

    Tensor<1, 2>
    gradient(const Point<2> &p, const unsigned int) const override
    {
      Tensor<1, 2> g;
      const double factor = std::cos(omega * this->get_time());
      g[0] = dealii::numbers::PI * std::cos(dealii::numbers::PI * p[0]) * std::sin(dealii::numbers::PI * p[1]) * factor;
      g[1] = dealii::numbers::PI * std::sin(dealii::numbers::PI * p[0]) * std::cos(dealii::numbers::PI * p[1]) * factor;
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
      return -omega * std::sin(dealii::numbers::PI * p[0]) * std::sin(dealii::numbers::PI * p[1]) *
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
      const double spatial  = std::sin(dealii::numbers::PI * p[0]) * std::sin(dealii::numbers::PI * p[1]);
      const double temporal = std::cos(omega * this->get_time());
      const double utt      = -omega * omega * spatial * temporal;
      const double laplace  = -2.0 * dealii::numbers::PI * dealii::numbers::PI * spatial * temporal;
      return utt - wave_speed * wave_speed * laplace;
    }

  private:
    const double wave_speed;
    const double omega;
  };

  class StandingWave5x5Function : public Function<2>
  {
  public:
    explicit StandingWave5x5Function(const double wave_speed)
      : Function<2>()
      , omega(wave_speed * dealii::numbers::PI * std::sqrt(2.0) / 5.0)
    {}

    double
    value(const Point<2> &p, const unsigned int) const override
    {
      return std::sin((dealii::numbers::PI / 5.0) * p[0]) * std::sin((dealii::numbers::PI / 5.0) * p[1]) *
             std::cos(omega * this->get_time());
    }

    Tensor<1, 2>
    gradient(const Point<2> &p, const unsigned int) const override
    {
      Tensor<1, 2> g;
      const double factor = std::cos(omega * this->get_time());
      g[0]                = (dealii::numbers::PI / 5.0) * std::cos((dealii::numbers::PI / 5.0) * p[0]) *
             std::sin((dealii::numbers::PI / 5.0) * p[1]) * factor;
      g[1] = (dealii::numbers::PI / 5.0) * std::sin((dealii::numbers::PI / 5.0) * p[0]) *
             std::cos((dealii::numbers::PI / 5.0) * p[1]) * factor;
      return g;
    }

  private:
    const double omega;
  };

  class StandingWave5x5VelocityFunction : public Function<2>
  {
  public:
    explicit StandingWave5x5VelocityFunction(const double wave_speed)
      : Function<2>()
      , omega(wave_speed * dealii::numbers::PI * std::sqrt(2.0) / 5.0)
    {}

    double
    value(const Point<2> &p, const unsigned int) const override
    {
      return -omega * std::sin((dealii::numbers::PI / 5.0) * p[0]) * std::sin((dealii::numbers::PI / 5.0) * p[1]) *
             std::sin(omega * this->get_time());
    }

  private:
    const double omega;
  };

  class StandingWave5x5ForcingFunction : public Function<2>
  {
  public:
    explicit StandingWave5x5ForcingFunction(const double wave_speed)
      : Function<2>()
      , wave_speed(wave_speed)
      , omega(wave_speed * dealii::numbers::PI * std::sqrt(2.0) / 5.0)
    {}

    double
    value(const Point<2> &p, const unsigned int) const override
    {
      const double spatial =
        std::sin((dealii::numbers::PI / 5.0) * p[0]) * std::sin((dealii::numbers::PI / 5.0) * p[1]);
      const double temporal = std::cos(omega * this->get_time());
      const double utt      = -omega * omega * spatial * temporal;
      const double laplace  = -2.0 * (dealii::numbers::PI * dealii::numbers::PI / 25.0) * spatial * temporal;
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
      const double dx = p[0] - 2.5;
      const double dy = p[1] - 2.5;
      return std::exp(-8.0 * (dx * dx + dy * dy));
    }
  };

  class PeriodicCenterSource5x5Function : public Function<2>
  {
  public:
    double
    value(const Point<2> &p, const unsigned int) const override
    {
      const double dx = p[0] - 2.5;
      const double dy = p[1] - 2.5;
      const double r2 = dx * dx + dy * dy;

      const double period      = 0.5;
      const double pulse_width = 0.2;
      const double t           = this->get_time();
      const double t_mod       = std::fmod(std::max(0.0, t), period);

      double temporal = 0.0;
      if (t_mod <= pulse_width)
        {
          const double s    = t_mod / pulse_width;
          const double sine = std::sin(dealii::numbers::PI * s);
          temporal          = sine * sine;
        }

      const double spatial = std::exp(-26.0 * r2);
      return 240.0 * temporal * spatial;
    }
  };
} // namespace

std::shared_ptr<dealii::Function<2>>
make_named_function(const std::string &name, const double wave_speed)
{
  if (name == "zero" || name == "zero_dirichlet")
    return std::make_shared<ZeroScalarFunction>();
  if (name == "absorbing")
    return std::make_shared<ZeroScalarFunction>();
  if (name == "standing_wave" || name == "standing_wave_exact" || name == "standing_wave_boundary")
    return std::make_shared<StandingWaveFunction>(wave_speed);
  if (name == "standing_wave_5x5" || name == "standing_wave_5x5_exact" || name == "standing_wave_5x5_boundary")
    return std::make_shared<StandingWave5x5Function>(wave_speed);
  if (name == "standing_wave_velocity")
    return std::make_shared<StandingWaveVelocityFunction>(wave_speed);
  if (name == "standing_wave_5x5_velocity")
    return std::make_shared<StandingWave5x5VelocityFunction>(wave_speed);
  if (name == "standing_wave_forcing")
    return std::make_shared<StandingWaveForcingFunction>(wave_speed);
  if (name == "standing_wave_5x5_forcing")
    return std::make_shared<StandingWave5x5ForcingFunction>(wave_speed);
  if (name == "gaussian_pulse")
    return std::make_shared<GaussianPulseFunction>();
  if (name == "periodic_center_source" || name == "periodic_center_source_5x5")
    return std::make_shared<PeriodicCenterSource5x5Function>();

  throw std::runtime_error("Unknown scenario function name: " + name);
}
