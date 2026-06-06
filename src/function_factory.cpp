#include "function_factory.hpp"

#include <deal.II/base/numbers.h>

#include <cmath>
#include <stdexcept>

namespace
{
  using dealii::Function;
  using dealii::Point;
  using dealii::Tensor;

  static constexpr unsigned int dim = WaveSolver::dim;

  class ZeroScalarFunction : public Function<dim>
  {
  public:
    double
    value(const Point<dim> &, const unsigned int) const override
    {
      return 0.0;
    }

    Tensor<1, dim>
    gradient(const Point<dim> &, const unsigned int) const override
    {
      return {};
    }
  };

  class StandingWave5x5Function : public Function<dim>
  {
  public:
    explicit StandingWave5x5Function(const double wave_speed)
      : Function<dim>()
      , omega(wave_speed * dealii::numbers::PI * std::sqrt(2.0) / 5.0)
    {}

    double
    value(const Point<dim> &p, const unsigned int) const override
    {
      return std::sin((dealii::numbers::PI / 5.0) * p[0]) * std::sin((dealii::numbers::PI / 5.0) * p[1]) *
             std::cos(omega * this->get_time());
    }

    Tensor<1, dim>
    gradient(const Point<dim> &p, const unsigned int) const override
    {
      Tensor<1, dim> g;
      const double   factor = std::cos(omega * this->get_time());
      g[0]                  = (dealii::numbers::PI / 5.0) * std::cos((dealii::numbers::PI / 5.0) * p[0]) *
                              std::sin((dealii::numbers::PI / 5.0) * p[1]) * factor;
      g[1]                  = (dealii::numbers::PI / 5.0) * std::sin((dealii::numbers::PI / 5.0) * p[0]) *
                              std::cos((dealii::numbers::PI / 5.0) * p[1]) * factor;
      return g;
    }

  private:
    const double omega;
  };

  class StandingWave5x5VelocityFunction : public Function<dim>
  {
  public:
    explicit StandingWave5x5VelocityFunction(const double wave_speed)
      : Function<dim>()
      , omega(wave_speed * dealii::numbers::PI * std::sqrt(2.0) / 5.0)
    {}

    double
    value(const Point<dim> &p, const unsigned int) const override
    {
      return -omega * std::sin((dealii::numbers::PI / 5.0) * p[0]) * std::sin((dealii::numbers::PI / 5.0) * p[1]) *
             std::sin(omega * this->get_time());
    }

  private:
    const double omega;
  };

  class StandingWave5x5ForcingFunction : public Function<dim>
  {
  public:
    explicit StandingWave5x5ForcingFunction(const double wave_speed)
      : Function<dim>()
      , wave_speed(wave_speed)
      , omega(wave_speed * dealii::numbers::PI * std::sqrt(2.0) / 5.0)
    {}

    double
    value(const Point<dim> &p, const unsigned int) const override
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

  class GaussianPulseFunction : public Function<dim>
  {
  public:
    double
    value(const Point<dim> &p, const unsigned int) const override
    {
      const double dx = p[0] - 2.5;
      const double dy = p[1] - 2.5;
      return std::exp(-8.0 * (dx * dx + dy * dy));
    }
  };

} // namespace

std::shared_ptr<dealii::Function<WaveSolver::dim>>
make_named_function(const std::string &name, const double wave_speed)
{
  if (name == "zero" || name == "zero_dirichlet")
    return std::make_shared<ZeroScalarFunction>();
  if (name == "standing_wave_5x5" || name == "standing_wave_5x5_exact")
    return std::make_shared<StandingWave5x5Function>(wave_speed);
  if (name == "standing_wave_5x5_velocity")
    return std::make_shared<StandingWave5x5VelocityFunction>(wave_speed);
  if (name == "standing_wave_5x5_forcing")
    return std::make_shared<StandingWave5x5ForcingFunction>(wave_speed);
  if (name == "gaussian_pulse")
    return std::make_shared<GaussianPulseFunction>();

  throw std::runtime_error("Unknown function name: " + name);
}
