#include "function_factory.hpp"

#include <deal.II/base/numbers.h>

#include <cmath>
#include <stdexcept>

namespace
{
  using dealii::Function;
  using dealii::Point;
  using dealii::Tensor;

  static constexpr unsigned int dim                = WaveSolver::dim;
  static constexpr double       domain_length      = 5.0;
  static constexpr double       circle_center      = 2.5;
  static constexpr double       circle_radius      = 2.5;
  static constexpr double       donut_inner_radius = 0.9;
  static constexpr double       donut_outer_radius = 2.5;

  bool
  is_left_boundary(const Point<dim> &p)
  {
    return std::abs(p[0]) < 1e-10;
  }

  bool
  is_right_boundary(const Point<dim> &p)
  {
    return std::abs(p[0] - domain_length) < 1e-10;
  }

  double
  transverse_profile(const Point<dim> &p)
  {
    return std::sin(dealii::numbers::PI * p[1] / domain_length);
  }

  double
  circle_radius_squared(const Point<dim> &p)
  {
    const double dx = p[0] - circle_center;
    const double dy = p[1] - circle_center;
    return dx * dx + dy * dy;
  }

  double
  circle_spatial_part(const Point<dim> &p)
  {
    const double q = circle_radius_squared(p) / (circle_radius * circle_radius);
    return (1.0 - q) * (1.0 - q);
  }

  Tensor<1, dim>
  circle_spatial_gradient(const Point<dim> &p)
  {
    Tensor<1, dim> g;
    const double   r2     = circle_radius_squared(p);
    const double   factor = -4.0 * (1.0 - r2 / (circle_radius * circle_radius)) / (circle_radius * circle_radius);

    g[0] = factor * (p[0] - circle_center);
    g[1] = factor * (p[1] - circle_center);
    return g;
  }

  double
  circle_spatial_laplace(const Point<dim> &p)
  {
    const double r2 = circle_radius_squared(p);
    return -8.0 / (circle_radius * circle_radius) +
           16.0 * r2 / (circle_radius * circle_radius * circle_radius * circle_radius);
  }

  double
  donut_spatial_part(const Point<dim> &p)
  {
    const double dx    = p[0] - circle_center;
    const double dy    = p[1] - circle_center;
    const double r     = std::sqrt(dx * dx + dy * dy);
    const double theta = std::atan2(dy, dx);
    const double s     = (r - donut_inner_radius) / (donut_outer_radius - donut_inner_radius);
    return std::sin(dealii::numbers::PI * s) * std::cos(3.0 * theta);
  }

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

  class CircleStandingWaveFunction : public Function<dim>
  {
  public:
    explicit CircleStandingWaveFunction(const double wave_speed)
      : Function<dim>()
      , omega(wave_speed * dealii::numbers::PI / circle_radius)
    {}

    double
    value(const Point<dim> &p, const unsigned int) const override
    {
      return circle_spatial_part(p) * std::cos(omega * this->get_time());
    }

    Tensor<1, dim>
    gradient(const Point<dim> &p, const unsigned int) const override
    {
      Tensor<1, dim> g      = circle_spatial_gradient(p);
      const double   factor = std::cos(omega * this->get_time());
      for (unsigned int d = 0; d < dim; ++d)
        g[d] *= factor;
      return g;
    }

  protected:
    const double omega;
  };

  class CircleStandingWaveVelocityFunction : public CircleStandingWaveFunction
  {
  public:
    explicit CircleStandingWaveVelocityFunction(const double wave_speed)
      : CircleStandingWaveFunction(wave_speed)
    {}

    double
    value(const Point<dim> &p, const unsigned int) const override
    {
      return -omega * circle_spatial_part(p) * std::sin(omega * this->get_time());
    }

    Tensor<1, dim>
    gradient(const Point<dim> &p, const unsigned int) const override
    {
      Tensor<1, dim> g      = circle_spatial_gradient(p);
      const double   factor = -omega * std::sin(omega * this->get_time());
      for (unsigned int d = 0; d < dim; ++d)
        g[d] *= factor;
      return g;
    }
  };

  class CircleStandingWaveForcingFunction : public Function<dim>
  {
  public:
    explicit CircleStandingWaveForcingFunction(const double wave_speed)
      : Function<dim>()
      , wave_speed(wave_speed)
      , omega(wave_speed * dealii::numbers::PI / circle_radius)
    {}

    double
    value(const Point<dim> &p, const unsigned int) const override
    {
      const double temporal = std::cos(omega * this->get_time());
      const double utt      = -omega * omega * circle_spatial_part(p) * temporal;
      const double laplace  = circle_spatial_laplace(p) * temporal;
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

  class RectangleDrivenLeftBoundaryFunction : public Function<dim>
  {
  public:
    double
    value(const Point<dim> &p, const unsigned int) const override
    {
      if (!is_left_boundary(p))
        return 0.0;

      const double omega = 2.0 * dealii::numbers::PI;
      return 0.5 * transverse_profile(p) * std::sin(omega * this->get_time());
    }
  };

  class BoundaryPulseLeftFunction : public Function<dim>
  {
  public:
    double
    value(const Point<dim> &p, const unsigned int) const override
    {
      if (!is_left_boundary(p))
        return 0.0;

      const double pulse_duration = 0.5;
      const double t              = this->get_time();
      if (t < 0.0 || t > pulse_duration)
        return 0.0;

      const double dy       = p[1] - 0.5 * domain_length;
      const double width    = 0.6;
      const double spatial  = transverse_profile(p) * std::exp(-0.5 * dy * dy / (width * width));
      const double temporal = std::pow(std::sin(dealii::numbers::PI * t / pulse_duration), 2.0);
      return spatial * temporal;
    }
  };

  class TwoSidedDrivenBoundaryFunction : public Function<dim>
  {
  public:
    double
    value(const Point<dim> &p, const unsigned int) const override
    {
      double side_factor = 0.0;
      if (is_left_boundary(p))
        side_factor = 1.0;
      else if (is_right_boundary(p))
        side_factor = -1.0;
      else
        return 0.0;

      const double omega = 2.0 * dealii::numbers::PI;
      return 0.5 * side_factor * transverse_profile(p) * std::sin(omega * this->get_time());
    }
  };

  class DonutAngularModeFunction : public Function<dim>
  {
  public:
    double
    value(const Point<dim> &p, const unsigned int) const override
    {
      return donut_spatial_part(p);
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
  if (name == "circle_standing_wave")
    return std::make_shared<CircleStandingWaveFunction>(wave_speed);
  if (name == "circle_standing_wave_velocity")
    return std::make_shared<CircleStandingWaveVelocityFunction>(wave_speed);
  if (name == "circle_standing_wave_forcing")
    return std::make_shared<CircleStandingWaveForcingFunction>(wave_speed);
  if (name == "gaussian_pulse")
    return std::make_shared<GaussianPulseFunction>();
  if (name == "rectangle_driven_left_boundary")
    return std::make_shared<RectangleDrivenLeftBoundaryFunction>();
  if (name == "boundary_pulse_left")
    return std::make_shared<BoundaryPulseLeftFunction>();
  if (name == "two_sided_driven_boundary")
    return std::make_shared<TwoSidedDrivenBoundaryFunction>();
  if (name == "donut_angular_mode")
    return std::make_shared<DonutAngularModeFunction>();

  throw std::runtime_error("Unknown function name: " + name);
}
