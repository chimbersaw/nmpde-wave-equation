#pragma once

class WaveSolver;

class NewmarkTimeIntegrator
{
public:
  NewmarkTimeIntegrator(WaveSolver &solver, double beta, double gamma);

  void
  run();

private:
  WaveSolver &solver;
  double      beta;
  double      gamma;
};
