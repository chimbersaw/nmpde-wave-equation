#pragma once

class WaveSolver;

class ThetaTimeIntegrator
{
public:
  ThetaTimeIntegrator(WaveSolver &solver, double theta);

  void
  run();

private:
  WaveSolver &solver;
  double      theta;
};
