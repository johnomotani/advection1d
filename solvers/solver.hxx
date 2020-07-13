#ifndef __SOLVER_H__
#define __SOLVER_H__

#include <memory>

#include "../array.hxx"
#include "../models/model.hxx"
#include "../output.hxx"

class Parameters;

/// Generic base class for explicit time-steppers
class Solver {
public:
  Solver(const Parameters &parameters, const Model *const model,
         Output &output);
  virtual ~Solver() = default;

  void run();

protected:
  virtual void updatef() = 0;

  void writeOutput();

  const Model *const model;

  double t;
  double dt;
  double t_out;
  int N_out;
  size_t Nz_plus_1;

  // state vector
  Array f;

private:
  Output &output;
};

std::unique_ptr<Solver> createSolver(const Parameters &parameters,
                                     const Model *const model, Output &output);

#endif // __SOLVER_H__
