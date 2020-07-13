#ifndef __SOLVER_H__
#define __SOLVER_H__

#include <memory>

#include "../array.hxx"
#include "../models/model.hxx"
#include "../output.hxx"

#include "../models/centred.hxx"
#include "../models/upwind.hxx"

class Parameters;

/// Generic base class for explicit time-steppers
class Solver {
public:
  virtual ~Solver() = default;

  virtual void run() = 0;
};

template <typename M> class SolverBase : public Solver {
public:
  SolverBase(const Parameters &parameters, Output &output);
  virtual ~SolverBase() = default;

  void run() final;

protected:
  virtual void updatef() = 0;

  void writeOutput();

  const M model;

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
                                     Output &output);

#endif // __SOLVER_H__
