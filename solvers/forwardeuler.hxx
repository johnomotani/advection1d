#ifndef __FORWARDEULER_H__
#define __FORWARDEULER_H__

#include "solver.hxx"

/// Forward-Euler time-stepper
/// f[n+1] = f[n] + dt * rhs
class ForwardEuler : public Solver {
public:
  ForwardEuler(const Parameters &parameters, const Model *const model,
               Output &output);
  ~ForwardEuler() final = default;

protected:
  void updatef() final;

private:
  // rhs vector
  Array k;
};

#endif // __FORWARDEULER_H__
