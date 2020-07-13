#ifndef __FORWARDEULER_H__
#define __FORWARDEULER_H__

#include "solver.hxx"

/// Forward-Euler time-stepper
/// f[n+1] = f[n] + dt * rhs
template <typename M> class ForwardEuler : public SolverBase<M> {
public:
  ForwardEuler(const Parameters &parameters, Output &output);
  ~ForwardEuler() final = default;

protected:
  using SolverBase<M>::model;
  using SolverBase<M>::t;
  using SolverBase<M>::dt;
  using SolverBase<M>::t_out;
  using SolverBase<M>::N_out;
  using SolverBase<M>::Nz_plus_1;
  using SolverBase<M>::f;

  void updatef() final;

private:
  // rhs vector
  Array k;
};

#endif // __FORWARDEULER_H__
