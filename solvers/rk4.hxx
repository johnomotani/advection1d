#ifndef __RK4_H__
#define __RK4_H__

#include "solver.hxx"

/// 4th order Runge-Kutta time-stepper
/// https://en.wikipedia.org/wiki/Runge-Kutta_methods
template <typename M> class RK4 : public SolverBase<M> {
public:
  RK4(const Parameters &parameters, Output &output);
  ~RK4() final = default;

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
  // rhs vectors
  Array k1;
  Array k2;
  Array k3;
  Array k4;

  // temporary state vector
  Array f_temp;

  void AequalBplussTimesC(const Array &B, const double s, const Array &C,
                          Array &A) {
    for (size_t i = 1; i < Nz_plus_1; i++) {
      A[i] = B[i] + s * C[i];
    }
  }
};

#endif // __RK4_H__
