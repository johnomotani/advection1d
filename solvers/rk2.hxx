#ifndef __RK2_H__
#define __RK2_H__

#include "solver.hxx"

/// 2nd order Runge-Kutta time-stepper (Heun's method)
/// https://en.wikipedia.org/wiki/Runge-Kutta_methods
template <typename M> class RK2 : public SolverBase<M> {
public:
  RK2(const Parameters &parameters, Output &output);
  ~RK2() final = default;

protected:
  using SolverBase<M>::model;
  using SolverBase<M>::t;
  using SolverBase<M>::dt;
  using SolverBase<M>::t_out;
  using SolverBase<M>::N_out;
  using SolverBase<M>::Nz;
  using SolverBase<M>::f;

  void updatef() final;

private:
  // rhs vectors
  Array k1;
  Array k2;

  // temporary state vector
  Array f_temp;

  void AequalBplussTimesC(const Array &B, const double s, const Array &C,
                          Array &A) {
    for (size_t i = 0; i < Nz; i++) {
      A[i] = B[i] + s * C[i];
    }
  }
};

#endif // __RK2_H__
