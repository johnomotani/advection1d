#ifndef __SSPRK3_H__
#define __SSPRK3_H__

#include "solver.hxx"

/// 3rd order, strong-stability-preserving Runge-Kutta time-stepper
/// https://en.wikipedia.org/wiki/List_of_Runge-Kutta_methods#Third-order_Strong_Stability_Preserving_Runge-Kutta_(SSPRK3)
template <typename M> class SSPRK3 : public SolverBase<M> {
public:
  SSPRK3(const Parameters &parameters, Output &output);
  ~SSPRK3() final = default;

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
  Array k3;

  // temporary state vector
  Array f_temp;

  void AequalBplussTimesC(const Array &B, const double s, const Array &C,
                          Array &A) {
    for (size_t i = 0; i < Nz; i++) {
      A[i] = B[i] + s * C[i];
    }
  }

  void AequalBplussTimesCPlusD(const Array &B, const double s, const Array &C,
                               const Array &D, Array &A) {
    for (size_t i = 0; i < Nz; i++) {
      A[i] = B[i] + s * (C[i] + D[i]);
    }
  }
};

#endif // __SSPRK3_H__
