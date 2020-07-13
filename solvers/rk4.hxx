#ifndef __RK4_H__
#define __RK4_H__

#include "solver.hxx"

/// 4th order Runge-Kutta time-stepper
/// https://en.wikipedia.org/wiki/Runge-Kutta_methods
class RK4 : public Solver {
public:
  RK4(const Parameters &parameters, const Model &model, Output &output);
  ~RK4() final = default;

protected:
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
