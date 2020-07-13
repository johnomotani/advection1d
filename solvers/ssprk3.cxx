#include "ssprk3.hxx"

template<typename M>
SSPRK3<M>::SSPRK3(const Parameters &parameters, Output &output)
    : SolverBase<M>(parameters, output) {

  k1 = createArray(Nz_plus_1);
  k2 = createArray(Nz_plus_1);
  k3 = createArray(Nz_plus_1);
  f_temp = createArray(Nz_plus_1);
}

template<typename M>
void SSPRK3<M>::updatef() {
  // k1 = rhs(t[n], f[n])
  model.rhs(t, f, k1);

  // f[n] + dt*k1
  AequalBplussTimesC(f, dt, k1, f_temp);
  // k2 = rhs(t[n] + dt, f[n] + dt*k1
  model.rhs(t + dt, f_temp, k2);

  // f[n] + dt/2*k2
  AequalBplussTimesC(f, 0.5 * dt, k2, f_temp);
  AequalBplussTimesCPlusD(f, 0.25 * dt, k1, k2, f_temp);
  // k3 = rhs(t[n] + dt/2, f[n] + dt/4*(k1 + k2)
  model.rhs(t + 0.5 * dt, f_temp, k3);

  for (size_t i = 1; i < Nz_plus_1; i++) {
    f[i] = f[i] + 0.16666666666666666 * dt * (k1[i] + k2[i] + 4.0 * k3[i]);
  }
}

INSTANTIATE_FOR_MODELS(SSPRK3)
