#include "rk2.hxx"

template <typename M>
RK2<M>::RK2(const Parameters &parameters, Output &output)
    : SolverBase<M>(parameters, output) {

  k1 = createArray(Nz);
  k2 = createArray(Nz);
  f_temp = createArray(Nz);
}

template <typename M> void RK2<M>::updatef() {
  // k1 = rhs(t[n], f[n])
  model.rhs(t, f, k1);

  // f[n] + dt*k1
  AequalBplussTimesC(f, dt, k1, f_temp);
  // k2 = rhs(t[n] + dt, f[n] + dt*k1
  model.rhs(t + dt, f_temp, k2);

  for (size_t i = 0; i < Nz; i++) {
    f[i] = f[i] + 0.5 * dt * (k1[i] + k2[i]);
  }
}

INSTANTIATE_FOR_MODELS(RK2)
