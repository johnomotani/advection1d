#include "rk4.hxx"

template <typename M>
RK4<M>::RK4(const Parameters &parameters, Output &output)
    : SolverBase<M>(parameters, output) {

  k1 = createArray(Nz_with_ghosts);
  k2 = createArray(Nz_with_ghosts);
  k3 = createArray(Nz_with_ghosts);
  k4 = createArray(Nz_with_ghosts);
  f_temp = createArray(Nz_with_ghosts);
}

template <typename M> void RK4<M>::updatef() {
  // k1 = rhs(t[n], f[n])
  model.rhs(t, f, k1);

  // f[n] + dt/2*k1
  AequalBplussTimesC(f, 0.5 * dt, k1, f_temp);
  // k2 = rhs(t[n] + dt/2, f[n] + dt/2*k1
  model.rhs(t + 0.5 * dt, f_temp, k2);

  // f[n] + dt/2*k2
  AequalBplussTimesC(f, 0.5 * dt, k2, f_temp);
  // k3 = rhs(t[n] + dt/2, f[n] + dt/2*k2
  model.rhs(t + 0.5 * dt, f_temp, k3);

  // f[n] + dt*k3
  AequalBplussTimesC(f, dt, k3, f_temp);
  // k4 = rhs(t[n] + dt, f[n] + dt*k3
  model.rhs(t + dt, f_temp, k4);

  for (size_t i = 1; i < Nz_plus_1; i++) {
    f[i] = f[i] + 0.16666666666666666 * dt *
                      (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
  }
}

INSTANTIATE_FOR_MODELS(RK4)
