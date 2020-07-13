#include "forwardeuler.hxx"

template<typename M>
ForwardEuler<M>::ForwardEuler(const Parameters &parameters, Output &output)
    : SolverBase<M>(parameters, output) {

  k = createArray(Nz_plus_1);
}

template<typename M>
void ForwardEuler<M>::updatef() {
  model.rhs(t, f, k);

  for (size_t i = 1; i < Nz_plus_1; i++) {
    f[i] = f[i] + dt * k[i];
  }
}

INSTANTIATE_FOR_MODELS(ForwardEuler)
