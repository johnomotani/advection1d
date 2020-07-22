#include "forwardeuler.hxx"

template <typename M>
ForwardEuler<M>::ForwardEuler(const Parameters &parameters, Output &output)
    : SolverBase<M>(parameters, output) {

  k = createArray(Nz);
}

template <typename M> void ForwardEuler<M>::updatef() {
  model.rhs(t, f, k);

  for (size_t i = 0; i < Nz; i++) {
    f[i] = f[i] + dt * k[i];
  }
}

INSTANTIATE_FOR_MODELS(ForwardEuler)
