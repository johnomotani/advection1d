#include "forwardeuler.hxx"

ForwardEuler::ForwardEuler(const Parameters &parameters, const Model &model,
                           Output &output)
    : Solver(parameters, model, output) {

  k = createArray(Nz_plus_1);
}

void ForwardEuler::updatef() {
  model.rhs(t, f, k);

  for (size_t i = 1; i < Nz_plus_1; i++) {
    f[i] = f[i] + dt * k[i];
  }
}
