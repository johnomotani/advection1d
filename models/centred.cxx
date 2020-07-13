#include "centred.hxx"

Centred::Centred(const Parameters &parameters) : Model(parameters) {}

void Centred::rhs(const double t, Array &f, Array &k) const {
  // Calculate -v*df/dz and store the result in k

  applyBoundary(t, f);

  for (size_t i = 1; i < f.size(); i++) {
    k[i] = -v(t, i) * (f[i + 1] - f[i - 1]) / (2.0 * dz);
  }
}
