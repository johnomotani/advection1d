#include "upwind.hxx"

Upwind::Upwind(const Parameters &parameters) : Model(parameters) {}

void Upwind::rhs(const double t, Array &f, Array &k) const {
  // Calculate -v*df/dz and store the result in k

  applyBoundary(t, f);

  for (size_t i = 1; i < f.size(); i++) {
    k[i] = -v(t, i) * (f[i] - f[i - 1]) / dz;
  }
}
