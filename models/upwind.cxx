#include <cmath>

#include "../parameters.hxx"

#include "upwind.hxx"

Upwind::Upwind(const Parameters &parameters)
    : Nz_with_ghosts(parameters.Nz + 1), Nz(parameters.Nz), L(parameters.L),
      dz(parameters.L / parameters.Nz), bc(stringToBC(parameters.bc)),
      z(createZValues()) {}

void Upwind::rhs(const double t, Array &f, Array &k) const {
  // Calculate -v*df/dz and store the result in k

  applyBoundary(t, f);

  for (size_t i = 1; i < Nz_with_ghosts; i++) {
    k[i] = -v(t, i) * (f[i] - f[i - 1]) / dz;
  }
}

double Upwind::v(const double t, const int i) const {
  return 0.1;
  // return 0.1 + 0.05 * sin(2.0 * pi * (i - 1) / Nz);
}

void Upwind::applyBoundary(const double t, Array &f) const {
  switch (bc) {
  case BC::periodic:
    f[0] = f[Nz];
    break;
  case BC::Dirichlet:
    f[0] = fLower(t);
    break;
  default:
    throw "Unrecognised boundary condition";
  }
}

double Upwind::fLower(const double t) const {
  return sin(t) + cos(0.9 * t) + 0.05 * t;
}

void Upwind::initialisef(Array &f) const {
  for (size_t i = 0; i < Nz_with_ghosts; ++i) {
    const double zhat = (i - Nz / 2.0) / Nz;
    f[i] = exp(-16.0 * zhat * zhat);
  }
}
