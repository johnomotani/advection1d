#include <cmath>

#include "../parameters.hxx"

#include "centred.hxx"

Centred::Centred(const Parameters &parameters)
    : Nz_with_ghosts(parameters.Nz + 2), Nz(parameters.Nz), L(parameters.L),
      dz(parameters.L / parameters.Nz), bc(stringToBC(parameters.bc)) {}

void Centred::rhs(const double t, Array &f, Array &k) const {
  // Calculate -v*df/dz and store the result in k

  applyBoundary(t, f);

  for (size_t i = 1; i <= Nz; i++) {
    k[i] = -v(t, i) * (f[i + 1] - f[i - 1]) / (2.0 * dz);
  }
}

double Centred::v(const double t, const int i) const {
  return 0.1;
  // return 0.1 + 0.05 * sin(2.0 * pi * (i - 1) / Nz);
}

void Centred::applyBoundary(const double t, Array &f) const {
  switch (bc) {
  case BC::periodic:
    f[0] = f[Nz];
    f[Nz + 1] = f[1];
    break;
  case BC::Dirichlet:
    f[0] = fLower(t);
    // Neumann outgoing boundary condition
    f[Nz + 1] = f[Nz];
    break;
  default:
    throw "Unrecognised boundary condition";
  }
}

double Centred::fLower(const double t) const {
  return sin(t) + cos(0.9 * t) + 0.05 * t;
}

void Centred::initialisef(Array &f) const {
  for (size_t i = 0; i < f.size(); ++i) {
    const double zhat = (i - Nz / 2.0) / Nz;
    f[i] = exp(-16.0 * zhat * zhat);
  }
}
