#include <cmath>
#include <cblas.h>

#include "../parameters.hxx"

#include "chebyshev.hxx"

Chebyshev::Chebyshev(const Parameters &parameters)
    : Nz_with_ghosts(parameters.Nz + 1), Nz(parameters.Nz), L(parameters.L),
      z(getZValues()), bc(stringToBC(parameters.bc)) {

  // Initialise as identity for trying out blaspp
  // Column-major order
  z_deriv_coefficients = createArray(Nz_with_ghosts*Nz_with_ghosts);
  for (size_t j = 0; j < Nz_with_ghosts; ++j) {
    for (size_t i = 0; i < Nz_with_ghosts; ++i) {
      auto ind = j*Nz_with_ghosts + i;
      z_deriv_coefficients[ind] = i == j ? 1.0 : 0.0;
    }
  }
}

void Chebyshev::rhs(const double t, Array &f, Array &k) const {
  // Calculate -v*df/dz and store the result in k

  applyBoundary(t, f);

  // dfdz = alpha * z_deriv_coefficients.f + beta * dfdz
  cblas_dgemv(
      CblasColMajor,                // Matrix layout
      CblasNoTrans,                 // Do not transpose matrix
      Nz_with_ghosts,               // Number of rows
      Nz_with_ghosts,               // Number of columns
      1.0,                          // Scalar coefficient multiplying matrix
      z_deriv_coefficients.data(),  // Derivative matrix coefficients
      Nz_with_ghosts,               // Leading dimension of A (=number of rows)
      f.data(),                     // Pointer to f's data
      1,                            // Stride between elements of f
      0.0,                          // Scalar beta multiplying dfdz on rhs
      k.data(),                     // Pointer to dfdz's data
      1                             // Stride between elements of dfdz
  );

  for (size_t i = 1; i < Nz_with_ghosts; i++) {
    k[i] *= -v(t, i);
  }
}

double Chebyshev::v(const double t, const int i) const {
  return 0.1;
  // return 0.1 + 0.05 * sin(2.0 * pi * (i - 1) / Nz);
}

void Chebyshev::applyBoundary(const double t, Array &f) const {
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

double Chebyshev::fLower(const double t) const {
  return sin(t) + cos(0.9 * t) + 0.05 * t;
}

void Chebyshev::initialisef(Array &f) const {
  for (size_t i = 0; i < Nz_with_ghosts; ++i) {
    const double zhat = (i - Nz / 2.0) / Nz;
    f[i] = exp(-16.0 * zhat * zhat);
  }
}
