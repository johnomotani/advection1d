#ifndef __CHEBYSHEV_H__
#define __CHEBYSHEV_H__

#include <cmath>

#include "../array.hxx"
#include "model.hxx"

class Parameters;

/// Provides Chebyshev-discretised right-hand-side function
/// rhs = - v * df/dz
class Chebyshev {
public:
  Chebyshev(const Parameters &parameters);

  void rhs(const double t, Array &f, Array &k) const;

  void applyBoundary(const double t, Array &f) const;

  void initialisef(Array &f) const;

private:
  double v(const double t, const int i) const;

  double fLower(const double t) const;

  /// number of grid points, including end points
  const size_t Nz;

  /// 'physical' length of grid
  const double L;

  /// Matrix coefficients for d/dz operator
  Array z_deriv_coefficients;

  /// Maximum order of Chebyshev polynomials being used
  /// One less than number of grid points (including end points)
  const int N;

  const Array createZValues() const {
    auto z_values = createArray(Nz);

    for (size_t i = 0; i < Nz; ++i) {
      /// rescaled from x in [-1, 1] to z in [0, L]
      z_values[i] = (cos(pi * double(N - i) / N) + 1.0) * 0.5 * L;
    }

    return z_values;
  }

  const BC bc = BC::periodic;

public:
  const Array z;
};

#endif // __CHEBYSHEV_H__
