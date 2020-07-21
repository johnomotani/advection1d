#ifndef __CHEBYSHEV_H__
#define __CHEBYSHEV_H__

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

  const size_t Nz_with_ghosts;

private:
  double v(const double t, const int i) const;

  double fLower(const double t) const;

  const size_t Nz;
  const double L;
  const Array z;
  Array z_deriv_coefficients;

  const Array getZValues() const {
    auto z_values = createArray(Nz_with_ghosts);
    return z_values;
  }

  const BC bc = BC::periodic;
};

#endif // __CHEBYSHEV_H__
