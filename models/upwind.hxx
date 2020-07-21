#ifndef __UPWIND_H__
#define __UPWIND_H__

#include "../array.hxx"
#include "model.hxx"

class Parameters;

/// Provides upwind-discretised right-hand-side function
/// rhs = - v * df/dz
/// rhs[i] = v(t, z[i]) * (f[i] - f[i-1]) / dz
class Upwind {
public:
  Upwind(const Parameters &parameters);

  void rhs(const double t, Array &f, Array &k) const;

  void applyBoundary(const double t, Array &f) const;

  void initialisef(Array &f) const;

  const size_t Nz_with_ghosts;

private:
  double v(const double t, const int i) const;

  double fLower(const double t) const;

  const size_t Nz;
  const double L;
  const double dz;

  const Array createZValues() const {
    auto z_values = createArray(Nz_with_ghosts);
    for (size_t i = 0; i < Nz_with_ghosts; ++i) {
      z_values[i] = i * dz;
    }
    return z_values;
  }

  const BC bc = BC::periodic;

public:
  const Array z;
};

#endif // __UPWIND_H__
