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

private:
  double v(const double t, const int i) const;

  double fLower(const double t) const;

  const int Nz;
  const double L;
  const double dz;

  const BC bc = BC::periodic;
};

#endif // __UPWIND_H__
