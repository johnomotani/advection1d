#ifndef __CENTRED_H__
#define __CENTRED_H__

#include "../array.hxx"
#include "model.hxx"

class Parameters;

/// Provides centred-discretised right-hand-side function
/// rhs = - v * df/dz
/// rhs[i] = v(t, z[i]) * (f[i+1] - f[i-1]) / (2*dz)
class Centred {
public:
  Centred(const Parameters &parameters);

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

  const BC bc = BC::periodic;
};

#endif // __CENTRED_H__
