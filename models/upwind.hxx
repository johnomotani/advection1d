#ifndef __UPWIND_H__
#define __UPWIND_H__

#include "model.hxx"

/// Provides upwind-discretised right-hand-side function
/// rhs = - v * df/dz
/// rhs[i] = v(t, z[i]) * (f[i] - f[i-1]) / dz
class Upwind : public Model {
public:
  Upwind(const Parameters &parameters);
  ~Upwind() final = default;

  void rhs(const double t, Array &f, Array &k) const final;
};

#endif // __UPWIND_H__
