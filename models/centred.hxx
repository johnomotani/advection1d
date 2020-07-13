#ifndef __CENTRED_H__
#define __CENTRED_H__

#include "model.hxx"

/// Provides centred-discretised right-hand-side function
/// rhs = - v * df/dz
/// rhs[i] = v(t, z[i]) * (f[i+1] - f[i-1]) / (2*dz)
class Centred : public Model {
public:
  Centred(const Parameters &parameters);
  ~Centred() final = default;

  void rhs(const double t, Array &f, Array &k) const final;
};

#endif // __CENTRED_H__
