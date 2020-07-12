#include "array.hxx"

class Parameters; // Forward-declare to avoid including header here

/// Provides upwind-discretised right-hand-side function
/// rhs = - v * df/dz
/// rhs[i] = v(t, z[i]) * (f[i] - f[i-1]) / dz
class Model {
public:
  Model(Parameters);

  void rhs(const double t, Array &f, Array &k) const;

  double v(const double t, const int i) const;

  void applyBoundary(const double t, Array &f) const;

  double fLower(const double t) const;

  void initialisef(Array &f) const;

private:
  const int Nz;
  const double L;
  const double dz;

  enum class BC { periodic, Dirichlet };
  const BC bc = BC::periodic;
  BC stringToBC(std::string input);
};
