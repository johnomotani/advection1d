#ifndef __MODEL_H__
#define __MODEL_H__

#include <memory>
#include <string>

#include "../array.hxx"

class Parameters; // Forward-declare to avoid including header here

/// Base class for implementations providing right-hand-side function
class Model {
public:
  Model(const Parameters &parameters);
  virtual ~Model() = default;

  virtual void rhs(const double t, Array &f, Array &k) const = 0;

  double v(const double t, const int i) const;

  void applyBoundary(const double t, Array &f) const;

  double fLower(const double t) const;

  void initialisef(Array &f) const;

protected:
  const int Nz;
  const double L;
  const double dz;

  enum class BC { periodic, Dirichlet };
  const BC bc = BC::periodic;
  BC stringToBC(std::string input);
};

std::unique_ptr<const Model> createModel(const Parameters &parameters);

#endif // __MODEL_H__
