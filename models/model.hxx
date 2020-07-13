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

#define FOR_MODEL(model, x, xtype)                                 \
  if (model == "upwind") {                                         \
    x(Upwind, xtype)                                               \
  } else if (model == "centred") {                                 \
    x(Centred, xtype)                                              \
  } else {                                                         \
    std::ostringstream message;                                    \
    message << "Unrecognised spatial_type " << model << std::endl; \
    throw std::runtime_error(message.str());                       \
  }

#define INSTANTIATE_FOR_MODELS(thisclass) \
  template class thisclass<Upwind>;       \
  template class thisclass<Centred>;      \

#endif // __MODEL_H__
