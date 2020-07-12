#include <cmath>
#include <sstream>
#include <stdexcept>

#include "parameters.hxx"

#include "model.hxx"

Model::Model(Parameters parameters)
    : Nz(parameters.Nz), L(parameters.L), dz(parameters.L / parameters.Nz),
      bc(stringToBC(parameters.bc)) {}

void Model::rhs(const double t, Array &f, Array &k) const {
  // Calculate -v*df/dz and store the result in k

  applyBoundary(t, f);

  for (size_t i = 1; i < f.size(); i++) {
    k[i] = -v(t, i) * (f[i] - f[i - 1]) / dz;
  }
}

double Model::v(const double t, const int i) const { return 0.1; }

void Model::applyBoundary(const double t, Array &f) const {
  switch (bc) {
  case BC::periodic:
    f[0] = f[Nz];
    break;
  case BC::Dirichlet:
    f[0] = fLower(t);
    break;
  default:
    throw "Unrecognised boundary condition";
  }
}

double Model::fLower(const double t) const {
  return sin(t) + cos(0.9 * t) + 0.05 * t;
}

void Model::initialisef(Array &f) const {
  for (size_t i = 0; i < f.size(); ++i) {
    const double zhat = (i - Nz / 2.0) / Nz;
    f[i] = exp(-16.0 * zhat * zhat);
  }
}

Model::BC Model::stringToBC(std::string input) {
  if (input == "periodic") {
    return BC::periodic;
  } else if (input == "Dirichlet" or input == "diriclet") {
    return BC::Dirichlet;
  } else {
    std::ostringstream message;
    message << "Unrecognised value " << input << " for boundary condition";
    throw std::runtime_error(message.str());
  }
}
