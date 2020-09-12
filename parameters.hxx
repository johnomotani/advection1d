#ifndef __PARAMETERS_H__
#define __PARAMETERS_H__

#include <string>

constexpr double pi = 3.1415926535897932;

struct Parameters {
  /// Order of polynomials
  int N = 32;

  /// Size of spatial domain
  double L = 1.0;

  /// internal time-step
  double dt = .01;

  /// time between outputs
  double t_out = 1.0;

  /// number of outputs
  int N_out = 100;

  /// Incoming boundary condition
  std::string bc = "periodic";

  /// Scheme to use for time-stepper
  std::string solver_type = "euler";

  /// Scheme to use for spatial discretisation
  std::string spatial_type = "chebyshevmatrix";
};

const Parameters createParameters();

#endif // __PARAMETERS_H__
