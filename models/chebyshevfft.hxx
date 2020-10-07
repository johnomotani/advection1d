#ifndef __CHEBYSHEV_FFT_H__
#define __CHEBYSHEV_FFT_H__

#include <cmath>
#include <complex.h>
#include <fftw3.h>

#include "../array.hxx"
#include "model.hxx"

class Parameters;

/// Provides Chebyshev-discretised right-hand-side function
/// rhs = - v * df/dz
class ChebyshevFFT {
public:
  ChebyshevFFT(const Parameters &parameters);

  ~ChebyshevFFT() {
    fftw_destroy_plan(transform_plan);
    fftw_destroy_plan(inverse_plan);
    fftw_free(dct);
    fftw_free(doubled_f);
  }

  void rhs(const double t, Array &f, Array &k);

  void applyBoundary(const double t, Array &f) const;

  void initialisef(Array &f) const;

private:
  double v(const double t, const int i) const;

  double fLower(const double t) const;

  void dfdz(Array &f, Array &k);

  void applyDdtBoundary(const double t, Array &f) const;

  /// number of grid points, including end points
  const size_t Nz;

  /// 'physical' length of grid
  const double L;

  /// Matrix coefficients for d/dz operator
  Array z_deriv_coefficients;

  /// Maximum order of Chebyshev polynomials being used
  /// One less than number of grid points (including end points)
  const int N;

  const Array createZValues() const {
    auto z_values = createArray(Nz);

    for (size_t i = 0; i < Nz; ++i) {
      /// rescaled from x in [-1, 1] to z in [0, L]
      z_values[i] = (cos(pi * double(N - i) / N) + 1.0) * 0.5 * L;
    }

    return z_values;
  }

  const BC bc = BC::periodic;

  fftw_complex *dct;
  fftw_complex *doubled_f;
  fftw_plan transform_plan;
  fftw_plan inverse_plan;

public:
  const Array z;
};

#endif // __CHEBYSHEV_FFT_H__
