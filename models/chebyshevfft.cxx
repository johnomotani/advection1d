#include <cblas.h>
#include <cmath>

#include "../parameters.hxx"

#include "chebyshevfft.hxx"

ChebyshevFFT::ChebyshevFFT(const Parameters &parameters)
    : Nz(parameters.N + 1), L(parameters.L), N(parameters.N),
      bc(stringToBC(parameters.bc)), dct(createArray(Nz)), z(createZValues()) {

  // Create plan for discrete cosine transform (Real-Even-Discrete-Fourier-Transform)
  auto temp = createArray(Nz);
  transform_plan = fftw_plan_r2r_1d(Nz, &temp[0], &dct[0], FFTW_REDFT00, FFTW_EXHAUSTIVE);
}

void ChebyshevFFT::rhs(const double t, Array &f, Array &k) {
  // Calculate -v*df/dz and store the result in k

  applyBoundary(t, f);

  // k = df/dz
  dfdz(f, k);

  for (size_t i = 1; i < Nz; i++) {
    k[i] *= -v(t, i);
  }

  applyDdtBoundary(t, k);
}

double ChebyshevFFT::v(const double t, const int i) const {
  return 0.1;
  // return 0.1 + 0.05 * sin(2.0 * pi * z[i] / L);
}

void ChebyshevFFT::applyBoundary(const double t, Array &f) const {
  switch (bc) {
  case BC::periodic:
    f[0] = f[Nz - 1];
    break;
  case BC::Dirichlet:
    f[0] = fLower(t);
    break;
  default:
    throw "Unrecognised boundary condition";
  }
}

void ChebyshevFFT::applyDdtBoundary(const double t, Array &k) const {
  switch (bc) {
  case BC::periodic:
  {
    //const auto mean = 0.5 * (k[0] + k[Nz - 1]);
    //std::cout << k[0] << " " << k[Nz - 1] << " " << mean << std::endl;
    //k[0] = mean;
    //k[Nz - 1] = mean;
    break;
  }
  case BC::Dirichlet:
    break;
  default:
    throw "Unrecognised boundary condition";
  }
}

double ChebyshevFFT::fLower(const double t) const {
  return 0.0;
  //return sin(t) + cos(0.9 * t) + 0.05 * t;
}

void ChebyshevFFT::initialisef(Array &f) const {
  for (size_t i = 0; i < Nz; ++i) {
    const double zhat = z[i] - 0.5 * L;
    f[i] = exp(-64.0 * zhat * zhat / (L * L));
  }
}

/// Calculate df/dz by
/// Chebyshev transform -> spectral derivative -> inverse Chebyshev Transform
void ChebyshevFFT::dfdz(Array &f, Array& k) {
  // Transform to Chebyshev-coefficient space
  // The output of the transform is dct and
  // dct[k]/N = c[k]*a[k]
  // where c[k] are the coefficients defined under Boyd's (A.15) and a[k] are the
  // Chebyshev coefficients defined in Boyd's (2.76)
  fftw_execute_r2r(transform_plan, &f[0], &dct[0]);

  // Calculate derivative in Chebyshev space using Boyd's (A.15).
  // Calculate in-place in dct to avoid allocating another array
  dct[N] = 0.0;

  auto temp = dct[N - 1];
  dct[N - 1] = 0.0;

  for (size_t j = N - 1; j > 0; --j) {
    // temp is input-dct[j]
    const auto newval = 2.0 * double(j) * temp + dct[j + 1];
    temp = dct[j - 1];
    dct[j - 1] = newval;
  }

  // contributions to prefactor:
  //   1/N  convert input dct[k] to c[k]*a[k]
  //   1/2  convert c[k]*a1[k] to input for inverse transform that outputs df/dx[j]
  //   -1   account for sign in x = -cos(theta)
  //   2/L  dx/dz converts from df/dx (derivative on [-1,1] grid) to df/dz (derivative
  //        on physical grid [0, L]
  // Due to definitions of FFTW3's REDFT00 transform and Boyd's definition (2.76) of the
  // spectral coefficients, the c[k] factor in the recursion equation for the derivative
  // actually just provides the conversion of Boyd's a[0] to the 0'th coefficient of the
  // input to the REDFT00.
  const auto prefactor = -1.0 / (double(N) * L);
  for (auto &value : dct) {
    value *= prefactor;
  }

  // The inverse transform would transform c[k]*a1[k] to 2*f[j]
  fftw_execute_r2r(transform_plan, &dct[0], &k[0]);
}
