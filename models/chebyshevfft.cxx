#include <cblas.h>
#include <cmath>

#include "../parameters.hxx"

#include "chebyshevfft.hxx"

ChebyshevFFT::ChebyshevFFT(const Parameters &parameters)
    : Nz(parameters.N + 1), L(parameters.L), N(parameters.N),
      bc(stringToBC(parameters.bc)), z(createZValues()) {

  dct = fftw_alloc_complex(2*N);
  doubled_f = fftw_alloc_complex(2*N);

  // Create plan for discrete cosine transform (Real-Even-Discrete-Fourier-Transform)
  transform_plan = fftw_plan_dft_1d(2*Nz - 2, doubled_f, dct, FFTW_FORWARD, FFTW_EXHAUSTIVE);
  inverse_plan = fftw_plan_dft_1d(2*Nz - 2, dct, doubled_f, FFTW_BACKWARD, FFTW_EXHAUSTIVE);
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
    f[i] = exp(-128.0 * zhat * zhat / (L * L));
  }
}

/// Calculate df/dz by
/// Chebyshev transform -> spectral derivative -> inverse Chebyshev Transform
void ChebyshevFFT::dfdz(Array &f, Array& k) {
  // create even input array on 0->2pi
  for (size_t i = 0; i < N; i++) {
    doubled_f[i][0] = f[i];
    doubled_f[i][1] = 0.0;
    doubled_f[2*N - 1 - i][0] = f[i + 1];
    doubled_f[2*N - 1 - i][1] = 0.0;
  }

  // Transform to Chebyshev-coefficient space
  // The output of the transform is dct and
  // dct[k]/N = c[k]*a[k]
  // where c[k] are the coefficients defined under Boyd's (A.15) and a[k] are the
  // Chebyshev coefficients defined in Boyd's (2.76)
  fftw_execute_dft(transform_plan, doubled_f, dct);

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

  // Calculate derivative in Chebyshev space using Boyd's (A.15).
  // Calculate in-place in dct to avoid allocating another array
  const auto before_loop_temp = prefactor * dct[N][0];
  dct[N][0] = 0.0;

  auto temp = prefactor * dct[N - 1][0];
  dct[N - 1][0] = 2.0 * double(N) * before_loop_temp;
  dct[N + 1][0] = dct[N - 1][0];

  for (size_t j = N - 1; j > 1; --j) {
    // temp is input-dct[j]
    const auto newval = 2.0 * double(j) * temp + dct[j + 1][0];
    temp = prefactor * dct[j - 1][0];
    dct[j - 1][0] = newval;
    dct[2*N - j + 1][0] = newval;
  }
  dct[0][0] = 2.0 * temp + dct[2][0];

  // The inverse transform would transform c[k]*a1[k] to 2*f[j]
  fftw_execute_dft(inverse_plan, dct, doubled_f);

  // copy result into k
  for (size_t i = 0; i < Nz; i++) {
    k[i] = doubled_f[i][0];
  }
}
