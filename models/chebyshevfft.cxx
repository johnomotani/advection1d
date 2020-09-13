#include <cblas.h>
#include <chrono>
#include <cmath>
#include <iostream>

#include "../parameters.hxx"

#include "chebyshevfft.hxx"

ChebyshevFFT::ChebyshevFFT(const Parameters &parameters)
    : Nz(parameters.N + 1), L(parameters.L), N(parameters.N),
      bc(stringToBC(parameters.bc)), dct(createArray(Nz)), z(createZValues()) {

  auto temp = createArray(Nz);
  auto start = std::chrono::high_resolution_clock::now();

  // Create plan for discrete cosine transform (Real-Even-Discrete-Fourier-Transform)
  // Done here so we don't include the timing of the first plan-creation in the main loop.
  // Some solvers store several state vectors, so need to create a new plan each time rhs is
  // called, but new plans are identical to the inital one, so this should be cheap.
  auto temp_plan = fftw_plan_r2r_1d(Nz, &temp[0], &dct[0], FFTW_REDFT00, FFTW_EXHAUSTIVE);
  //auto temp_plan = fftw_plan_r2r_1d(Nz, &temp[0], &dct[0], FFTW_REDFT00, FFTW_PATIENT);
  fftw_destroy_plan(temp_plan);

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cout << std::endl
            << "FFTW plan creation took " << elapsed.count() << "s" << std::endl;

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
  auto start = std::chrono::high_resolution_clock::now();
  auto end = std::chrono::high_resolution_clock::now();
  static std::chrono::duration<double> planning=start-end;
  static std::chrono::duration<double> execute=start-end;
  static std::chrono::duration<double> loop=start-end;

  // Transform to Chebyshev-coefficient space
  // The output of the transform is dct and
  // dct[k]/N = c[k]*a[k]
  // where c[k] are the coefficients defined under Boyd's (A.15) and a[k] are the
  // Chebyshev coefficients defined in Boyd's (2.76)
  start = std::chrono::high_resolution_clock::now();
  auto forward_plan = fftw_plan_r2r_1d(Nz, &f[0], &dct[0], FFTW_REDFT00, FFTW_WISDOM_ONLY);
  end = std::chrono::high_resolution_clock::now();
  planning += end - start;
  start = std::chrono::high_resolution_clock::now();
  fftw_execute(forward_plan);
  fftw_destroy_plan(forward_plan);
  end = std::chrono::high_resolution_clock::now();
  execute += end - start;

  start = std::chrono::high_resolution_clock::now();
  // Calculate derivative in Chebyshev space using Boyd's (A.15).
  // Calculate in-place in dct to avoid allocating another array
  dct[N] = 0.0;

  auto temp = dct[N - 1];
  dct[N - 1] = 0.0;

  // contributions to prefactor:
  //   1/N  convert input dct[k] to c[k]*a[k]
  //   1/2  convert c[k]*a1[k] to input for inverse transform that outputs df/dx[j]
  //   2/L  dx/dz converts from df/dx (derivative on [-1,1] grid) to df/dz (derivative
  //        on physical grid [0, L]
  // Due to definitions of FFTW3's REDFT00 transform and Boyd's definition (2.76) of the
  // spectral coefficients, the c[k] factor in the recursion equation for the derivative
  // actually just provides the conversion of Boyd's a[0] to the 0'th coefficient of the
  // input to the REDFT00.
  const auto prefactor = 1.0 / (double(N) * L);

  for (size_t k = N - 1; k > 0; --k) {
    // temp is input-dct[k]
    const auto newval = prefactor * (2.0 * double(k) * temp + dct[k + 1]);
    temp = dct[k - 1];
    dct[k - 1] = newval;
  }
  end = std::chrono::high_resolution_clock::now();
  loop += end - start;

  // The inverse transform would transform c[k]*a1[k] to 2*f[j]
  start = std::chrono::high_resolution_clock::now();
  auto inverse_plan = fftw_plan_r2r_1d(Nz, &dct[0], &k[0], FFTW_REDFT00, FFTW_WISDOM_ONLY);
  end = std::chrono::high_resolution_clock::now();
  planning += end - start;
  start = std::chrono::high_resolution_clock::now();
  fftw_execute(inverse_plan);
  fftw_destroy_plan(inverse_plan);
  end = std::chrono::high_resolution_clock::now();
  execute += end - start;

  std::cout<<"planning="<<planning.count()<<" execute="<<execute.count()<<" loop="<<loop.count()<<std::endl;
}
