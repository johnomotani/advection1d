#include <cblas.h>
#include <cmath>
#include <iostream>

#include "../parameters.hxx"

#include "chebyshev.hxx"

Chebyshev::Chebyshev(const Parameters &parameters)
    : Nz_with_ghosts(parameters.Nz + 1), Nz(parameters.Nz), L(parameters.L),
      N(Nz_with_ghosts/* - 1*/), bc(stringToBC(parameters.bc)), z(createZValues()) {

  // Derivative coefficients using Cardinal functions for Chebyshev polynomials
  // on a grid including end points (Boyd Appendix F) x-coordinate of Boyd is in
  // range [-1, 1]. z-coordinate here is in range [0, L]. So z = L/2 * (x+1),
  // dz/dx = L/2
  //
  // i, j here are (N - i) (N - j) for Boyd. Conventions chosen here to have z increasing
  // with i.
  //
  // Note: not optimal to have if statements inside loop body, but this is only
  // done once during initialisation, so nicer to have clear code that looks
  // like eq. (F.45) of Boyd
  z_deriv_coefficients = createArray(Nz_with_ghosts * Nz_with_ghosts);
  const double dxdz = 2.0 / L;
  for (size_t j = 0; j < Nz_with_ghosts; ++j) {
    //const double x_j = cos(pi * double(N - j) / N);
    const double x_j = cos(pi * (2.0 * double(N - j) - 1.0) / (2.0 * N));
    for (size_t i = 0; i < Nz_with_ghosts; ++i) {
      //const double x_i = cos(pi * double(N - i) / N);
      const double x_i = cos(pi * (2.0 * double(N - i) - 1.0) / (2.0 * N));
      const auto ind = j * Nz_with_ghosts + i;
      //if ((N - i) == 0 and (N - j) == 0) {
      //  z_deriv_coefficients[ind] = dxdz * (1 + 2 * N * N) / 6.0;
      //} else if ((N - i) == N and (N - j) == N) {
      //  z_deriv_coefficients[ind] = -dxdz * (1 + 2 * N * N) / 6.0;
      //} else if ((N - i) == (N - j)) {
      //  z_deriv_coefficients[ind] = -dxdz * x_j / (2.0 * (1.0 - x_j * x_j));
      //} else {
      //  const double sign = ((N - i) + (N - j)) % 2 == 0 ? 1.0 : -1.0;
      //  const double p_i = ((N - i) == 0 or (N - i) == N) ? 2.0 : 1.0;
      //  const double p_j = ((N - j) == 0 or (N - j) == N) ? 2.0 : 1.0;
      //  z_deriv_coefficients[ind] = sign * p_i / (p_j * (x_i - x_j));
      //}
      if (i == j) {
        z_deriv_coefficients[ind] = 0.5 * x_j / (1.0 - x_j * x_j);
        if (i==N-1 and j==N-1) {
          std::cout<<"here "<<x_j<<" "<<0.5 * x_j / (1.0 - x_j * x_j)<<" "<<0.5 * x_j << " " << (1.0 - x_j * x_j)<<std::endl;
        }
      } else {
        const double sign = ((N - 1 - i) + (N - 1 - j)) % 2 == 0 ? 1.0 : -1.0;
        z_deriv_coefficients[ind] = sign * sqrt((1.0 - x_j * x_j) / (1.0 - x_i * x_i)) / (x_i - x_j);
        if (i==1 and j==0) {
          //std::cout<<"here "<<x_i<<" "<<x_j<<" "<<sign<<" "<<(1.0 - x_j * x_j)<<" "<<(1.0 - x_i * x_i)<<" "<<(x_i - x_j)<<std::endl;
        }
      }
    }
  }
}

Array initialdfdz(int N, double L, const Array &z) {
  auto f = createArray(N);
  for (size_t i = 0; i < N; ++i) {
    const double zhat = z[i] - 0.5 * L;
    f[i] = -32.0 * zhat * exp(-16.0 * zhat * zhat);
  }
  return f;
}

void Chebyshev::rhs(const double t, Array &f, Array &k) const {
  // Calculate -v*df/dz and store the result in k

  //applyBoundary(t, f);
  //f[N] = 0.0;

  // dfdz = alpha * z_deriv_coefficients.f + beta * dfdz
  cblas_dgemv(CblasColMajor,  // Matrix layout
              CblasNoTrans,   // Do not transpose matrix
              Nz_with_ghosts, // Number of rows
              Nz_with_ghosts, // Number of columns
              1.0,            // Scalar coefficient multiplying matrix
              z_deriv_coefficients.data(), // Derivative matrix coefficients
              Nz_with_ghosts, // Leading dimension of A (=number of rows)
              f.data(),       // Pointer to f's data
              1,              // Stride between elements of f
              0.0,            // Scalar beta multiplying dfdz on rhs
              k.data(),       // Pointer to dfdz's data
              1               // Stride between elements of dfdz
  );
  auto dfdz = initialdfdz(Nz_with_ghosts, L, z);
  for (size_t i = 0; i < Nz_with_ghosts; i++) {
    std::cout << i << "  " << z[i] << "  \t" << f[i] << "  \t" << dfdz[i] << "  \t" << k[i] << "  \t" << dfdz[i] - k[i] << std::endl;
  }
  std::cout << std::endl;
  for (size_t i = 0; i < Nz_with_ghosts; i++) {
    std::cout << i << "  ";
    for (size_t j = 0; j < Nz_with_ghosts; j++) {
      const auto ind = j * Nz_with_ghosts + i;
      std::cout << z_deriv_coefficients[ind] << "  \t";
    }
    std::cout << std::endl;
  }
  exit(0);

  for (size_t i = 1; i < Nz_with_ghosts; i++) {
    k[i] *= -v(t, i);
  }
}

double Chebyshev::v(const double t, const int i) const {
  return 0.1;
  // return 0.1 + 0.05 * sin(2.0 * pi * z[i] / L);
}

void Chebyshev::applyBoundary(const double t, Array &f) const {
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

double Chebyshev::fLower(const double t) const {
  return sin(t) + cos(0.9 * t) + 0.05 * t;
}

void Chebyshev::initialisef(Array &f) const {
  for (size_t i = 0; i < Nz_with_ghosts; ++i) {
    const double zhat = z[i] - 0.5 * L;
    f[i] = exp(-16.0 * zhat * zhat);
  }
}
