#ifndef __MODEL_H__
#define __MODEL_H__

#include <sstream>
#include <stdexcept>

#define FOR_MODEL(model, x, xtype)                                             \
  if (model == "chebyshevmatrix") {                                            \
    x(ChebyshevMatrix, xtype)                                                  \
  } else if (model == "chebyshevfft") {                                        \
    x(ChebyshevFFT, xtype)                                                     \
  } else if (model == "chebyshevfft_r2r") {                                    \
    x(ChebyshevFFT_r2r, xtype)                                                 \
  } else if (model == "chebyshevfft_r2c") {                                    \
    x(ChebyshevFFT_r2c, xtype)                                                 \
  } else {                                                                     \
    std::ostringstream message;                                                \
    message << "Unrecognised spatial_type " << model << std::endl;             \
    throw std::runtime_error(message.str());                                   \
  }

#define INSTANTIATE_FOR_MODELS(thisclass)                                      \
  template class thisclass<ChebyshevMatrix>;                                   \
  template class thisclass<ChebyshevFFT>;                                      \
  template class thisclass<ChebyshevFFT_r2r>;                                  \
  template class thisclass<ChebyshevFFT_r2c>;

enum class BC { periodic, Dirichlet };

inline BC stringToBC(std::string input) {
  if (input == "periodic") {
    return BC::periodic;
  } else if (input == "Dirichlet" or input == "dirichlet") {
    return BC::Dirichlet;
  } else {
    std::ostringstream message;
    message << "Unrecognised value " << input << " for boundary condition";
    throw std::runtime_error(message.str());
  }
}

#endif // __MODEL_H__
