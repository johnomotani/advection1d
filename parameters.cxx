#include "submodules/toml11/toml.hpp"

#include "parameters.hxx"

const Parameters createParameters() {
  Parameters result;

  toml::value input;
  try {
    input = toml::parse("input.toml");
  } catch (const std::runtime_error) {
    // No input file - keep default values
    return result;
  }

  // Get values, using values set in parameters.hxx as defaults
  result.Nz = toml::find_or(input, "Nz", result.Nz);
  result.L = toml::find_or(input, "L", result.L);
  result.dt = toml::find_or(input, "dt", result.dt);
  result.t_out = toml::find_or(input, "t_out", result.t_out);
  result.N_out = toml::find_or(input, "N_out", result.N_out);

  return result;
}
