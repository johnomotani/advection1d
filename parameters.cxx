#include "submodules/toml11/toml.hpp"

#include "parameters.hxx"

#ifdef GET_OPTION
#error "GET_OPTION already defined. This should not happen"
#endif
#define GET_OPTION(x) result.x = toml::find_or(input, "x", result.x);

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
  GET_OPTION(Nz);
  GET_OPTION(L);
  GET_OPTION(dt);
  GET_OPTION(t_out);
  GET_OPTION(N_out);
  GET_OPTION(bc);
  GET_OPTION(solver_type);

  return result;
}

#undef GET_OPTION
