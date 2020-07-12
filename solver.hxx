#include "array.hxx"

class Model;
class Parameters;
class Output;

/// Forward-Euler time-stepper
/// f[n+1] = f[n] + dt * rhs
class Solver {
public:
  Solver(const Parameters parameters, const Model model, Output &output);

  void run();

private:
  void updatef();

  void writeOutput();

  Model model;
  Output &output;

  double t;
  double dt;
  double t_out;
  int N_out;
  size_t Nz_plus_1;

  // state vector
  Array f;

  // rhs vector
  Array k;
};
