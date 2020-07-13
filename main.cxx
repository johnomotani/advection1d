#include <iostream>

#include "models/model.hxx"
#include "output.hxx"
#include "parameters.hxx"
#include "solvers/solver.hxx"

int main() {

  std::cout << "running 1d advection..." << std::endl;

  const auto parameters = createParameters();

  auto model = createModel(parameters);

  Output output;

  auto solver = createSolver(parameters, model.get(), output);

  solver->run();

  return 0;
}
