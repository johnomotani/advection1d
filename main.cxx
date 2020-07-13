#include <iostream>

#include "models/model.hxx"
#include "output.hxx"
#include "parameters.hxx"
#include "solvers/solver.hxx"

int main() {

  std::cout << "running 1d advection..." << std::endl;

  const auto parameters = createParameters();

  Output output;

  auto solver = createSolver(parameters, output);

  solver->run();

  return 0;
}
