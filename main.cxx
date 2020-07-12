#include <iostream>

#include "model.hxx"
#include "output.hxx"
#include "parameters.hxx"
#include "solvers/solver.hxx"

int main() {

  std::cout << "running 1d advection..." << std::endl;

  const auto parameters = createParameters();

  const Model model(parameters);

  Output output;

  auto solver = createSolver(parameters, model, output);

  solver->run();

  return 0;
}
