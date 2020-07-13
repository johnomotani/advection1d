#include <chrono>
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

  auto start = std::chrono::high_resolution_clock::now();
  solver->run();
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cout << std::endl
            << "Main solve took " << elapsed.count() << "s" << std::endl;

  return 0;
}
