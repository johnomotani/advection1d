#include <cmath>
#include <iostream>

#include "model.hxx"
#include "output.hxx"
#include "parameters.hxx"

#include "solver.hxx"

Solver::Solver(const Parameters parameters, const Model model, Output &output)
    : model(model), output(output), t(0.0), dt(parameters.dt),
      t_out(parameters.t_out), N_out(parameters.N_out),
      Nz_plus_1(parameters.Nz + 1) {

  auto Nz = parameters.Nz;

  f = createArray(Nz + 1);
  k = createArray(Nz + 1);

  model.initialisef(f);
}

void Solver::run() {
  // write initial state
  writeOutput();

  for (int output_step = 0; output_step < N_out; ++output_step) {
    std::cout << output_step << " " << t << std::endl;

    const int N_internal = round(((output_step + 1) * t_out - t) / dt);

    for (int n = 0; n < N_internal; ++n) {
      updatef();
      t += dt;
    }

    model.applyBoundary(t, f);
    writeOutput();
  }

  std::cout << N_out << " " << t << std::endl;
}

void Solver::updatef() {
  model.rhs(t, f, k);

  for (size_t i = 1; i < Nz_plus_1; i++) {
    f[i] = f[i] + dt * k[i];
  }
}

void Solver::writeOutput() { output.writeStep(t, f); }
