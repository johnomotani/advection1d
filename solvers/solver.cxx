#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "../parameters.hxx"

#include "forwardeuler.hxx"
#include "rk4.hxx"
#include "solver.hxx"
#include "ssprk3.hxx"

Solver::Solver(const Parameters &parameters, const Model *const model,
               Output &output)
    : model(model), t(0.0), dt(parameters.dt), t_out(parameters.t_out),
      N_out(parameters.N_out), Nz_plus_1(parameters.Nz + 1), output(output) {

  f = createArray(Nz_plus_1);

  model->initialisef(f);
}

void Solver::run() {
  // write initial state
  writeOutput();

  for (int output_step = 0; output_step < N_out; ++output_step) {
    // std::cout << output_step << " " << t << std::endl;

    const int N_internal = round(((output_step + 1) * t_out - t) / dt);

    for (int n = 0; n < N_internal; ++n) {
      updatef();
      t += dt;
    }

    model->applyBoundary(t, f);
    writeOutput();
  }

  std::cout << N_out << " " << t << std::endl;
}

void Solver::writeOutput() { output.writeStep(t, f); }

std::unique_ptr<Solver> createSolver(const Parameters &parameters,
                                     const Model *const model, Output &output) {
  const auto type = parameters.solver_type;
  if (type == "euler") {
    return std::unique_ptr<Solver>(new ForwardEuler(parameters, model, output));
  } else if (type == "rk4") {
    return std::unique_ptr<Solver>(new RK4(parameters, model, output));
  } else if (type == "ssprk3") {
    return std::unique_ptr<Solver>(new SSPRK3(parameters, model, output));
  } else {
    std::ostringstream message;
    message << "Unrecognised time-step scheme option " << type << std::endl;
    throw std::runtime_error(message.str());
  }
}
