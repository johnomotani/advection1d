#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "../parameters.hxx"

#include "forwardeuler.hxx"
#include "rk4.hxx"
#include "solver.hxx"
#include "ssprk3.hxx"

template <typename M>
SolverBase<M>::SolverBase(const Parameters &parameters, Output &output)
    : model(parameters), t(0.0), dt(parameters.dt), t_out(parameters.t_out),
      N_out(parameters.N_out), Nz_plus_1(parameters.Nz + 1), output(output) {

  f = createArray(Nz_plus_1);

  model.initialisef(f);
}

template <typename M> void SolverBase<M>::run() {
  // write initial state
  writeOutput();

  for (int output_step = 0; output_step < N_out; ++output_step) {
    // std::cout << output_step << " " << t << std::endl;

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

template <typename M> void SolverBase<M>::writeOutput() {
  output.writeStep(t, f);
}

// Explicit instantiation of SolverBase for each model
template class SolverBase<Upwind>;
template class SolverBase<Centred>;

#define RETURN_SOLVER(model, solver)                                           \
  return std::unique_ptr<Solver>(new solver<model>(parameters, output));

std::unique_ptr<Solver> createSolver(const Parameters &parameters,
                                     Output &output) {
  const auto solver_type = parameters.solver_type;
  const auto spatial_type = parameters.spatial_type;
  if (solver_type == "euler") {
    FOR_MODEL(spatial_type, RETURN_SOLVER, ForwardEuler)
  } else if (solver_type == "rk4") {
    FOR_MODEL(spatial_type, RETURN_SOLVER, RK4)
  } else if (solver_type == "ssprk3") {
    FOR_MODEL(spatial_type, RETURN_SOLVER, SSPRK3)
  } else {
    std::ostringstream message;
    message << "Unrecognised time-step scheme option " << solver_type
            << std::endl;
    throw std::runtime_error(message.str());
  }
}
