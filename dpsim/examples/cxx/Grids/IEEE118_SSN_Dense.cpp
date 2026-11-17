// SPDX-License-Identifier: MPL-2.0
// Three-phase EMT simulation with selectable SSN grouping and DenseLU.
#define IEEE_SSN_BUS_COUNT 118
#include "IEEE_SSN_BenchmarkOutput.h"

int main(int argc, char **argv) {
  using namespace IEEEBenchmark;
  try {
    const String name = "IEEE118_SSN_Dense";
    const Options options(argc, argv, name, 2000);
    const auto preparationStart = BenchmarkClock::now();
    auto network = makeNetwork(options.formulation);

    // Identical operating point, machines, passive circuit and physical outputs.
    // The performance runs have no switching event.
    Simulation sim(name, Logger::Level::off);
    sim.setSystem(network.system);
    sim.setDomain(Domain::EMT);
    sim.setSolverType(Solver::Type::MNA);
    sim.setDirectLinearSolverImplementation(DirectLinearSolverImpl::DenseLU);
    DirectLinearSolverConfiguration solver;
    solver.setPartialRefactorizationMethod(
        PARTIAL_REFACTORIZATION_METHOD::NO_PARTIAL_REFACTORIZATION);
    solver.setStaticPreconditioner(false);
    sim.setDirectLinearSolverConfiguration(solver);
    sim.setTimeStep(timeStep);
    sim.setFinalTime(options.steps * timeStep);

    // Timing runs publish outputs in memory. Separate profile runs measure
    // electrical calculation stages; validation runs also write trajectories.
    Measurements measurements(sim, network, options, name, "DenseLU",
                              preparationStart);
    sim.initialize();
    measurements.bindOutputs();
    sim.start();
    measurements.beginLoop();
    for (size_t step = 1; step <= options.steps; ++step) {
      const Real time = sim.time();
      sim.step();
      measurements.observe(step, time);
    }
    measurements.endLoop();
    sim.stop();
    measurements.write(sim, {{"events", "none"}});
    return 0;
  } catch (const std::exception &error) {
    SPDLOG_LOGGER_ERROR(Logger::get("ssn-benchmark"), "{}", error.what());
    return 1;
  }
}
