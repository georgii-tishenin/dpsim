// SPDX-License-Identifier: MPL-2.0
// Three-phase EMT response to a load step or a cleared fault.
#define IEEE_SSN_BUS_COUNT 9
#include "IEEE_SSN_BenchmarkOutput.h"

int main(int argc, char **argv) {
  using namespace IEEEBenchmark;
  try {
    const String name = "IEEE9_SSN_Accuracy";
    const Options options(argc, argv, name, 4000, "validate", "network",
                          "load");
    if (!options.validate() || options.steps * timeStep <= .10)
      throw std::invalid_argument("Accuracy requires validate and >2000 steps");
    const auto preparationStart = BenchmarkClock::now();
    auto network = makeNetwork(options.formulation);

    // Both events act at a retained generator bus, outside the passive group.
    // The load's nominal power is 20% of total demand; actual power follows V²/R.
    // The fault uses 0.001 pu resistance and is cleared after three 60 Hz cycles.
    const bool fault = options.event == "fault";
    const Real eventTime = .05, clearingTime = .10;
    const Real voltage = std::abs(
        network.powerflow.node<SP::SimNode>("N" + std::to_string(eventBus))
            ->singleVoltage());
    const Real resistance = fault ? .001 * zBase(busFor(eventBus))
                                  : voltage * voltage / loadStepPower();
    auto event = EMT::Ph3::Switch::make("Event");
    event->setParameters(parameter3(1e12), parameter3(resistance), false);
    event->connect({EMT::SimNode::GND, network.terminal(eventBus)});
    network.system.addComponent(event);

    Simulation sim(name, Logger::Level::off);
    sim.setSystem(network.system);
    sim.setDomain(Domain::EMT);
    sim.setSolverType(Solver::Type::MNA);
    sim.setDirectLinearSolverImplementation(DirectLinearSolverImpl::KLU);
    DirectLinearSolverConfiguration solver;
    solver.setPartialRefactorizationMethod(
        PARTIAL_REFACTORIZATION_METHOD::NO_PARTIAL_REFACTORIZATION);
    solver.setStaticPreconditioner(false);
    sim.setDirectLinearSolverConfiguration(solver);
    sim.setTimeStep(timeStep);
    sim.setFinalTime(options.steps * timeStep);
    sim.addEvent(SwitchEvent3Ph::make(eventTime, event, true));
    if (fault)
      sim.addEvent(SwitchEvent3Ph::make(clearingTime, event, false));

    Measurements measurements(sim, network, options, name, "KLU",
                              preparationStart, event);
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
    nlohmann::json parameters = {{"event_type", options.event},
                                 {"event_time_s", eventTime},
                                 {"event_bus_line_voltage_V", voltage}};
    if (fault)
      parameters.update({{"fault_clear_time_s", clearingTime},
                         {"fault_resistance_ohm", resistance}});
    else
      parameters.update({{"load_step_power_W", loadStepPower()},
                         {"load_step_resistance_ohm", resistance}});
    measurements.write(sim, parameters);
    return 0;
  } catch (const std::exception &error) {
    SPDLOG_LOGGER_ERROR(Logger::get("ssn-benchmark"), "{}", error.what());
    return 1;
  }
}
