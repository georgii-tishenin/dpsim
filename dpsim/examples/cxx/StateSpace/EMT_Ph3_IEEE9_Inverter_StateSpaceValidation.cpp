// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#define DPSIM_IEEE9_INVERTER_MIX_LIBRARY
#include "../Circuits/EMT_Ph3_IEEE9_SSN_InverterMix.cpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <optional>
#include <vector>

namespace {

void writeModalResults(const DPsim::StateSpaceModalAnalysis &modal,
                       const std::filesystem::path &outputDirectory) {
  std::filesystem::create_directories(outputDirectory);

  std::ofstream modes(outputDirectory / "modes.csv");
  modes << std::setprecision(std::numeric_limits<CPS::Real>::max_digits10)
        << "mode,z_real,z_imag,lambda_real,lambda_imag,frequency_hz\n";
  const auto &z = modal.getDiscreteEigenvalues();
  const auto &lambda = modal.getContinuousEigenvalues();
  for (Eigen::Index idx = 0; idx < lambda.rows(); ++idx) {
    modes << idx << ',' << z(idx).real() << ',' << z(idx).imag() << ','
          << lambda(idx).real() << ',' << lambda(idx).imag() << ','
          << std::abs(lambda(idx).imag()) / (2.0 * PI) << '\n';
  }

  std::ofstream participation(outputDirectory / "participation.csv");
  participation
      << std::setprecision(std::numeric_limits<CPS::Real>::max_digits10)
      << "mode,state,state_name,p_real,p_imag,p_abs\n";
  const auto &factors = modal.getParticipationFactors();
  const auto &stateNames = modal.getStateNames();
  for (Eigen::Index mode = 0; mode < factors.cols(); ++mode) {
    for (Eigen::Index state = 0; state < factors.rows(); ++state) {
      participation << mode << ',' << state << ',' << stateNames[state] << ','
                    << factors(state, mode).real() << ','
                    << factors(state, mode).imag() << ','
                    << std::abs(factors(state, mode)) << '\n';
    }
  }

  std::ofstream summary(outputDirectory / "summary.csv");
  summary << std::setprecision(std::numeric_limits<CPS::Real>::max_digits10)
          << "state_count,auxiliary_reduction_residual,"
             "auxiliary_reduction_pole_error,zero_sequence_coupling_residual\n"
          << modal.getStateNames().size() << ','
          << modal.getAuxiliaryReductionResidual() << ','
          << modal.getAuxiliaryReductionPoleError() << ','
          << modal.getZeroSequenceCouplingResidual() << '\n';
}

struct OperatingPointSignal {
  String name;
  String unit;
  Real base;
  CPS::Attribute<Real>::Ptr attribute;
  Real initial = 0.0;
  std::optional<Real> preEvent;
  Real final = 0.0;
  Real finalPeriodMin = std::numeric_limits<Real>::infinity();
  Real finalPeriodMax = -std::numeric_limits<Real>::infinity();

  void sampleInitial() { initial = **attribute; }
  void samplePreEvent() { preEvent = **attribute; }
  void sampleFinalPeriod() {
    const Real value = **attribute;
    finalPeriodMin = std::min(finalPeriodMin, value);
    finalPeriodMax = std::max(finalPeriodMax, value);
  }
  void sampleFinal() { final = **attribute; }
};

void writeOperatingPointDiagnostics(
    const std::vector<OperatingPointSignal> &signals,
    const std::filesystem::path &outputDirectory) {
  std::filesystem::create_directories(outputDirectory);
  std::ofstream output(outputDirectory / "operating_point.csv");
  output << std::setprecision(std::numeric_limits<CPS::Real>::max_digits10)
         << "signal,unit,base,initial,pre_event,final,final_period_min,"
            "final_period_max,initial_to_final_pu,final_period_span_pu\n";
  for (const auto &signal : signals) {
    const Real scale = std::max(std::abs(signal.base),
                                std::numeric_limits<Real>::epsilon());
    output << signal.name << ',' << signal.unit << ',' << signal.base << ','
           << signal.initial << ',';
    if (signal.preEvent)
      output << *signal.preEvent;
    output << ',' << signal.final << ',' << signal.finalPeriodMin << ','
           << signal.finalPeriodMax << ','
           << (signal.final - signal.initial) / scale << ','
           << (signal.finalPeriodMax - signal.finalPeriodMin) / scale << '\n';
  }
}

} // namespace

int main(int argc, char *argv[]) {
  DPsim::CommandLineArgs args(
      argc, argv, "EMT_Ph3_IEEE9_Inverter_StateSpaceValidation", 50e-6, 0.1,
      60, -1, CPS::Logger::Level::info, CPS::Logger::Level::off, false, false,
      false, CPS::Domain::EMT);

  // Defaults specific to this validation; command-line options remain able to
  // override each value.
  args.options.try_emplace("network_equivalent", "true");
  args.options.try_emplace("transformer_magnetizing", "true");
  args.options.try_emplace("transformer_p0_pu", "0.001");
  args.options.try_emplace("transformer_q0_pu", "0.01");
  args.options.try_emplace("gfl_formulation", "variable");
  args.options.try_emplace("gfl_current_cross_coupling", "true");
  // The network is grounded through its sources and loads, and the power-flow
  // model neglects line leakage. Do not add the PiLine legacy 1 uS grounding
  // conductances only on the EMT side.
  args.options.try_emplace("line_default_conductance", "false");
  args.options.try_emplace("line96_breaker", "false");
  // 0.15 ohm is approximately 0.093 pu on the 18 kV / 200 MVA GFM base. The virtual
  // resistance damps the very-slow grid-connected voltage/synchronization
  // mode without changing the physical network impedance.
  args.options.try_emplace("gfm_rv", "0.15");

  const String scenario =
      args.options.find("scenario") != args.options.end()
          ? args.getOptionString("scenario")
          : "baseline";
  const Real eventTime =
      args.options.find("event_time_s") != args.options.end()
          ? args.getOptionReal("event_time_s")
          : 0.4 * args.duration;
  const Real loadStepPowerPu =
      args.options.find("load_step_p_pu") != args.options.end()
          ? args.getOptionReal("load_step_p_pu")
          : 0.1;
  const Real pulseDuration =
      args.options.find("pulse_duration_s") != args.options.end()
          ? args.getOptionReal("pulse_duration_s")
          : 0.02;
  const String disturbanceBus =
      args.options.find("disturbance_bus") != args.options.end()
          ? args.getOptionString("disturbance_bus")
          : "BUS8";

  CPS::Logger::setLogDir("./logs/" + args.name);
  const Bool logTimeDomain =
      args.options.find("log") != args.options.end() &&
      args.getOptionBool("log");

  std::shared_ptr<DPsim::DataLoggerInterface> logger;
  if (logTimeDomain) {
    std::filesystem::path filename =
        "./logs/" + args.name + "/" + args.name + ".csv";
    logger = DPsim::RealTimeDataLogger::make(filename, args.duration,
                                             args.timeStep);
  }

  auto system = buildTopology(args, logger);

  const auto gen2 = system.component<CPS::EMT::Ph3::SSN_GFM>("GEN2");
  const auto gen3 = system.component<CPS::IdentifiedObject>("GEN3");
  if (!gen2 || !gen3)
    throw std::logic_error(
        "IEEE 9-bus validation requires GEN2 (GFM) and GEN3 (GFL).");
  std::vector<OperatingPointSignal> operatingPointSignals{
      {"GEN2.p", "W", 200.0e6, gen2->attributeTyped<Real>("p_inst")},
      {"GEN2.q", "var", 200.0e6, gen2->attributeTyped<Real>("q_inst")},
      {"GEN2.omega", "rad/s", 2.0 * PI * args.sysFreq,
       gen2->attributeTyped<Real>("omega_gfm")},
      {"GEN3.p", "W", 100.0e6, gen3->attributeTyped<Real>("p_inst")},
      {"GEN3.q", "var", 100.0e6, gen3->attributeTyped<Real>("q_inst")},
      {"GEN3.omega", "rad/s", 2.0 * PI * args.sysFreq,
       gen3->attributeTyped<Real>("omega_pll")},
  };

  std::vector<std::shared_ptr<DPsim::SwitchEvent3Ph>> scenarioEvents;
  if (scenario == "load_pickup" || scenario == "load_pulse") {
    if (!(eventTime >= 0.0 && eventTime < args.duration))
      throw std::invalid_argument(
          "event_time_s must be nonnegative and precede the simulation end.");
    if (!(loadStepPowerPu > 0.0))
      throw std::invalid_argument("load_step_p_pu must be positive.");

    auto disturbanceNode =
        system.node<CPS::SimNode<Real>>(disturbanceBus);
    if (!disturbanceNode)
      throw std::invalid_argument("Unknown disturbance_bus: " +
                                  disturbanceBus);
    const Real loadResistance =
        std::norm(disturbanceNode->initialSingleVoltage()) /
        (loadStepPowerPu * 100.0e6);
    auto loadSwitch = CPS::EMT::Ph3::Switch::make(
        disturbanceBus + "_LOAD_DISTURBANCE", CPS::Logger::Level::off);
    loadSwitch->setParameters(CPS::Matrix::Identity(3, 3) * 1.0e9,
                              CPS::Matrix::Identity(3, 3) * loadResistance);
    loadSwitch->openSwitch();
    system.addComponent(loadSwitch);
    system.connectComponentToNodes<Real>(
        loadSwitch, {CPS::EMT::SimNode::GND, disturbanceNode});
    scenarioEvents.push_back(
        DPsim::SwitchEvent3Ph::make(eventTime, loadSwitch, true));
    if (scenario == "load_pulse") {
      if (!(pulseDuration > args.timeStep) ||
          eventTime + pulseDuration >= args.duration) {
        throw std::invalid_argument(
            "pulse_duration_s must exceed one step and end before the "
            "simulation duration.");
      }
      scenarioEvents.push_back(DPsim::SwitchEvent3Ph::make(
          eventTime + pulseDuration, loadSwitch, false));
    }
  } else if (scenario == "line_outage") {
    if (!(eventTime > args.timeStep && eventTime < args.duration))
      throw std::invalid_argument(
          "A line-outage event must occur after the first time step and "
          "before the simulation end.");
    const auto lineBreaker =
        system.component<CPS::EMT::Ph3::Switch>("LINE96_BREAKER");
    if (!lineBreaker)
      throw std::invalid_argument(
          "scenario=line_outage requires line96_breaker=true.");
    scenarioEvents.push_back(
        DPsim::SwitchEvent3Ph::make(eventTime, lineBreaker, false));
  } else if (scenario != "baseline") {
    throw std::invalid_argument("Unknown scenario: " + scenario);
  }

  DPsim::Simulation simulation(args.name, args);
  simulation.setSystem(system);
  simulation.setDomain(CPS::Domain::EMT);
  simulation.doSystemMatrixRecomputation(true);
  simulation.doInitFromNodesAndTerminals(true);
  simulation.doStateSpaceExtraction(true);
  if (logger)
    simulation.addLogger(logger);
  for (const auto &event : scenarioEvents)
    simulation.addEvent(event);

  const auto writeSnapshot = [&](const String &label) {
    DPsim::StateSpaceModalAnalysis modal(simulation.getStateSpaceExtractor());
    modal.setAnalysisFrame(DPsim::StateSpaceAnalysisFrame::GlobalDQ0);
    modal.setGlobalDq0Frame(2.0 * PI * args.sysFreq);
    modal.setExcludeDecoupledZeroSequenceStates(true);
    if (args.getOptionString("gfl_formulation") == "split") {
      modal.setPoleMapping(DPsim::StateSpacePoleMapping::Logarithmic);
    } else {
      modal.setReduceAuxiliaryStates(true);
    }
    modal.update();
    writeModalResults(modal, std::filesystem::path("./logs") / args.name /
                                 "modal" / label);
    CPS::Logger::get(args.name)->info(
        "Wrote {} modal states for snapshot '{}' at t={:.6f} s "
        "(zero-sequence coupling residual {:.3e}).",
        modal.getStateNames().size(), label, simulation.time(),
        modal.getZeroSequenceCouplingResidual());
  };

  simulation.initialize();
  for (auto &signal : operatingPointSignals)
    signal.sampleInitial();
  simulation.start();
  // A snapshot at exactly t=0 precedes the first valid extraction timestamp.
  // Zero-time disturbances use the separately generated baseline snapshot.
  Bool preEventWritten = scenarioEvents.empty() || eventTime <= 0.0;
  while (simulation.time() + 0.5 * args.timeStep < args.duration) {
    if (!scenarioEvents.empty() && !preEventWritten &&
        simulation.time() + 0.5 * args.timeStep >= eventTime) {
      for (auto &signal : operatingPointSignals)
        signal.samplePreEvent();
      writeSnapshot("pre_event");
      preEventWritten = true;
    }
    simulation.step();
    if (simulation.time() >= args.duration - 1.0 / args.sysFreq)
      for (auto &signal : operatingPointSignals)
        signal.sampleFinalPeriod();
  }
  for (auto &signal : operatingPointSignals)
    signal.sampleFinal();
  writeSnapshot(!scenarioEvents.empty() ? "post_event" : "baseline");
  writeOperatingPointDiagnostics(operatingPointSignals,
                                 std::filesystem::path("./logs") / args.name);
  simulation.stop();
}
