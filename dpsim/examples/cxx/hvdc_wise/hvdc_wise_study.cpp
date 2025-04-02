#include <DPsim.h>
using namespace DPsim;
using namespace CPS;

struct SimulationParameters {
  double timeStep = 1e-5;
  double eventTime = 0.1;
  double finalTime = 0.2; 
};

struct PowerSystemParameters {
  double frequency = 50;
  double voltage = 10e3;
  double infeed_resistance = 1;
  double infeed_inductance = 0.5;
  double line1_resistance = 1;
  double line1_inductance = 0.1;
  double line1_capacitance = 3e-6;
  double line2_resistance = 2;
  double line2_inductance = 0.2;
  double line2_capacitance = 2e-6;
  double circuit_breaker_closed_resistance = 1e-4;
  double circuit_breaker_open_resistance = 1e6;
  double load_resistance1 = 100;
  double load_resistance2 = 50;
};

void simulate_EMT(const SimulationParameters &simParams, const PowerSystemParameters &psParams) {
  String simName = "EMT_simulation";
  Logger::setLogDir("logs/" + simName);

  // nodes
  auto node1 = EMT::SimNode::make("node1", PhaseType::ABC);
  auto node2 = EMT::SimNode::make("node2", PhaseType::ABC);
  auto node3 = EMT::SimNode::make("node3", PhaseType::ABC);
  auto node4 = EMT::SimNode::make("node4", PhaseType::ABC);
  auto node5 = EMT::SimNode::make("node5", PhaseType::ABC);
  auto node6 = EMT::SimNode::make("node6", PhaseType::ABC);
  auto node7 = EMT::SimNode::make("node7", PhaseType::ABC);

  // components
  auto infeed_source = EMT::Ph3::VoltageSource::make("infeed_source");
  infeed_source->setParameters(CPS::Math::singlePhaseVariableToThreePhase(CPS::Math::polar(psParams.voltage, 0.0)), psParams.frequency);
  auto infeed_impedance = EMT::Ph3::PiLine::make("infeed_impedance");
  infeed_impedance->setParameters(CPS::Math::singlePhaseParameterToThreePhase(psParams.infeed_resistance), CPS::Math::singlePhaseParameterToThreePhase(psParams.infeed_inductance),
                                  CPS::Math::singlePhaseParameterToThreePhase(0));
  auto converter1 = EMT::Ph3::VoltageSource::make("converter1");
  converter1->setParameters(CPS::Math::singlePhaseVariableToThreePhase(CPS::Math::polar(psParams.voltage, 0.0)), psParams.frequency);
  auto line1 = EMT::Ph3::PiLine::make("line1");
  line1->setParameters(CPS::Math::singlePhaseParameterToThreePhase(psParams.line1_resistance), CPS::Math::singlePhaseParameterToThreePhase(psParams.line1_inductance),
                       CPS::Math::singlePhaseParameterToThreePhase(psParams.line1_capacitance));
  auto converter2 = EMT::Ph3::VoltageSource::make("converter2");
  converter2->setParameters(CPS::Math::singlePhaseVariableToThreePhase(CPS::Math::polar(psParams.voltage, 0.0)), psParams.frequency);
  auto line2 = EMT::Ph3::PiLine::make("line2");
  line2->setParameters(CPS::Math::singlePhaseParameterToThreePhase(psParams.line2_resistance), CPS::Math::singlePhaseParameterToThreePhase(psParams.line2_inductance),
                       CPS::Math::singlePhaseParameterToThreePhase(psParams.line2_capacitance));
  auto circuit_breaker = EMT::Ph3::Switch::make("circuit_breaker");
  circuit_breaker->setParameters(CPS::Math::singlePhaseParameterToThreePhase(psParams.circuit_breaker_open_resistance),
                                 CPS::Math::singlePhaseParameterToThreePhase(psParams.circuit_breaker_closed_resistance), true);
  auto load1 = EMT::Ph3::Resistor::make("load1");
  load1->setParameters(CPS::Math::singlePhaseParameterToThreePhase(psParams.load_resistance1));
  auto load1_switch = EMT::Ph3::Switch::make("load1_switch");
  load1_switch->setParameters(CPS::Math::singlePhaseParameterToThreePhase(psParams.circuit_breaker_open_resistance),
                                 CPS::Math::singlePhaseParameterToThreePhase(psParams.circuit_breaker_closed_resistance), true);
  auto load2 = EMT::Ph3::Resistor::make("load2");
  load2->setParameters(CPS::Math::singlePhaseParameterToThreePhase(psParams.load_resistance2));
  auto load2_switch = EMT::Ph3::Switch::make("load2_switch");
  load2_switch->setParameters(CPS::Math::singlePhaseParameterToThreePhase(psParams.circuit_breaker_open_resistance),
                                 CPS::Math::singlePhaseParameterToThreePhase(psParams.circuit_breaker_closed_resistance), false);

  // topology
  infeed_source->connect({EMT::SimNode::GND, node1});
  infeed_impedance->connect({node1, node4});
  converter1->connect({EMT::SimNode::GND, node2});
  converter2->connect({EMT::SimNode::GND, node3});
  line1->connect({node2, node4});
  line2->connect({node3, node4});
  circuit_breaker->connect({node4, node5});
  load1_switch->connect({node5, node6});
  load1->connect({node6, EMT::SimNode::GND});
  load2_switch->connect({node5, node7});
  load2->connect({node7, EMT::SimNode::GND});
  auto systemNodeList = SystemNodeList{node1, node2, node3, node4, node5, node6, node7};
  auto componentList = SystemComponentList{infeed_source, infeed_impedance, converter1, line1, converter2, line2, circuit_breaker, load1, load1_switch, load2, load2_switch};
  auto systemTopology = SystemTopology(psParams.frequency, systemNodeList, componentList);

  // events
  auto disconnect_load1 = DPsim::SwitchEvent3Ph::make(simParams.eventTime, load1_switch, false);
  auto connect_load2 = DPsim::SwitchEvent3Ph::make(simParams.eventTime, load2_switch, true);

  // logging
  auto logger = DataLogger::make(simName);
  logger->logAttribute("vInfeed", node1->attribute("v"));
  logger->logAttribute("vLoad", node5->attribute("v"));
  logger->logAttribute("iLoad", circuit_breaker->attribute("i_intf"));

  // simulation
  Simulation sim(simName, Logger::Level::info);
  sim.setSystem(systemTopology);
  sim.setTimeStep(simParams.timeStep);
  sim.setFinalTime(simParams.finalTime);
  sim.setDomain(Domain::EMT);
  sim.addLogger(logger);
  sim.addEvent(disconnect_load1);
  sim.addEvent(connect_load2);
  sim.run();
}

void simulate_DP(const SimulationParameters &simParams, const PowerSystemParameters &psParams) {
  String simName = "DP_simulation";
  Logger::setLogDir("logs/" + simName);

  // nodes
  auto node1 = DP::SimNode::make("node1", PhaseType::Single);
  auto node2 = DP::SimNode::make("node2", PhaseType::Single);
  auto node3 = DP::SimNode::make("node3", PhaseType::Single);
  auto node4 = DP::SimNode::make("node4", PhaseType::Single);
  auto node5 = DP::SimNode::make("node5", PhaseType::Single);
  auto node6 = DP::SimNode::make("node6", PhaseType::Single);
  auto node7 = DP::SimNode::make("node7", PhaseType::Single);

  // components
  auto infeed_source = DP::Ph1::VoltageSource::make("infeed_source");
  infeed_source->setParameters(CPS::Math::polar(psParams.voltage, 0.0));
  auto infeed_impedance = DP::Ph1::PiLine::make("infeed_impedance");
  infeed_impedance->setParameters(psParams.infeed_resistance, psParams.infeed_inductance, 0);
  auto converter1 = DP::Ph1::VoltageSource::make("converter1");
  converter1->setParameters(CPS::Math::polar(psParams.voltage, 0.0));
  auto line1 = DP::Ph1::PiLine::make("line1");
  line1->setParameters(psParams.line1_resistance, psParams.line1_inductance, psParams.line1_capacitance);
  auto converter2 = DP::Ph1::VoltageSource::make("converter2");
  converter2->setParameters(CPS::Math::polar(psParams.voltage, 0.0));
  auto line2 = DP::Ph1::PiLine::make("line2");
  line2->setParameters(psParams.line2_resistance, psParams.line2_inductance, psParams.line2_capacitance);
  auto circuit_breaker = DP::Ph1::Switch::make("circuit_breaker");
  circuit_breaker->setParameters(psParams.circuit_breaker_open_resistance, psParams.circuit_breaker_closed_resistance, true);
  auto load1 = DP::Ph1::Resistor::make("load");
  load1->setParameters(psParams.load_resistance1);
  auto load1_switch = DP::Ph1::Switch::make("load1_switch");
  load1_switch->setParameters(psParams.circuit_breaker_open_resistance, psParams.circuit_breaker_closed_resistance, true);
  auto load2 = DP::Ph1::Resistor::make("load2");
  load2->setParameters(psParams.load_resistance2);
  auto load2_switch = DP::Ph1::Switch::make("load2_switch");
  load2_switch->setParameters(psParams.circuit_breaker_open_resistance, psParams.circuit_breaker_closed_resistance, false);

  // topology
  infeed_source->connect({DP::SimNode::GND, node1});
  infeed_impedance->connect({node1, node4});
  converter1->connect({DP::SimNode::GND, node2});
  converter2->connect({DP::SimNode::GND, node3});
  line1->connect({node2, node4});
  line2->connect({node3, node4});
  circuit_breaker->connect({node4, node5});
  load1_switch->connect({node5, node6});
  load1->connect({node6, DP::SimNode::GND});
  load2_switch->connect({node5, node7});
  load2->connect({node7, DP::SimNode::GND});
  auto systemNodeList = SystemNodeList{node1, node2, node3, node4, node5, node6, node7};
  auto componentList = SystemComponentList{infeed_source, infeed_impedance, converter1, line1, converter2, line2, circuit_breaker, load1, load1_switch, load2, load2_switch};
  auto systemTopology = SystemTopology(psParams.frequency, systemNodeList, componentList);

  // events
  auto disconnect_load1 = DPsim::SwitchEvent::make(simParams.eventTime, load1_switch, false);
  auto connect_load2 = DPsim::SwitchEvent::make(simParams.eventTime, load2_switch, true);

  // logging
  auto logger = DataLogger::make(simName);
  logger->logAttribute("vInfeed", node1->attribute("v"));
  logger->logAttribute("vLoad", node5->attribute("v"));
  logger->logAttribute("iLoad", circuit_breaker->attribute("i_intf"));

  // simulation
  Simulation sim(simName, Logger::Level::info);
  sim.setSystem(systemTopology);
  sim.setTimeStep(simParams.timeStep);
  sim.setFinalTime(simParams.finalTime);
  sim.setDomain(Domain::DP);
  sim.addLogger(logger);
  sim.addEvent(disconnect_load1);
  sim.addEvent(connect_load2);
  sim.run();
}

void simulate_SP(const SimulationParameters &simParams, const PowerSystemParameters &psParams) {
  String simName = "SP_simulation";
  Logger::setLogDir("logs/" + simName);

  // nodes
  auto node1 = SP::SimNode::make("node1", PhaseType::Single);
  auto node2 = SP::SimNode::make("node2", PhaseType::Single);
  auto node3 = SP::SimNode::make("node3", PhaseType::Single);
  auto node4 = SP::SimNode::make("node4", PhaseType::Single);
  auto node5 = SP::SimNode::make("node5", PhaseType::Single);
  auto node6 = SP::SimNode::make("node6", PhaseType::Single);
  auto node7 = SP::SimNode::make("node7", PhaseType::Single);

  // components
  auto infeed_source = SP::Ph1::VoltageSource::make("infeed_source");
  infeed_source->setParameters(CPS::Math::polar(psParams.voltage, 0.0));
  auto infeed_impedance = SP::Ph1::PiLine::make("infeed_impedance");
  infeed_impedance->setParameters(psParams.infeed_resistance, psParams.infeed_inductance, 0);
  auto converter1 = SP::Ph1::VoltageSource::make("converter1");
  converter1->setParameters(CPS::Math::polar(psParams.voltage, 0.0));
  auto line1 = SP::Ph1::PiLine::make("line1");
  line1->setParameters(psParams.line1_resistance, psParams.line1_inductance, psParams.line1_capacitance);
  auto converter2 = SP::Ph1::VoltageSource::make("converter2");
  converter2->setParameters(CPS::Math::polar(psParams.voltage, 0.0));
  auto line2 = SP::Ph1::PiLine::make("line2");
  line2->setParameters(psParams.line2_resistance, psParams.line2_inductance, psParams.line2_capacitance);
  auto circuit_breaker = SP::Ph1::Switch::make("circuit_breaker");
  circuit_breaker->setParameters(psParams.circuit_breaker_open_resistance, psParams.circuit_breaker_closed_resistance, true);
  auto load1 = SP::Ph1::Resistor::make("load");
  load1->setParameters(psParams.load_resistance1);
  auto load1_switch = SP::Ph1::Switch::make("load1_switch");
  load1_switch->setParameters(psParams.circuit_breaker_open_resistance, psParams.circuit_breaker_closed_resistance, true);
  auto load2 = SP::Ph1::Resistor::make("load2");
  load2->setParameters(psParams.load_resistance2);
  auto load2_switch = SP::Ph1::Switch::make("load2_switch");
  load2_switch->setParameters(psParams.circuit_breaker_open_resistance, psParams.circuit_breaker_closed_resistance, false);

  // topology
  infeed_source->connect({SP::SimNode::GND, node1});
  infeed_impedance->connect({node1, node4});
  converter1->connect({SP::SimNode::GND, node2});
  converter2->connect({SP::SimNode::GND, node3});
  line1->connect({node2, node4});
  line2->connect({node3, node4});
  circuit_breaker->connect({node4, node5});
  load1_switch->connect({node5, node6});
  load1->connect({node6, SP::SimNode::GND});
  load2_switch->connect({node5, node7});
  load2->connect({node7, SP::SimNode::GND});
  auto systemNodeList = SystemNodeList{node1, node2, node3, node4, node5, node6, node7};
  auto componentList = SystemComponentList{infeed_source, infeed_impedance, converter1, line1, converter2, line2, circuit_breaker, load1, load1_switch, load2, load2_switch};
  auto systemTopology = SystemTopology(psParams.frequency, systemNodeList, componentList);

  // events
  auto disconnect_load1 = DPsim::SwitchEvent::make(simParams.eventTime, load1_switch, false);
  auto connect_load2 = DPsim::SwitchEvent::make(simParams.eventTime, load2_switch, true);

  // logging
  auto logger = DataLogger::make(simName);
  logger->logAttribute("vInfeed", node1->attribute("v"));
  logger->logAttribute("vLoad", node5->attribute("v"));
  logger->logAttribute("iLoad", circuit_breaker->attribute("i_intf"));

  // simulation
  Simulation sim(simName, Logger::Level::info);
  sim.setSystem(systemTopology);
  sim.setTimeStep(simParams.timeStep);
  sim.setFinalTime(simParams.finalTime);
  sim.setDomain(Domain::SP);
  sim.addEvent(disconnect_load1);
  sim.addEvent(connect_load2);
  sim.addLogger(logger);
  sim.run();
}

int main() {
  SimulationParameters simParams;
  PowerSystemParameters psParams;

  simulate_EMT(simParams, psParams);
  simulate_DP(simParams, psParams);
  simulate_SP(simParams, psParams);
  return 0;
}