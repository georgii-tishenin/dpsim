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
  double infeedResistance = 1;
  double infeedInductance = 0.5;
  double line1Resistance = 1;
  double line1Inductance = 0.1;
  double line1Capacitance = 3e-6;
  double line2Resistance = 2;
  double line2Inductance = 0.2;
  double line2Capacitance = 2e-6;
  double switchClosedResistance = 1e-4;
  double switchOpenResistance = 1e6;
  double loadResistance1 = 100;
  double loadResistance2 = 50;
};

void simulateEMT(const SimulationParameters &simParams, const PowerSystemParameters &psParams) {
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
  auto infeedSource = EMT::Ph3::VoltageSource::make("infeed_source");
  infeedSource->setParameters(CPS::Math::singlePhaseVariableToThreePhase(CPS::Math::polar(psParams.voltage, 0.0)), psParams.frequency);
  infeedSource->connect({EMT::SimNode::GND, node1});

  auto infeedImpedance = EMT::Ph3::PiLine::make("infeed_impedance");
  infeedImpedance->setParameters(CPS::Math::singlePhaseParameterToThreePhase(psParams.infeedResistance), CPS::Math::singlePhaseParameterToThreePhase(psParams.infeedInductance),
                                 CPS::Math::singlePhaseParameterToThreePhase(0));
  infeedImpedance->connect({node1, node4});

  auto converter1 = EMT::Ph3::VoltageSource::make("converter1");
  converter1->setParameters(CPS::Math::singlePhaseVariableToThreePhase(CPS::Math::polar(psParams.voltage, 0.0)), psParams.frequency);
  converter1->connect({EMT::SimNode::GND, node2});

  auto line1 = EMT::Ph3::PiLine::make("line1");
  line1->setParameters(CPS::Math::singlePhaseParameterToThreePhase(psParams.line1Resistance), CPS::Math::singlePhaseParameterToThreePhase(psParams.line1Inductance),
                       CPS::Math::singlePhaseParameterToThreePhase(psParams.line1Capacitance));
  line1->connect({node2, node4});

  auto converter2 = EMT::Ph3::VoltageSource::make("converter2");
  converter2->setParameters(CPS::Math::singlePhaseVariableToThreePhase(CPS::Math::polar(psParams.voltage, 0.0)), psParams.frequency);
  converter2->connect({EMT::SimNode::GND, node3});

  auto line2 = EMT::Ph3::PiLine::make("line2");
  line2->setParameters(CPS::Math::singlePhaseParameterToThreePhase(psParams.line2Resistance), CPS::Math::singlePhaseParameterToThreePhase(psParams.line2Inductance),
                       CPS::Math::singlePhaseParameterToThreePhase(psParams.line2Capacitance));
  line2->connect({node3, node4});

  auto circuitBreaker = EMT::Ph3::Switch::make("circuit_breaker");
  circuitBreaker->setParameters(CPS::Math::singlePhaseParameterToThreePhase(psParams.switchOpenResistance), CPS::Math::singlePhaseParameterToThreePhase(psParams.switchClosedResistance), true);
  circuitBreaker->connect({node4, node5});

  auto load1 = EMT::Ph3::Resistor::make("load1");
  load1->setParameters(CPS::Math::singlePhaseParameterToThreePhase(psParams.loadResistance1));
  load1->connect({node6, EMT::SimNode::GND});

  auto load1Switch = EMT::Ph3::Switch::make("load1_switch");
  load1Switch->setParameters(CPS::Math::singlePhaseParameterToThreePhase(psParams.switchOpenResistance), CPS::Math::singlePhaseParameterToThreePhase(psParams.switchClosedResistance), true);
  load1Switch->connect({node5, node6});

  auto load2 = EMT::Ph3::Resistor::make("load2");
  load2->setParameters(CPS::Math::singlePhaseParameterToThreePhase(psParams.loadResistance2));
  load2->connect({node7, EMT::SimNode::GND});

  auto load2Switch = EMT::Ph3::Switch::make("load2_switch");
  load2Switch->setParameters(CPS::Math::singlePhaseParameterToThreePhase(psParams.switchOpenResistance), CPS::Math::singlePhaseParameterToThreePhase(psParams.switchClosedResistance), false);
  load2Switch->connect({node5, node7});

  // topology
  auto systemNodeList = SystemNodeList{node1, node2, node3, node4, node5, node6, node7};
  auto componentList = SystemComponentList{infeedSource, infeedImpedance, converter1, line1, converter2, line2, circuitBreaker, load1, load1Switch, load2, load2Switch};
  auto systemTopology = SystemTopology(psParams.frequency, systemNodeList, componentList);

  // events
  auto disconnectLoad1 = DPsim::SwitchEvent3Ph::make(simParams.eventTime, load1Switch, false);
  auto connectLoad2 = DPsim::SwitchEvent3Ph::make(simParams.eventTime, load2Switch, true);

  // logging
  auto logger = DataLogger::make(simName);
  logger->logAttribute("vInfeed", node1->attribute("v"));
  logger->logAttribute("vLoad", node5->attribute("v"));
  logger->logAttribute("iLoad", circuitBreaker->attribute("i_intf"));

  // simulation
  Simulation sim(simName, Logger::Level::info);
  sim.setSystem(systemTopology);
  sim.setTimeStep(simParams.timeStep);
  sim.setFinalTime(simParams.finalTime);
  sim.setDomain(Domain::EMT);
  sim.addLogger(logger);
  sim.addEvent(disconnectLoad1);
  sim.addEvent(connectLoad2);
  sim.run();
}

void simulateDP(const SimulationParameters &simParams, const PowerSystemParameters &psParams) {
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
  auto infeedSource = DP::Ph1::VoltageSource::make("infeed_source");
  infeedSource->setParameters(CPS::Math::polar(psParams.voltage, 0.0));
  infeedSource->connect({DP::SimNode::GND, node1});

  auto infeedImpedance = DP::Ph1::PiLine::make("infeed_impedance");
  infeedImpedance->setParameters(psParams.infeedResistance, psParams.infeedInductance, 0);
  infeedImpedance->connect({node1, node4});

  auto converter1 = DP::Ph1::VoltageSource::make("converter1");
  converter1->setParameters(CPS::Math::polar(psParams.voltage, 0.0));
  converter1->connect({DP::SimNode::GND, node2});

  auto line1 = DP::Ph1::PiLine::make("line1");
  line1->setParameters(psParams.line1Resistance, psParams.line1Inductance, psParams.line1Capacitance);
  line1->connect({node2, node4});

  auto converter2 = DP::Ph1::VoltageSource::make("converter2");
  converter2->setParameters(CPS::Math::polar(psParams.voltage, 0.0));
  converter2->connect({DP::SimNode::GND, node3});

  auto line2 = DP::Ph1::PiLine::make("line2");
  line2->setParameters(psParams.line2Resistance, psParams.line2Inductance, psParams.line2Capacitance);
  line2->connect({node3, node4});

  auto circuitBreaker = DP::Ph1::Switch::make("circuit_breaker");
  circuitBreaker->setParameters(psParams.switchOpenResistance, psParams.switchClosedResistance, true);
  circuitBreaker->connect({node4, node5});

  auto load1 = DP::Ph1::Resistor::make("load");
  load1->setParameters(psParams.loadResistance1);
  load1->connect({node6, DP::SimNode::GND});

  auto load1Switch = DP::Ph1::Switch::make("load1_switch");
  load1Switch->setParameters(psParams.switchOpenResistance, psParams.switchClosedResistance, true);
  load1Switch->connect({node5, node6});

  auto load2 = DP::Ph1::Resistor::make("load2");
  load2->setParameters(psParams.loadResistance2);
  load2->connect({node7, DP::SimNode::GND});

  auto load2Switch = DP::Ph1::Switch::make("load2_switch");
  load2Switch->setParameters(psParams.switchOpenResistance, psParams.switchClosedResistance, false);
  load2Switch->connect({node5, node7});

  // topology
  auto systemNodeList = SystemNodeList{node1, node2, node3, node4, node5, node6, node7};
  auto componentList = SystemComponentList{infeedSource, infeedImpedance, converter1, line1, converter2, line2, circuitBreaker, load1, load1Switch, load2, load2Switch};
  auto systemTopology = SystemTopology(psParams.frequency, systemNodeList, componentList);

  // events
  auto disconnectLoad1 = DPsim::SwitchEvent::make(simParams.eventTime, load1Switch, false);
  auto connectLoad2 = DPsim::SwitchEvent::make(simParams.eventTime, load2Switch, true);

  // logging
  auto logger = DataLogger::make(simName);
  logger->logAttribute("vInfeed", node1->attribute("v"));
  logger->logAttribute("vLoad", node5->attribute("v"));
  logger->logAttribute("iLoad", circuitBreaker->attribute("i_intf"));

  // simulation
  Simulation sim(simName, Logger::Level::info);
  sim.setSystem(systemTopology);
  sim.setTimeStep(simParams.timeStep);
  sim.setFinalTime(simParams.finalTime);
  sim.setDomain(Domain::DP);
  sim.addLogger(logger);
  sim.addEvent(disconnectLoad1);
  sim.addEvent(connectLoad2);
  sim.run();
}

void simulateSP(const SimulationParameters &simParams, const PowerSystemParameters &psParams) {
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
  auto infeedSource = SP::Ph1::VoltageSource::make("infeed_source");
  infeedSource->setParameters(CPS::Math::polar(psParams.voltage, 0.0));
  infeedSource->connect({SP::SimNode::GND, node1});

  auto infeedImpedance = SP::Ph1::PiLine::make("infeed_impedance");
  infeedImpedance->setParameters(psParams.infeedResistance, psParams.infeedInductance, 0);
  infeedImpedance->connect({node1, node4});

  auto converter1 = SP::Ph1::VoltageSource::make("converter1");
  converter1->setParameters(CPS::Math::polar(psParams.voltage, 0.0));
  converter1->connect({SP::SimNode::GND, node2});

  auto line1 = SP::Ph1::PiLine::make("line1");
  line1->setParameters(psParams.line1Resistance, psParams.line1Inductance, psParams.line1Capacitance);
  line1->connect({node2, node4});

  auto converter2 = SP::Ph1::VoltageSource::make("converter2");
  converter2->setParameters(CPS::Math::polar(psParams.voltage, 0.0));
  converter2->connect({SP::SimNode::GND, node3});

  auto line2 = SP::Ph1::PiLine::make("line2");
  line2->setParameters(psParams.line2Resistance, psParams.line2Inductance, psParams.line2Capacitance);
  line2->connect({node3, node4});

  auto circuitBreaker = SP::Ph1::Switch::make("circuit_breaker");
  circuitBreaker->setParameters(psParams.switchOpenResistance, psParams.switchClosedResistance, true);
  circuitBreaker->connect({node4, node5});

  auto load1 = SP::Ph1::Resistor::make("load");
  load1->setParameters(psParams.loadResistance1);
  load1->connect({node6, SP::SimNode::GND});

  auto load1Switch = SP::Ph1::Switch::make("load1_switch");
  load1Switch->setParameters(psParams.switchOpenResistance, psParams.switchClosedResistance, true);
  load1Switch->connect({node5, node6});

  auto load2 = SP::Ph1::Resistor::make("load2");
  load2->setParameters(psParams.loadResistance2);
  load2->connect({node7, SP::SimNode::GND});

  auto load2Switch = SP::Ph1::Switch::make("load2_switch");
  load2Switch->setParameters(psParams.switchOpenResistance, psParams.switchClosedResistance, false);
  load2Switch->connect({node5, node7});

  // topology
  auto systemNodeList = SystemNodeList{node1, node2, node3, node4, node5, node6, node7};
  auto componentList = SystemComponentList{infeedSource, infeedImpedance, converter1, line1, converter2, line2, circuitBreaker, load1, load1Switch, load2, load2Switch};
  auto systemTopology = SystemTopology(psParams.frequency, systemNodeList, componentList);

  // events
  auto disconnectLoad1 = DPsim::SwitchEvent::make(simParams.eventTime, load1Switch, false);
  auto connectLoad2 = DPsim::SwitchEvent::make(simParams.eventTime, load2Switch, true);

  // logging
  auto logger = DataLogger::make(simName);
  logger->logAttribute("vInfeed", node1->attribute("v"));
  logger->logAttribute("vLoad", node5->attribute("v"));
  logger->logAttribute("iLoad", circuitBreaker->attribute("i_intf"));

  // simulation
  Simulation sim(simName, Logger::Level::info);
  sim.setSystem(systemTopology);
  sim.setTimeStep(simParams.timeStep);
  sim.setFinalTime(simParams.finalTime);
  sim.setDomain(Domain::SP);
  sim.addEvent(disconnectLoad1);
  sim.addEvent(connectLoad2);
  sim.addLogger(logger);
  sim.run();
}

int main() {
  SimulationParameters simParams;
  PowerSystemParameters psParams;

  simulateEMT(simParams, psParams);
  simulateDP(simParams, psParams);
  simulateSP(simParams, psParams);
  return 0;
}