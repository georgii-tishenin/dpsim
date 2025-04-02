#include <DPsim.h>
using namespace DPsim;
using namespace CPS;

struct SimulationParameters {
  double timeStep = 0.0001;
  double finalTime = 1.0;
};

struct PowerSystemParameters {
  double frequency = 50;
  double voltage = 10e3;
  double infeed_resistance = 1;
  double infeed_inductance = 0.01;
  double line1_resistance = 1;
  double line1_inductance = 0.01;
  double line2_resistance = 1;
  double line2_inductance = 0.01;
  double circuit_breaker_closed_resistance = 1e-4;
  double circuit_breaker_open_resistance = 1e6;
  double load_resistance = 100;
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
                       CPS::Math::singlePhaseParameterToThreePhase(0));
  auto converter2 = EMT::Ph3::VoltageSource::make("converter2");
  converter2->setParameters(CPS::Math::singlePhaseVariableToThreePhase(CPS::Math::polar(psParams.voltage, 0.0)), psParams.frequency);
  auto line2 = EMT::Ph3::PiLine::make("line2");
  line2->setParameters(CPS::Math::singlePhaseParameterToThreePhase(psParams.line2_resistance), CPS::Math::singlePhaseParameterToThreePhase(psParams.line2_inductance),
                       CPS::Math::singlePhaseParameterToThreePhase(0));
  auto circuit_breaker = EMT::Ph3::Switch::make("circuit_breaker");
  circuit_breaker->setParameters(CPS::Math::singlePhaseParameterToThreePhase(psParams.circuit_breaker_open_resistance),
                                 CPS::Math::singlePhaseParameterToThreePhase(psParams.circuit_breaker_closed_resistance), true);
  auto load = EMT::Ph3::Resistor::make("load");
  load->setParameters(CPS::Math::singlePhaseParameterToThreePhase(psParams.load_resistance));

  // topology
  infeed_source->connect({EMT::SimNode::GND, node1});
  infeed_impedance->connect({node1, node4});
  converter1->connect({EMT::SimNode::GND, node2});
  converter2->connect({EMT::SimNode::GND, node3});
  line1->connect({node2, node4});
  line2->connect({node3, node4});
  circuit_breaker->connect({node4, node5});
  load->connect({node5, EMT::SimNode::GND});
  auto systemNodeList = SystemNodeList{node1, node2, node3, node4, node5};
  auto componentList = SystemComponentList{infeed_source, infeed_impedance, converter1, line1, converter2, line2, circuit_breaker, load};
  auto systemTopology = SystemTopology(psParams.frequency, systemNodeList, componentList);

  // logging
  auto logger = DataLogger::make(simName);
  logger->logAttribute("vInfeed", node1->attribute("v"));
  logger->logAttribute("vLoad", node5->attribute("v"));
  logger->logAttribute("iLoad", load->attribute("i_intf"));

  // simulation
  Simulation sim(simName, Logger::Level::info);
  sim.setSystem(systemTopology);
  sim.setTimeStep(simParams.timeStep);
  sim.setFinalTime(simParams.finalTime);
  sim.setDomain(Domain::EMT);
  sim.addLogger(logger);
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

  // components
  auto infeed_source = DP::Ph1::VoltageSource::make("infeed_source");
  infeed_source->setParameters(CPS::Math::polar(psParams.voltage, 0.0));
  auto infeed_impedance = DP::Ph1::PiLine::make("infeed_impedance");
  infeed_impedance->setParameters(psParams.infeed_resistance, psParams.infeed_inductance, 0);
  auto converter1 = DP::Ph1::VoltageSource::make("converter1");
  converter1->setParameters(CPS::Math::polar(psParams.voltage, 0.0));
  auto line1 = DP::Ph1::PiLine::make("line1");
  line1->setParameters(psParams.line1_resistance, psParams.line1_inductance, 0);
  auto converter2 = DP::Ph1::VoltageSource::make("converter2");
  converter2->setParameters(CPS::Math::polar(psParams.voltage, 0.0));
  auto line2 = DP::Ph1::PiLine::make("line2");
  line2->setParameters(psParams.line2_resistance, psParams.line2_inductance, 0);
  auto circuit_breaker = DP::Ph1::Switch::make("circuit_breaker");
  circuit_breaker->setParameters(psParams.circuit_breaker_open_resistance, psParams.circuit_breaker_closed_resistance, true);
  auto load = DP::Ph1::Resistor::make("load");
  load->setParameters(psParams.load_resistance);

  // topology
  infeed_source->connect({DP::SimNode::GND, node1});
  infeed_impedance->connect({node1, node4});
  converter1->connect({DP::SimNode::GND, node2});
  converter2->connect({DP::SimNode::GND, node3});
  line1->connect({node2, node4});
  line2->connect({node3, node4});
  circuit_breaker->connect({node4, node5});
  load->connect({node5, DP::SimNode::GND});
  auto systemNodeList = SystemNodeList{node1, node2, node3, node4, node5};
  auto componentList = SystemComponentList{infeed_source, infeed_impedance, converter1, line1, converter2, line2, circuit_breaker, load};
  auto systemTopology = SystemTopology(psParams.frequency, systemNodeList, componentList);

  // logging
  auto logger = DataLogger::make(simName);
  logger->logAttribute("vInfeed", node1->attribute("v"));
  logger->logAttribute("vLoad", node5->attribute("v"));
  logger->logAttribute("iLoad", load->attribute("i_intf"));

  // simulation
  Simulation sim(simName, Logger::Level::info);
  sim.setSystem(systemTopology);
  sim.setTimeStep(simParams.timeStep);
  sim.setFinalTime(simParams.finalTime);
  sim.setDomain(Domain::DP);
  sim.addLogger(logger);
  sim.run();
}

int main() {
  SimulationParameters simParams;
  PowerSystemParameters psParams;

  simulate_EMT(simParams, psParams);
  simulate_DP(simParams, psParams);
  return 0;
}