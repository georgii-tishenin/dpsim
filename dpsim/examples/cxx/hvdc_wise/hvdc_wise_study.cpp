#include <DPsim.h>
using namespace DPsim;
using namespace CPS;

void simulate_EMT() {
  String simName = "EMT_simulation";
  Logger::setLogDir("logs/" + simName);

  // simulation parameters
  double timeStep = 0.0001;
  double finalTime = 0.1;

  // power system parameters
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

  // nodes
  auto node1 = EMT::SimNode::make("node1", PhaseType::ABC);
  auto node2 = EMT::SimNode::make("node2", PhaseType::ABC);
  auto node3 = EMT::SimNode::make("node3", PhaseType::ABC);
  auto node4 = EMT::SimNode::make("node4", PhaseType::ABC);
  auto node5 = EMT::SimNode::make("node5", PhaseType::ABC);

  // components
  auto infeed_source = EMT::Ph3::VoltageSource::make("infeed_source");
  infeed_source->setParameters(CPS::Math::singlePhaseVariableToThreePhase(CPS::Math::polar(voltage, 0.0)), frequency);
  auto infeed_impedance = EMT::Ph3::PiLine::make("infeed_impedance");
  infeed_impedance->setParameters(CPS::Math::singlePhaseParameterToThreePhase(infeed_resistance), CPS::Math::singlePhaseParameterToThreePhase(infeed_inductance),
                                  CPS::Math::singlePhaseParameterToThreePhase(0));
  auto converter1 = EMT::Ph3::VoltageSource::make("converter1");
  converter1->setParameters(CPS::Math::singlePhaseVariableToThreePhase(CPS::Math::polar(voltage, 0.0)), frequency);
  auto line1 = EMT::Ph3::PiLine::make("line1");
  line1->setParameters(CPS::Math::singlePhaseParameterToThreePhase(line1_resistance), CPS::Math::singlePhaseParameterToThreePhase(line1_inductance), CPS::Math::singlePhaseParameterToThreePhase(0));
  auto converter2 = EMT::Ph3::VoltageSource::make("converter2");
  converter2->setParameters(CPS::Math::singlePhaseVariableToThreePhase(CPS::Math::polar(voltage, 0.0)), frequency);
  auto line2 = EMT::Ph3::PiLine::make("line2");
  line2->setParameters(CPS::Math::singlePhaseParameterToThreePhase(line2_resistance), CPS::Math::singlePhaseParameterToThreePhase(line2_inductance), CPS::Math::singlePhaseParameterToThreePhase(0));
  auto circuit_breaker = EMT::Ph3::Switch::make("circuit_breaker");
  circuit_breaker->setParameters(CPS::Math::singlePhaseParameterToThreePhase(circuit_breaker_open_resistance), CPS::Math::singlePhaseParameterToThreePhase(circuit_breaker_closed_resistance), true);
  auto load = EMT::Ph3::Resistor::make("load");
  load->setParameters(CPS::Math::singlePhaseParameterToThreePhase(load_resistance));

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
  auto systemTopology = SystemTopology(frequency, systemNodeList, componentList);

  // logging
  auto logger = DataLogger::make(simName);
  logger->logAttribute("infeedVoltage", node1->attribute("v"));
  logger->logAttribute("loadVoltage", node5->attribute("v"));
  logger->logAttribute("loadCurrent", load->attribute("i_intf"));

  // simulation
  Simulation sim(simName, Logger::Level::info);
  sim.setSystem(systemTopology);
  sim.setTimeStep(timeStep);
  sim.setFinalTime(finalTime);
  sim.setDomain(Domain::EMT);
  sim.addLogger(logger);
  sim.run();
}

int main() {
  simulate_EMT();
  return 0;
}
