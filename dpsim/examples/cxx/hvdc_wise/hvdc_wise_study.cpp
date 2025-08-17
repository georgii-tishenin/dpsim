#include <DPsim.h>
using namespace DPsim;
using namespace CPS;

namespace VariableNames {
constexpr const char *vInfeed = "vInfeed";
constexpr const char *iInfeed = "iInfeed";
constexpr const char *vConverter1 = "vConverter1";
constexpr const char *iConverter1 = "iConverter1";
constexpr const char *vConverter2 = "vConverter2";
constexpr const char *iConverter2 = "iConverter2";
constexpr const char *vLoad = "vLoad";
constexpr const char *iLoad = "iLoad";
} // namespace VariableNames

namespace SwitchConstants {
constexpr double closedResistance = 1e-4;
constexpr double openResistance = 1e6;
} // namespace SwitchConstants

namespace AttributeNames {
constexpr const char *v = "v";
constexpr const char *i = "i_intf";
} // namespace AttributeNames

struct SimulationParameters {
  double timeStep = 1e-5;
  double eventTime = 0.1;
  double finalTime = 0.2;
};

struct PowerSystemInputParameters {
  double frequency = 50;
  double baseVoltageLineToLine = 110e3;
  double baseThreePhasePower = 100e6;

  // infeed parameters
  double infeedResistanceInPerUnit = 0.01;
  double infeedReactanceInPerUnit = 0.1;

  // line1 parameters
  double line1LengthInKm = 10;
  double line1ResistancePerKm = 0.1;
  double line1ReactancePerKm = 0.4;
  double line1CapacitancePerKm = 0;

  // line2 parameters
  double line2LengthInKm = 10;
  double line2ResistancePerKm = 0.2;
  double line2ReactancePerKm = 0.4;
  double lin2CapacitancePerKm = 0;

  // load parameters
  double load1InPerUnit = 0.1;
  double load2InPerUnit = 0.2;
};

struct PowerSystemParameters {
  double frequency;
  double voltageLineToLine;
  double voltageLineToGround;
  double infeedResistance;
  double infeedInductance;
  double line1Resistance;
  double line1Inductance;
  double line1Capacitance;
  double line2Resistance;
  double line2Inductance;
  double line2Capacitance;
  double loadResistance1;
  double loadResistance2;

  PowerSystemParameters(double freq, double voltLineToLine,
                        double voltLineToGround, double infeedRes,
                        double infeedInd, double line1Res, double line1Ind,
                        double line1Cap, double line2Res, double line2Ind,
                        double line2Cap, double loadRes1, double loadRes2)
      : frequency(freq), voltageLineToLine(voltLineToLine),
        voltageLineToGround(voltLineToGround), infeedResistance(infeedRes),
        infeedInductance(infeedInd), line1Resistance(line1Res),
        line1Inductance(line1Ind), line1Capacitance(line1Cap),
        line2Resistance(line2Res), line2Inductance(line2Ind),
        line2Capacitance(line2Cap), loadResistance1(loadRes1),
        loadResistance2(loadRes2) {}
};

PowerSystemParameters
calculatePowerSystemParameters(const PowerSystemInputParameters &inputParams) {
  double voltageLineToGround = inputParams.baseVoltageLineToLine / sqrt(3);
  double baseImpedance = inputParams.baseVoltageLineToLine *
                         inputParams.baseVoltageLineToLine /
                         inputParams.baseThreePhasePower;
  double omega = 2 * M_PI * inputParams.frequency;

  double infeedResistance =
      inputParams.infeedResistanceInPerUnit * baseImpedance;
  double infeedInductance =
      inputParams.infeedReactanceInPerUnit * baseImpedance / omega;

  double line1Resistance =
      inputParams.line1ResistancePerKm * inputParams.line1LengthInKm;
  double line1Inductance =
      inputParams.line1ReactancePerKm * inputParams.line1LengthInKm / omega;
  double line1Capacitance =
      inputParams.line1CapacitancePerKm * inputParams.line1LengthInKm;

  double line2Resistance =
      inputParams.line2ResistancePerKm * inputParams.line2LengthInKm;
  double line2Inductance =
      inputParams.line2ReactancePerKm * inputParams.line2LengthInKm / omega;
  double line2Capacitance =
      inputParams.lin2CapacitancePerKm * inputParams.line2LengthInKm;

  double loadResistance1 = baseImpedance / inputParams.load1InPerUnit;
  double loadResistance2 = baseImpedance / inputParams.load2InPerUnit;

  return PowerSystemParameters(
      inputParams.frequency, inputParams.baseVoltageLineToLine,
      voltageLineToGround, infeedResistance, infeedInductance, line1Resistance,
      line1Inductance, line1Capacitance, line2Resistance, line2Inductance,
      line2Capacitance, loadResistance1, loadResistance2);
}

Simulation setupSimulation(const std::string &simName,
                           const SimulationParameters &simParams,
                           const SystemTopology &systemTopology,
                           const std::shared_ptr<DataLogger> &logger,
                           Domain domain) {
  Simulation sim(simName, Logger::Level::info);
  sim.setSystem(systemTopology);
  sim.setTimeStep(simParams.timeStep);
  sim.setFinalTime(simParams.finalTime);
  sim.setDomain(domain);
  sim.addLogger(logger);
  return sim;
}

void simulateEMT(const SimulationParameters &simParams,
                 const PowerSystemParameters &psParams) {
  String simName = "EMT_simulation";
  Logger::setLogDir("logs/" + simName);
  auto logger = DataLogger::make(simName);

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
  infeedSource->setParameters(
      CPS::Math::singlePhaseVariableToThreePhase(
          CPS::Math::polar(psParams.voltageLineToLine, 0.0)),
      psParams.frequency);
  infeedSource->connect({EMT::SimNode::GND, node1});
  auto infeedImpedance = EMT::Ph3::PiLine::make("infeed_impedance");
  infeedImpedance->setParameters(
      CPS::Math::singlePhaseParameterToThreePhase(psParams.infeedResistance),
      CPS::Math::singlePhaseParameterToThreePhase(psParams.infeedInductance),
      CPS::Math::singlePhaseParameterToThreePhase(0));
  infeedImpedance->connect({node1, node4});
  logger->logAttribute(VariableNames::vInfeed,
                       node4->attribute(AttributeNames::v));
  logger->logAttribute(VariableNames::iInfeed,
                       infeedImpedance->attribute(AttributeNames::i));

  auto converter1 = EMT::Ph3::VoltageSource::make("converter1");
  converter1->setParameters(
      CPS::Math::singlePhaseVariableToThreePhase(
          CPS::Math::polar(psParams.voltageLineToLine, 0.0)),
      psParams.frequency);
  converter1->connect({EMT::SimNode::GND, node2});
  logger->logAttribute(VariableNames::vConverter1,
                       node2->attribute(AttributeNames::v));
  logger->logAttribute(VariableNames::iConverter1,
                       converter1->attribute(AttributeNames::i));

  auto line1 = EMT::Ph3::PiLine::make("line1");
  line1->setParameters(
      CPS::Math::singlePhaseParameterToThreePhase(psParams.line1Resistance),
      CPS::Math::singlePhaseParameterToThreePhase(psParams.line1Inductance),
      CPS::Math::singlePhaseParameterToThreePhase(psParams.line1Capacitance));
  line1->connect({node2, node4});

  auto converter2 = EMT::Ph3::VoltageSource::make("converter2");
  converter2->setParameters(
      CPS::Math::singlePhaseVariableToThreePhase(
          CPS::Math::polar(psParams.voltageLineToLine, 0.0)),
      psParams.frequency);
  converter2->connect({EMT::SimNode::GND, node3});
  logger->logAttribute(VariableNames::vConverter2,
                       node3->attribute(AttributeNames::v));
  logger->logAttribute(VariableNames::iConverter2,
                       converter2->attribute(AttributeNames::i));

  auto line2 = EMT::Ph3::PiLine::make("line2");
  line2->setParameters(
      CPS::Math::singlePhaseParameterToThreePhase(psParams.line2Resistance),
      CPS::Math::singlePhaseParameterToThreePhase(psParams.line2Inductance),
      CPS::Math::singlePhaseParameterToThreePhase(psParams.line2Capacitance));
  line2->connect({node3, node4});

  auto circuitBreaker = EMT::Ph3::Switch::make("circuit_breaker");
  circuitBreaker->setParameters(CPS::Math::singlePhaseParameterToThreePhase(
                                    SwitchConstants::openResistance),
                                CPS::Math::singlePhaseParameterToThreePhase(
                                    SwitchConstants::closedResistance),
                                true);
  circuitBreaker->connect({node4, node5});
  logger->logAttribute(VariableNames::vLoad,
                       node5->attribute(AttributeNames::v));
  logger->logAttribute(VariableNames::iLoad,
                       circuitBreaker->attribute(AttributeNames::i));

  auto load1 = EMT::Ph3::Resistor::make("load1");
  load1->setParameters(
      CPS::Math::singlePhaseParameterToThreePhase(psParams.loadResistance1));
  load1->connect({node6, EMT::SimNode::GND});

  auto load1Switch = EMT::Ph3::Switch::make("load1_switch");
  load1Switch->setParameters(CPS::Math::singlePhaseParameterToThreePhase(
                                 SwitchConstants::openResistance),
                             CPS::Math::singlePhaseParameterToThreePhase(
                                 SwitchConstants::closedResistance),
                             true);
  load1Switch->connect({node5, node6});

  auto load2 = EMT::Ph3::Resistor::make("load2");
  load2->setParameters(
      CPS::Math::singlePhaseParameterToThreePhase(psParams.loadResistance2));
  load2->connect({node7, EMT::SimNode::GND});

  auto load2Switch = EMT::Ph3::Switch::make("load2_switch");
  load2Switch->setParameters(CPS::Math::singlePhaseParameterToThreePhase(
                                 SwitchConstants::openResistance),
                             CPS::Math::singlePhaseParameterToThreePhase(
                                 SwitchConstants::closedResistance),
                             false);
  load2Switch->connect({node5, node7});

  // topology
  auto systemNodeList =
      SystemNodeList{node1, node2, node3, node4, node5, node6, node7};
  auto componentList = SystemComponentList{
      infeedSource,   infeedImpedance, converter1,  line1, converter2, line2,
      circuitBreaker, load1,           load1Switch, load2, load2Switch};
  auto systemTopology =
      SystemTopology(psParams.frequency, systemNodeList, componentList);

  // events
  auto disconnectLoad1 =
      DPsim::SwitchEvent3Ph::make(simParams.eventTime, load1Switch, false);
  auto connectLoad2 =
      DPsim::SwitchEvent3Ph::make(simParams.eventTime, load2Switch, true);

  // simulation
  auto sim =
      setupSimulation(simName, simParams, systemTopology, logger, Domain::EMT);
  sim.addEvent(disconnectLoad1);
  sim.addEvent(connectLoad2);
  sim.run();
}

void simulateDP(const SimulationParameters &simParams,
                const PowerSystemParameters &psParams) {
  String simName = "DP_simulation";
  Logger::setLogDir("logs/" + simName);
  auto logger = DataLogger::make(simName);

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
  infeedSource->setParameters(
      CPS::Math::polar(psParams.voltageLineToGround, 0.0));
  infeedSource->connect({DP::SimNode::GND, node1});
  auto infeedImpedance = DP::Ph1::PiLine::make("infeed_impedance");
  infeedImpedance->setParameters(psParams.infeedResistance,
                                 psParams.infeedInductance, 0);
  infeedImpedance->connect({node1, node4});
  logger->logAttribute(VariableNames::vInfeed,
                       node4->attribute(AttributeNames::v));
  logger->logAttribute(VariableNames::iInfeed,
                       infeedImpedance->attribute(AttributeNames::i));

  auto converter1 = DP::Ph1::VoltageSource::make("converter1");
  converter1->setParameters(
      CPS::Math::polar(psParams.voltageLineToGround, 0.0));
  converter1->connect({DP::SimNode::GND, node2});
  logger->logAttribute(VariableNames::vConverter1,
                       node2->attribute(AttributeNames::v));
  logger->logAttribute(VariableNames::iConverter1,
                       converter1->attribute(AttributeNames::i));

  auto line1 = DP::Ph1::PiLine::make("line1");
  line1->setParameters(psParams.line1Resistance, psParams.line1Inductance,
                       psParams.line1Capacitance);
  line1->connect({node2, node4});

  auto converter2 = DP::Ph1::VoltageSource::make("converter2");
  converter2->setParameters(
      CPS::Math::polar(psParams.voltageLineToGround, 0.0));
  converter2->connect({DP::SimNode::GND, node3});
  logger->logAttribute(VariableNames::vConverter2,
                       node3->attribute(AttributeNames::v));
  logger->logAttribute(VariableNames::iConverter2,
                       converter2->attribute(AttributeNames::i));

  auto line2 = DP::Ph1::PiLine::make("line2");
  line2->setParameters(psParams.line2Resistance, psParams.line2Inductance,
                       psParams.line2Capacitance);
  line2->connect({node3, node4});

  auto circuitBreaker = DP::Ph1::Switch::make("circuit_breaker");
  circuitBreaker->setParameters(SwitchConstants::openResistance,
                                SwitchConstants::closedResistance, true);
  circuitBreaker->connect({node4, node5});
  logger->logAttribute(VariableNames::vLoad,
                       node5->attribute(AttributeNames::v));
  logger->logAttribute(VariableNames::iLoad,
                       circuitBreaker->attribute(AttributeNames::i));

  auto load1 = DP::Ph1::Resistor::make("load");
  load1->setParameters(psParams.loadResistance1);
  load1->connect({node6, DP::SimNode::GND});

  auto load1Switch = DP::Ph1::Switch::make("load1_switch");
  load1Switch->setParameters(SwitchConstants::openResistance,
                             SwitchConstants::closedResistance, true);
  load1Switch->connect({node5, node6});

  auto load2 = DP::Ph1::Resistor::make("load2");
  load2->setParameters(psParams.loadResistance2);
  load2->connect({node7, DP::SimNode::GND});

  auto load2Switch = DP::Ph1::Switch::make("load2_switch");
  load2Switch->setParameters(SwitchConstants::openResistance,
                             SwitchConstants::closedResistance, false);
  load2Switch->connect({node5, node7});

  // topology
  auto systemNodeList =
      SystemNodeList{node1, node2, node3, node4, node5, node6, node7};
  auto componentList = SystemComponentList{
      infeedSource,   infeedImpedance, converter1,  line1, converter2, line2,
      circuitBreaker, load1,           load1Switch, load2, load2Switch};
  auto systemTopology =
      SystemTopology(psParams.frequency, systemNodeList, componentList);

  // events
  auto disconnectLoad1 =
      DPsim::SwitchEvent::make(simParams.eventTime, load1Switch, false);
  auto connectLoad2 =
      DPsim::SwitchEvent::make(simParams.eventTime, load2Switch, true);

  // simulation
  auto sim =
      setupSimulation(simName, simParams, systemTopology, logger, Domain::DP);
  sim.addEvent(disconnectLoad1);
  sim.addEvent(connectLoad2);
  sim.run();
}

void simulateSP(const SimulationParameters &simParams,
                const PowerSystemParameters &psParams) {
  String simName = "SP_simulation";
  Logger::setLogDir("logs/" + simName);
  auto logger = DataLogger::make(simName);

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
  infeedSource->setParameters(
      CPS::Math::polar(psParams.voltageLineToGround, 0.0));
  infeedSource->connect({SP::SimNode::GND, node1});
  auto infeedImpedance = SP::Ph1::PiLine::make("infeed_impedance");
  infeedImpedance->setParameters(psParams.infeedResistance,
                                 psParams.infeedInductance, 0);
  infeedImpedance->connect({node1, node4});
  logger->logAttribute(VariableNames::vInfeed,
                       node4->attribute(AttributeNames::v));
  logger->logAttribute(VariableNames::iInfeed,
                       infeedImpedance->attribute(AttributeNames::i));

  auto converter1 = SP::Ph1::VoltageSource::make("converter1");
  converter1->setParameters(
      CPS::Math::polar(psParams.voltageLineToGround, 0.0));
  converter1->connect({SP::SimNode::GND, node2});
  logger->logAttribute(VariableNames::vConverter1,
                       node2->attribute(AttributeNames::v));
  logger->logAttribute(VariableNames::iConverter1,
                       converter1->attribute(AttributeNames::i));

  auto line1 = SP::Ph1::PiLine::make("line1");
  line1->setParameters(psParams.line1Resistance, psParams.line1Inductance,
                       psParams.line1Capacitance);
  line1->connect({node2, node4});

  auto converter2 = SP::Ph1::VoltageSource::make("converter2");
  converter2->setParameters(
      CPS::Math::polar(psParams.voltageLineToGround, 0.0));
  converter2->connect({SP::SimNode::GND, node3});
  logger->logAttribute(VariableNames::vConverter2,
                       node3->attribute(AttributeNames::v));
  logger->logAttribute(VariableNames::iConverter2,
                       converter2->attribute(AttributeNames::i));

  auto line2 = SP::Ph1::PiLine::make("line2");
  line2->setParameters(psParams.line2Resistance, psParams.line2Inductance,
                       psParams.line2Capacitance);
  line2->connect({node3, node4});

  auto circuitBreaker = SP::Ph1::Switch::make("circuit_breaker");
  circuitBreaker->setParameters(SwitchConstants::openResistance,
                                SwitchConstants::closedResistance, true);
  circuitBreaker->connect({node4, node5});
  logger->logAttribute(VariableNames::vLoad,
                       node5->attribute(AttributeNames::v));
  logger->logAttribute(VariableNames::iLoad,
                       circuitBreaker->attribute(AttributeNames::i));

  auto load1 = SP::Ph1::Resistor::make("load");
  load1->setParameters(psParams.loadResistance1);
  load1->connect({node6, SP::SimNode::GND});

  auto load1Switch = SP::Ph1::Switch::make("load1_switch");
  load1Switch->setParameters(SwitchConstants::openResistance,
                             SwitchConstants::closedResistance, true);
  load1Switch->connect({node5, node6});

  auto load2 = SP::Ph1::Resistor::make("load2");
  load2->setParameters(psParams.loadResistance2);
  load2->connect({node7, SP::SimNode::GND});

  auto load2Switch = SP::Ph1::Switch::make("load2_switch");
  load2Switch->setParameters(SwitchConstants::openResistance,
                             SwitchConstants::closedResistance, false);
  load2Switch->connect({node5, node7});

  // topology
  auto systemNodeList =
      SystemNodeList{node1, node2, node3, node4, node5, node6, node7};
  auto componentList = SystemComponentList{
      infeedSource,   infeedImpedance, converter1,  line1, converter2, line2,
      circuitBreaker, load1,           load1Switch, load2, load2Switch};
  auto systemTopology =
      SystemTopology(psParams.frequency, systemNodeList, componentList);

  // events
  auto disconnectLoad1 =
      DPsim::SwitchEvent::make(simParams.eventTime, load1Switch, false);
  auto connectLoad2 =
      DPsim::SwitchEvent::make(simParams.eventTime, load2Switch, true);

  // simulation
  auto sim =
      setupSimulation(simName, simParams, systemTopology, logger, Domain::SP);
  sim.addEvent(disconnectLoad1);
  sim.addEvent(connectLoad2);
  sim.run();
}

int main() {
  SimulationParameters simParams;
  PowerSystemInputParameters psInputParams;
  PowerSystemParameters psParams =
      calculatePowerSystemParameters(psInputParams);

  simulateEMT(simParams, psParams);
  simulateDP(simParams, psParams);
  simulateSP(simParams, psParams);
  return 0;
}