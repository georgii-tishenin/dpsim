#include "../Examples.h"
#include <DPsim.h>
using namespace DPsim;
using namespace CPS;

namespace VariableNames {
constexpr const char *vInfeed = "vInfeed";
constexpr const char *iInfeed = "iInfeed";
constexpr const char *vLoad = "vLoad";
constexpr const char *iLoad1 = "iLoad1";
constexpr const char *iLoad2 = "iLoad2";
constexpr const char *iLine1 = "iLine1";
constexpr const char *iLine2 = "iLine2";
} // namespace VariableNames

namespace SwitchConstants {
constexpr double closedResistance = 1e-4;
constexpr double openResistance = 1e6;
} // namespace SwitchConstants

namespace AttributeNames {
constexpr const char *v = "v";
constexpr const char *i = "i_intf";
constexpr const char *id = "Irc_d";
constexpr const char *iq = "Irc_q";
constexpr const char *vd = "Vc_d";
constexpr const char *vq = "Vc_q";
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
  sim.doInitFromNodesAndTerminals(true);
  sim.setDomain(domain);
  sim.addLogger(logger);
  return sim;
}

void createEMTConverter(const std::shared_ptr<DataLogger> &logger,
                        const PowerSystemParameters &psParams,
                        CPS::SystemTopology &systemTopology,
                        const std::shared_ptr<EMT::SimNode> &node,
                        const std::string &name) {
  CIM::Examples::Grids::SGIB::ScenarioConfig scenario;
  auto converter = EMT::Ph3::AvVoltageSourceInverterDQ::make(
      name, name, Logger::Level::debug, true);
  converter->setParameters(scenario.systemOmega, scenario.pvNominalVoltage,
                           scenario.pvNominalActivePower,
                           scenario.pvNominalReactivePower);
  converter->setControllerParameters(
      1 * scenario.KpPLL, 1 * scenario.KiPLL, 1 * scenario.KpPowerCtrl,
      1 * scenario.KiPowerCtrl, 1 * scenario.KpCurrCtrl,
      1 * scenario.KiCurrCtrl, scenario.OmegaCutoff);
  converter->setFilterParameters(scenario.Lf, scenario.Cf, scenario.Rf,
                                 scenario.Rc);
  converter->setTransformerParameters(
      psParams.voltageLineToLine, scenario.pvNominalVoltage,
      scenario.transformerNominalPower,
      psParams.voltageLineToLine / scenario.pvNominalVoltage, 0, 0,
      scenario.transformerInductance, scenario.systemOmega);
  // converter->setInitialStateValues(scenario.pvNominalActivePower,
  //                                  scenario.pvNominalReactivePower,
  //                                  scenario.phi_dInit, scenario.phi_qInit,
  //                                  scenario.gamma_dInit, scenario.gamma_qInit);
  converter->withControl(true);

  converter->connect({node});
  systemTopology.addComponent(converter);
  logger->logAttribute("v" + name, node->attribute(AttributeNames::v));
  logger->logAttribute("id" + name, converter->attribute(AttributeNames::id));
  logger->logAttribute("iq" + name, converter->attribute(AttributeNames::iq));
  logger->logAttribute("vd" + name, converter->attribute(AttributeNames::vd));
  logger->logAttribute("vq" + name, converter->attribute(AttributeNames::vq));
}

void createDPConverter(const std::shared_ptr<DataLogger> &logger,
                        const PowerSystemParameters &psParams,
                        CPS::SystemTopology &systemTopology,
                        const std::shared_ptr<DP::SimNode> &node,
                        const std::string &name) {
  CIM::Examples::Grids::SGIB::ScenarioConfig scenario;
  auto converter = DP::Ph1::AvVoltageSourceInverterDQ::make(
      name, name, Logger::Level::debug, true);
  converter->setParameters(scenario.systemOmega, scenario.pvNominalVoltage,
                           scenario.pvNominalActivePower,
                           scenario.pvNominalReactivePower);
  converter->setControllerParameters(
      1 * scenario.KpPLL, 1 * scenario.KiPLL, 1 * scenario.KpPowerCtrl,
      1 * scenario.KiPowerCtrl, 1 * scenario.KpCurrCtrl,
      1 * scenario.KiCurrCtrl, scenario.OmegaCutoff);
  converter->setFilterParameters(scenario.Lf, scenario.Cf, scenario.Rf,
                                 scenario.Rc);
  converter->setTransformerParameters(
      psParams.voltageLineToLine, scenario.pvNominalVoltage,
      scenario.transformerNominalPower,
      psParams.voltageLineToLine/ scenario.pvNominalVoltage, 0, 0,
      scenario.transformerInductance);
  // converter->setInitialStateValues(scenario.pvNominalActivePower,
  //                                  scenario.pvNominalReactivePower,
  //                                  scenario.phi_dInit, scenario.phi_qInit,
  //                                  scenario.gamma_dInit, scenario.gamma_qInit);
  converter->withControl(true);

  converter->connect({node});
  systemTopology.addComponent(converter);
  logger->logAttribute("v" + name, node->attribute(AttributeNames::v));
  logger->logAttribute("id" + name, converter->attribute(AttributeNames::id));
  logger->logAttribute("iq" + name, converter->attribute(AttributeNames::iq));
  logger->logAttribute("vd" + name, converter->attribute(AttributeNames::vd));
  logger->logAttribute("vq" + name, converter->attribute(AttributeNames::vq));
}

void createSPConverter(const std::shared_ptr<DataLogger> &logger,
                        const PowerSystemParameters &psParams,
                        CPS::SystemTopology &systemTopology,
                        const std::shared_ptr<SP::SimNode> &node,
                        const std::string &name) {
  CIM::Examples::Grids::SGIB::ScenarioConfig scenario;
  auto converter = SP::Ph1::AvVoltageSourceInverterDQ::make(
      name, name, Logger::Level::debug, true);
  converter->setParameters(scenario.systemOmega, scenario.pvNominalVoltage,
                           scenario.pvNominalActivePower,
                           scenario.pvNominalReactivePower);
  converter->setControllerParameters(
      1 * scenario.KpPLL, 1 * scenario.KiPLL, 1 * scenario.KpPowerCtrl,
      1 * scenario.KiPowerCtrl, 1 * scenario.KpCurrCtrl,
      1 * scenario.KiCurrCtrl, scenario.OmegaCutoff);
  converter->setFilterParameters(scenario.Lf, scenario.Cf, scenario.Rf,
                                 scenario.Rc);
  converter->setTransformerParameters(
      psParams.voltageLineToLine, scenario.pvNominalVoltage,
      scenario.transformerNominalPower,
      psParams.voltageLineToLine/ scenario.pvNominalVoltage, 0, 0,
      scenario.transformerInductance);
  // converter->setInitialStateValues(scenario.pvNominalActivePower,
  //                                  scenario.pvNominalReactivePower,
  //                                  scenario.phi_dInit, scenario.phi_qInit,
  //                                  scenario.gamma_dInit, scenario.gamma_qInit);
  converter->withControl(true);

  converter->connect({node});
  systemTopology.addComponent(converter);
  logger->logAttribute("v" + name, node->attribute(AttributeNames::v));
  logger->logAttribute("id" + name, converter->attribute(AttributeNames::id));
  logger->logAttribute("iq" + name, converter->attribute(AttributeNames::iq));
  logger->logAttribute("vd" + name, converter->attribute(AttributeNames::vd));
  logger->logAttribute("vq" + name, converter->attribute(AttributeNames::vq));
}

void createEMTConverterAsVoltageSource(
    const std::shared_ptr<DataLogger> &logger,
    const PowerSystemParameters &psParams, CPS::SystemTopology &systemTopology,
    const std::shared_ptr<EMT::SimNode> &node, const std::string &name) {
  auto converter = EMT::Ph3::VoltageSource::make(name);
  converter->setParameters(
      CPS::Math::singlePhaseVariableToThreePhase(
          CPS::Math::polar(psParams.voltageLineToLine, 0.0)),
      psParams.frequency);
  converter->connect({EMT::SimNode::GND, node});
  systemTopology.addComponent(converter);
  logger->logAttribute("v" + name, node->attribute(AttributeNames::v));
}

void createDPConverterAsVoltageSource(const std::shared_ptr<DataLogger> &logger,
                                      const PowerSystemParameters &psParams,
                                      CPS::SystemTopology &systemTopology,
                                      const std::shared_ptr<DP::SimNode> &node,
                                      const std::string &name) {
  auto converter = DP::Ph1::VoltageSource::make(name);
  converter->setParameters(CPS::Math::polar(psParams.voltageLineToGround, 0.0));
  converter->connect({DP::SimNode::GND, node});
  systemTopology.addComponent(converter);
  logger->logAttribute("v" + name, node->attribute(AttributeNames::v));
}

void createSPConverterAsVoltageSource(const std::shared_ptr<DataLogger> &logger,
                                      const PowerSystemParameters &psParams,
                                      CPS::SystemTopology &systemTopology,
                                      const std::shared_ptr<SP::SimNode> &node,
                                      const std::string &name) {
  auto converter = SP::Ph1::VoltageSource::make(name);
  converter->setParameters(CPS::Math::polar(psParams.voltageLineToGround, 0.0));
  converter->connect({SP::SimNode::GND, node});
  systemTopology.addComponent(converter);
  logger->logAttribute("v" + name, node->attribute(AttributeNames::v));
}

void simulateEMT(const SimulationParameters &simParams,
                 const PowerSystemParameters &psParams) {
  String simName = "EMT_simulation";
  Logger::setLogDir("logs/" + simName);
  auto logger = DataLogger::make(simName);

  // nodes
  auto node1 = EMT::SimNode::make("node1", PhaseType::ABC);
  auto node2 = EMT::SimNode::make("node2", PhaseType::ABC);
  node2->setInitialVoltage(
      CPS::Math::singlePhaseVariableToThreePhase(CPS::Math::polar(0.0, 0.0)));
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

  auto line1 = EMT::Ph3::PiLine::make("line1");
  line1->setParameters(
      CPS::Math::singlePhaseParameterToThreePhase(psParams.line1Resistance),
      CPS::Math::singlePhaseParameterToThreePhase(psParams.line1Inductance),
      CPS::Math::singlePhaseParameterToThreePhase(psParams.line1Capacitance));
  line1->connect({node2, node4});
  logger->logAttribute(VariableNames::iLine1,
                       line1->attribute(AttributeNames::i));

  auto line2 = EMT::Ph3::PiLine::make("line2");
  line2->setParameters(
      CPS::Math::singlePhaseParameterToThreePhase(psParams.line2Resistance),
      CPS::Math::singlePhaseParameterToThreePhase(psParams.line2Inductance),
      CPS::Math::singlePhaseParameterToThreePhase(psParams.line2Capacitance));
  line2->connect({node3, node4});
  logger->logAttribute(VariableNames::iLine2,
                       line2->attribute(AttributeNames::i));

  auto circuitBreaker = EMT::Ph3::Switch::make("circuit_breaker");
  circuitBreaker->setParameters(CPS::Math::singlePhaseParameterToThreePhase(
                                    SwitchConstants::openResistance),
                                CPS::Math::singlePhaseParameterToThreePhase(
                                    SwitchConstants::closedResistance),
                                true);
  circuitBreaker->connect({node4, node5});
  logger->logAttribute(VariableNames::vLoad,
                       node5->attribute(AttributeNames::v));

  auto load1 = EMT::Ph3::Resistor::make("load1");
  load1->setParameters(
      CPS::Math::singlePhaseParameterToThreePhase(psParams.loadResistance1));
  load1->connect({node6, EMT::SimNode::GND});
  logger->logAttribute(VariableNames::iLoad1,
                       load1->attribute(AttributeNames::i));

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
  logger->logAttribute(VariableNames::iLoad2,
                       load2->attribute(AttributeNames::i));

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
      infeedSource, infeedImpedance, line1, line2,      circuitBreaker,
      load1,        load1Switch,     load2, load2Switch};
  auto systemTopology =
      SystemTopology(psParams.frequency, systemNodeList, componentList);

  createEMTConverter(logger, psParams, systemTopology, node2, "Converter1");
//   createEMTConverterAsVoltageSource(logger, psParams, systemTopology, node2,
//                                     "Converter1");
  createEMTConverterAsVoltageSource(logger, psParams, systemTopology, node3,
                                    "Converter2");

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

  auto line1 = DP::Ph1::PiLine::make("line1");
  line1->setParameters(psParams.line1Resistance, psParams.line1Inductance,
                       psParams.line1Capacitance);
  line1->connect({node2, node4});
  logger->logAttribute(VariableNames::iLine1,
                       line1->attribute(AttributeNames::i));

  auto line2 = DP::Ph1::PiLine::make("line2");
  line2->setParameters(psParams.line2Resistance, psParams.line2Inductance,
                       psParams.line2Capacitance);
  line2->connect({node3, node4});
  logger->logAttribute(VariableNames::iLine2,
                       line2->attribute(AttributeNames::i));

  auto circuitBreaker = DP::Ph1::Switch::make("circuit_breaker");
  circuitBreaker->setParameters(SwitchConstants::openResistance,
                                SwitchConstants::closedResistance, true);
  circuitBreaker->connect({node4, node5});
  logger->logAttribute(VariableNames::vLoad,
                       node5->attribute(AttributeNames::v));

  auto load1 = DP::Ph1::Resistor::make("load");
  load1->setParameters(psParams.loadResistance1);
  load1->connect({node6, DP::SimNode::GND});
  logger->logAttribute(VariableNames::iLoad1,
                       load1->attribute(AttributeNames::i));

  auto load1Switch = DP::Ph1::Switch::make("load1_switch");
  load1Switch->setParameters(SwitchConstants::openResistance,
                             SwitchConstants::closedResistance, true);
  load1Switch->connect({node5, node6});

  auto load2 = DP::Ph1::Resistor::make("load2");
  load2->setParameters(psParams.loadResistance2);
  load2->connect({node7, DP::SimNode::GND});
  logger->logAttribute(VariableNames::iLoad2,
                       load2->attribute(AttributeNames::i));

  auto load2Switch = DP::Ph1::Switch::make("load2_switch");
  load2Switch->setParameters(SwitchConstants::openResistance,
                             SwitchConstants::closedResistance, false);
  load2Switch->connect({node5, node7});

  // topology
  auto systemNodeList =
      SystemNodeList{node1, node2, node3, node4, node5, node6, node7};
  auto componentList = SystemComponentList{
      infeedSource, infeedImpedance, line1, line2,      circuitBreaker,
      load1,        load1Switch,     load2, load2Switch};
  auto systemTopology =
      SystemTopology(psParams.frequency, systemNodeList, componentList);

//   createDPConverterAsVoltageSource(logger, psParams, systemTopology, node2,
//                                    "Converter1");
  createDPConverter(logger, psParams, systemTopology, node2, "Converter1");
  createDPConverterAsVoltageSource(logger, psParams, systemTopology, node3,
                                   "Converter2");

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

  auto line1 = SP::Ph1::PiLine::make("line1");
  line1->setParameters(psParams.line1Resistance, psParams.line1Inductance,
                       psParams.line1Capacitance);
  line1->connect({node2, node4});
  logger->logAttribute(VariableNames::iLine1,
                       line1->attribute(AttributeNames::i));

  auto line2 = SP::Ph1::PiLine::make("line2");
  line2->setParameters(psParams.line2Resistance, psParams.line2Inductance,
                       psParams.line2Capacitance);
  line2->connect({node3, node4});
  logger->logAttribute(VariableNames::iLine2,
                       line2->attribute(AttributeNames::i));

  auto circuitBreaker = SP::Ph1::Switch::make("circuit_breaker");
  circuitBreaker->setParameters(SwitchConstants::openResistance,
                                SwitchConstants::closedResistance, true);
  circuitBreaker->connect({node4, node5});
  logger->logAttribute(VariableNames::vLoad,
                       node5->attribute(AttributeNames::v));

  auto load1 = SP::Ph1::Resistor::make("load");
  load1->setParameters(psParams.loadResistance1);
  load1->connect({node6, SP::SimNode::GND});
  logger->logAttribute(VariableNames::iLoad1,
                       load1->attribute(AttributeNames::i));

  auto load1Switch = SP::Ph1::Switch::make("load1_switch");
  load1Switch->setParameters(SwitchConstants::openResistance,
                             SwitchConstants::closedResistance, true);
  load1Switch->connect({node5, node6});

  auto load2 = SP::Ph1::Resistor::make("load2");
  load2->setParameters(psParams.loadResistance2);
  load2->connect({node7, SP::SimNode::GND});
  logger->logAttribute(VariableNames::iLoad2,
                       load2->attribute(AttributeNames::i));

  auto load2Switch = SP::Ph1::Switch::make("load2_switch");
  load2Switch->setParameters(SwitchConstants::openResistance,
                             SwitchConstants::closedResistance, false);
  load2Switch->connect({node5, node7});

  // topology
  auto systemNodeList =
      SystemNodeList{node1, node2, node3, node4, node5, node6, node7};
  auto componentList = SystemComponentList{
      infeedSource, infeedImpedance, line1, line2,      circuitBreaker,
      load1,        load1Switch,     load2, load2Switch};
  auto systemTopology =
      SystemTopology(psParams.frequency, systemNodeList, componentList);

//   createSPConverterAsVoltageSource(logger, psParams, systemTopology, node2,
//                                    "Converter1");
  createSPConverter(logger, psParams, systemTopology, node2, "Converter1");
  createSPConverterAsVoltageSource(logger, psParams, systemTopology, node3,
                                   "Converter2");

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