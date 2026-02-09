#include "../Examples.h"
#include <DPsim.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>

using namespace DPsim;
using namespace CPS;

namespace HVDCWise {

namespace VariableNames {
constexpr const char *vInfeed = "vInfeed";
constexpr const char *iInfeedStrong = "iInfeedStrong";
constexpr const char *iInfeedWeak = "iInfeedWeak";
constexpr const char *fInfeed = "fInfeed";

constexpr const char *vLoad = "vLoad";
constexpr const char *iLoad1 = "iLoad1";
constexpr const char *iLoad2 = "iLoad2";
constexpr const char *iLine1 = "iLine1";
constexpr const char *iLine2 = "iLine2";
constexpr const char *iFault = "iFault";
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
constexpr const char *f = "f_src";
constexpr const char *pref = "P_ref";
constexpr const char *qref = "Q_ref";
constexpr const char *pllOut = "pll_output";
} // namespace AttributeNames

static inline double clamp01(double x) {
  if (x < 0.0)
    return 0.0;
  if (x > 1.0)
    return 1.0;
  return x;
}

struct SimulationParameters {
  double timeStep = 1e-4;
  double eventTime = 3.8;
  double finalTime = 4.2;

  // Startup ramp of converter P/Q references (to avoid numerical issues at t=0)
  bool enableStartupRampPQ = true;

  // hold Pref/Qref = 0 for this time BEFORE ramp starts (lets PLL settle)
  double startupPQZeroHoldTime = 0.5; // seconds, Pref=Qref=0 for t in [0, hold)

  // ramp from 0 -> target during [hold, hold + startupRampDuration]
  double startupRampDuration = 0.5; // seconds

  // frequency ramp parameters
  double frequencyRampDuration = 0.1;
  double rocof = -10; // in Hz/s

  // frequency step parameters
  double frequencyStepDelta = -1.0; // in Hz

  double prefStepFactor = 5.0;

  // ---------------- Load-bus fault parameters ----------------
  // Fault is modeled as a shunt branch at load bus (node5) to GND.
  double faultDuration = 0.02;  // seconds, clear at eventTime + faultDuration
  double faultResistance = 1.0; // Ohm when fault is "ON" (closed branch)

  // ---------------- Infeed SCR-step parameters ----------------
  // Two parallel infeed impedance branches (strong/weak),
  // each gated by a series switch. At eventTime we open strong and close weak.
  // If factor > 1 => weaker grid (lower SCR). If factor < 1 => stronger grid (higher SCR).
  double infeedImpedanceStepFactor = 3.0;

  // ---------------- Infeed voltage ANGLE-step parameters ----------------
  // Single angle step (degrees) applied at eventTime to NetworkInjection voltage phasor.
  // Positive => advance angle.
  double infeedVoltageAngleStepDeg = 10.0;

  // ---------------- Optional converter3 (same PCC as converter1) ----------------
  bool enableConverter3 = false; // if true: add Converter3 at node2 with same P/Q as Converter1
};

enum class PowerSystemEventType {
  None,
  LoadStep,
  LoadBusFault,
  InfeedSCRStep,
  InfeedFrequencyRamp,
  InfeedFrequencyStep,
  Converter1PrefStep,
  InfeedVoltageAngleStep
};

struct PowerSystemInputParameters {
  double frequency = 50;
  double baseVoltageLineToLine = 110e3;
  double baseThreePhasePower = 100e6;

  // infeed parameters
  double infeedResistanceInPerUnit = 0.01;
  double infeedReactanceInPerUnit = 0.1;

  // coefficients for line parameters
  double lineLengthCoefficient = 1;     // 3, 5;
  double lineResistanceCoefficient = 1; // 0.1;

  // line1 parameters
  double line1LengthInKm = 80 * lineLengthCoefficient;
  double line1ResistancePerKm = 0.1 * lineResistanceCoefficient;
  double line1ReactancePerKm = 0.4;
  double line1CapacitancePerKm = 1e-8;

  // line2 parameters
  double line2LengthInKm = 20 * lineLengthCoefficient;
  double line2ResistancePerKm = 0.2 * lineResistanceCoefficient;
  double line2ReactancePerKm = 0.4;
  double lin2CapacitancePerKm = 1e-8;

  // load parameters
  double load1InPerUnit = 0.5;
  double load2InPerUnit = 1.0;

  // converter1 parameters
  double converter1PinPerUnit = 0.1;
  double converter1QinPerUnit = 0;

  // converter2 parameters
  double converter2PinPerUnit = 0.05;
  double converter2QinPerUnit = 0;
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
  double converter1P;
  double converter1Q;
  double converter2P;
  double converter2Q;

  PowerSystemParameters(double freq, double voltLineToLine,
                        double voltLineToGround, double infeedRes,
                        double infeedInd, double line1Res, double line1Ind,
                        double line1Cap, double line2Res, double line2Ind,
                        double line2Cap, double loadRes1, double loadRes2,
                        double conv1P, double conv1Q, double conv2P,
                        double conv2Q)
      : frequency(freq), voltageLineToLine(voltLineToLine),
        voltageLineToGround(voltLineToGround), infeedResistance(infeedRes),
        infeedInductance(infeedInd), line1Resistance(line1Res),
        line1Inductance(line1Ind), line1Capacitance(line1Cap),
        line2Resistance(line2Res), line2Inductance(line2Ind),
        line2Capacitance(line2Cap), loadResistance1(loadRes1),
        loadResistance2(loadRes2), converter1P(conv1P), converter1Q(conv1Q),
        converter2P(conv2P), converter2Q(conv2Q) {}
};

PowerSystemParameters
calculatePowerSystemParameters(const PowerSystemInputParameters &inputParams) {
  double voltageLineToGround = inputParams.baseVoltageLineToLine / std::sqrt(3);
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

  double converter1P =
      inputParams.converter1PinPerUnit * inputParams.baseThreePhasePower;
  double converter1Q =
      inputParams.converter1QinPerUnit * inputParams.baseThreePhasePower;
  double converter2P =
      inputParams.converter2PinPerUnit * inputParams.baseThreePhasePower;
  double converter2Q =
      inputParams.converter2QinPerUnit * inputParams.baseThreePhasePower;

  return PowerSystemParameters(
      inputParams.frequency, inputParams.baseVoltageLineToLine,
      voltageLineToGround, infeedResistance, infeedInductance, line1Resistance,
      line1Inductance, line1Capacitance, line2Resistance, line2Inductance,
      line2Capacitance, loadResistance1, loadResistance2, converter1P,
      converter1Q, converter2P, converter2Q);
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

template <typename Hook>
static void runStepped(DPsim::Simulation &sim, Hook &&hook) {
  sim.initialize();
  sim.start();
  while (sim.time() < sim.finalTime()) {
    hook(sim); // hook runs at every step (before step())
    sim.step();
  }
  sim.stop();
}

// -------------------- Converter handles --------------------

struct EMTConverterHandle {
  std::shared_ptr<EMT::Ph3::AvVoltageSourceInverterDQ> conv;
  double sysOmega;
  double sysVoltNom;
  double pFinal;
  double qFinal;
};

struct DPConverterHandle {
  std::shared_ptr<DP::Ph1::AvVoltageSourceInverterDQ> conv;
  double sysOmega;
  double sysVoltNom;
  double pFinal;
  double qFinal;
};

struct SPConverterHandle {
  std::shared_ptr<SP::Ph1::AvVoltageSourceInverterDQ> conv;
  double sysOmega;
  double sysVoltNom;
  double pFinal;
  double qFinal;
};

// -------------------- Converter creation --------------------

EMTConverterHandle createEMTConverter(const std::shared_ptr<DataLogger> &logger,
                                      const PowerSystemParameters &psParams,
                                      CPS::SystemTopology &systemTopology,
                                      const std::shared_ptr<EMT::SimNode> &node,
                                      int converterNumber,
                                      bool startRampEnabled = false) {
  CIM::Examples::Grids::SGIB::ScenarioConfig scenario;

  double converterP_final = 0.0;
  double converterQ_final = 0.0;
  switch (converterNumber) {
  case 1:
    converterP_final = psParams.converter1P;
    converterQ_final = psParams.converter1Q;
    break;
  case 2:
    converterP_final = psParams.converter2P;
    converterQ_final = psParams.converter2Q;
    break;
  case 3: // same as converter1 by default
    converterP_final = psParams.converter1P;
    converterQ_final = psParams.converter1Q;
    break;
  default:
    throw std::invalid_argument("Unsupported converter number: " +
                                std::to_string(converterNumber));
  }

  auto converter = EMT::Ph3::AvVoltageSourceInverterDQ::make(
      "Converter" + std::to_string(converterNumber),
      "Converter" + std::to_string(converterNumber), Logger::Level::debug,
      true);

  const double sysOmega = 2.0 * M_PI * psParams.frequency;
  const double sysVoltNom = scenario.pvNominalVoltage;

  const double converterP_init = startRampEnabled ? 0.0 : converterP_final;
  const double converterQ_init = startRampEnabled ? 0.0 : converterQ_final;

  converter->setParameters(sysOmega, sysVoltNom, converterP_init,
                           converterQ_init);
  converter->setControllerParameters(
      1 * scenario.KpPLL, 1 * scenario.KiPLL, 1 * scenario.KpPowerCtrl,
      1 * scenario.KiPowerCtrl, 1 * scenario.KpCurrCtrl,
      1 * scenario.KiCurrCtrl, scenario.OmegaCutoff);
  converter->setFilterParameters(scenario.Lf, scenario.Cf, scenario.Rf,
                                 scenario.Rc);
  converter->setTransformerParameters(
      psParams.voltageLineToLine, sysVoltNom, scenario.transformerNominalPower,
      psParams.voltageLineToLine / scenario.pvNominalVoltage, 0, 0,
      scenario.transformerInductance, scenario.systemOmega);

  converter->withControl(true);
  converter->connect({node});
  systemTopology.addComponent(converter);

  logger->logAttribute("vConverter" + std::to_string(converterNumber),
                       node->attribute(AttributeNames::v));
  logger->logAttribute("idConverter" + std::to_string(converterNumber),
                       converter->attribute(AttributeNames::id));
  logger->logAttribute("iqConverter" + std::to_string(converterNumber),
                       converter->attribute(AttributeNames::iq));
  logger->logAttribute("vdConverter" + std::to_string(converterNumber),
                       converter->attribute(AttributeNames::vd));
  logger->logAttribute("vqConverter" + std::to_string(converterNumber),
                       converter->attribute(AttributeNames::vq));
  logger->logAttribute("PrefConverter" + std::to_string(converterNumber),
                       converter->attribute(AttributeNames::pref));
  logger->logAttribute("QrefConverter" + std::to_string(converterNumber),
                       converter->attribute(AttributeNames::qref));
  logger->logAttribute("pllOutputConverter" + std::to_string(converterNumber),
                       converter->attribute(AttributeNames::pllOut));

  return {converter, sysOmega, sysVoltNom, converterP_final, converterQ_final};
}

DPConverterHandle createDPConverter(const std::shared_ptr<DataLogger> &logger,
                                    const PowerSystemParameters &psParams,
                                    CPS::SystemTopology &systemTopology,
                                    const std::shared_ptr<DP::SimNode> &node,
                                    int converterNumber,
                                    bool startRampEnabled = false) {
  CIM::Examples::Grids::SGIB::ScenarioConfig scenario;

  double converterP_final = 0.0;
  double converterQ_final = 0.0;
  switch (converterNumber) {
  case 1:
    converterP_final = psParams.converter1P;
    converterQ_final = psParams.converter1Q;
    break;
  case 2:
    converterP_final = psParams.converter2P;
    converterQ_final = psParams.converter2Q;
    break;
  case 3: // same as converter1 by default
    converterP_final = psParams.converter1P;
    converterQ_final = psParams.converter1Q;
    break;
  default:
    throw std::invalid_argument("Unsupported converter number: " +
                                std::to_string(converterNumber));
  }

  auto converter = DP::Ph1::AvVoltageSourceInverterDQ::make(
      "Converter" + std::to_string(converterNumber),
      "Converter" + std::to_string(converterNumber), Logger::Level::debug,
      true);

  const double sysOmega = 2.0 * M_PI * psParams.frequency;
  const double sysVoltNom = scenario.pvNominalVoltage;

  const double converterP_init = startRampEnabled ? 0.0 : converterP_final;
  const double converterQ_init = startRampEnabled ? 0.0 : converterQ_final;

  converter->setParameters(sysOmega, sysVoltNom, converterP_init,
                           converterQ_init);
  converter->setControllerParameters(
      1 * scenario.KpPLL, 1 * scenario.KiPLL, 1 * scenario.KpPowerCtrl,
      1 * scenario.KiPowerCtrl, 1 * scenario.KpCurrCtrl,
      1 * scenario.KiCurrCtrl, scenario.OmegaCutoff);
  converter->setFilterParameters(scenario.Lf, scenario.Cf, scenario.Rf,
                                 scenario.Rc);
  converter->setTransformerParameters(psParams.voltageLineToLine, sysVoltNom,
                                      scenario.transformerNominalPower,
                                      psParams.voltageLineToLine / sysVoltNom,
                                      0, 0, scenario.transformerInductance);

  converter->withControl(true);
  converter->connect({node});
  systemTopology.addComponent(converter);

  logger->logAttribute("vConverter" + std::to_string(converterNumber),
                       node->attribute(AttributeNames::v));
  logger->logAttribute("idConverter" + std::to_string(converterNumber),
                       converter->attribute(AttributeNames::id));
  logger->logAttribute("iqConverter" + std::to_string(converterNumber),
                       converter->attribute(AttributeNames::iq));
  logger->logAttribute("vdConverter" + std::to_string(converterNumber),
                       converter->attribute(AttributeNames::vd));
  logger->logAttribute("vqConverter" + std::to_string(converterNumber),
                       converter->attribute(AttributeNames::vq));
  logger->logAttribute("PrefConverter" + std::to_string(converterNumber),
                       converter->attribute(AttributeNames::pref));
  logger->logAttribute("QrefConverter" + std::to_string(converterNumber),
                       converter->attribute(AttributeNames::qref));
  logger->logAttribute("pllOutputConverter" + std::to_string(converterNumber),
                       converter->attribute(AttributeNames::pllOut));

  return {converter, sysOmega, sysVoltNom, converterP_final, converterQ_final};
}

SPConverterHandle createSPConverter(const std::shared_ptr<DataLogger> &logger,
                                    const PowerSystemParameters &psParams,
                                    CPS::SystemTopology &systemTopology,
                                    const std::shared_ptr<SP::SimNode> &node,
                                    int converterNumber,
                                    bool startRampEnabled = false) {
  CIM::Examples::Grids::SGIB::ScenarioConfig scenario;

  double converterP_final = 0.0;
  double converterQ_final = 0.0;
  switch (converterNumber) {
  case 1:
    converterP_final = psParams.converter1P;
    converterQ_final = psParams.converter1Q;
    break;
  case 2:
    converterP_final = psParams.converter2P;
    converterQ_final = psParams.converter2Q;
    break;
  case 3: // same as converter1 by default
    converterP_final = psParams.converter1P;
    converterQ_final = psParams.converter1Q;
    break;
  default:
    throw std::invalid_argument("Unsupported converter number: " +
                                std::to_string(converterNumber));
  }

  auto converter = SP::Ph1::AvVoltageSourceInverterDQ::make(
      "Converter" + std::to_string(converterNumber),
      "Converter" + std::to_string(converterNumber), Logger::Level::debug,
      true);

  const double sysOmega = 2.0 * M_PI * psParams.frequency;
  const double sysVoltNom = scenario.pvNominalVoltage;

  const double converterP_init = startRampEnabled ? 0.0 : converterP_final;
  const double converterQ_init = startRampEnabled ? 0.0 : converterQ_final;

  converter->setParameters(sysOmega, sysVoltNom, converterP_init,
                           converterQ_init);
  converter->setControllerParameters(
      1 * scenario.KpPLL, 1 * scenario.KiPLL, 1 * scenario.KpPowerCtrl,
      1 * scenario.KiPowerCtrl, 1 * scenario.KpCurrCtrl,
      1 * scenario.KiCurrCtrl, scenario.OmegaCutoff);
  converter->setFilterParameters(scenario.Lf, scenario.Cf, scenario.Rf,
                                 scenario.Rc);
  converter->setTransformerParameters(psParams.voltageLineToLine, sysVoltNom,
                                      scenario.transformerNominalPower,
                                      psParams.voltageLineToLine / sysVoltNom,
                                      0, 0, scenario.transformerInductance);

  converter->withControl(true);
  converter->connect({node});
  systemTopology.addComponent(converter);

  logger->logAttribute("vConverter" + std::to_string(converterNumber),
                       node->attribute(AttributeNames::v));
  logger->logAttribute("idConverter" + std::to_string(converterNumber),
                       converter->attribute(AttributeNames::id));
  logger->logAttribute("iqConverter" + std::to_string(converterNumber),
                       converter->attribute(AttributeNames::iq));
  logger->logAttribute("vdConverter" + std::to_string(converterNumber),
                       converter->attribute(AttributeNames::vd));
  logger->logAttribute("vqConverter" + std::to_string(converterNumber),
                       converter->attribute(AttributeNames::vq));
  logger->logAttribute("PrefConverter" + std::to_string(converterNumber),
                       converter->attribute(AttributeNames::pref));
  logger->logAttribute("QrefConverter" + std::to_string(converterNumber),
                       converter->attribute(AttributeNames::qref));
  logger->logAttribute("pllOutputConverter" + std::to_string(converterNumber),
                       converter->attribute(AttributeNames::pllOut));

  return {converter, sysOmega, sysVoltNom, converterP_final, converterQ_final};
}

// --------- Simulation functions ---------

void simulateEMT(const SimulationParameters &simParams,
                 const PowerSystemParameters &psParams,
                 const SystemTopology &systemTopologyPF,
                 const PowerSystemEventType &psEvent) {
  String simName = "EMT_simulation";
  Logger::setLogDir("logs/" + simName);
  auto logger = DataLogger::make(simName);

  // ---------------- Nodes ----------------
  auto node1 = EMT::SimNode::make("node1", PhaseType::ABC);
  auto node2 = EMT::SimNode::make("node2", PhaseType::ABC);
  auto node3 = EMT::SimNode::make("node3", PhaseType::ABC);
  auto node4 = EMT::SimNode::make("node4", PhaseType::ABC);
  auto node5 = EMT::SimNode::make("node5", PhaseType::ABC);
  auto node6 = EMT::SimNode::make("node6", PhaseType::ABC);
  auto node7 = EMT::SimNode::make("node7", PhaseType::ABC);

  // extra nodes for robust SCR step (series switches need intermediate nodes)
  auto node1s = EMT::SimNode::make("node1_strong", PhaseType::ABC);
  auto node1w = EMT::SimNode::make("node1_weak", PhaseType::ABC);

  // ---------------- Components ----------------
  auto infeedSource = EMT::Ph3::NetworkInjection::make("infeed_source");
  infeedSource->connect({node1});

  // Robust SCR step: strong/weak infeed branches with series switches
  const double kZ = std::max(1e-9, simParams.infeedImpedanceStepFactor);
  const double Rstrong = psParams.infeedResistance;
  const double Lstrong = psParams.infeedInductance;
  const double Rweak = psParams.infeedResistance * kZ;
  const double Lweak = psParams.infeedInductance * kZ;

  auto infeedSwStrong = EMT::Ph3::Switch::make("infeed_sw_strong");
  infeedSwStrong->setParameters(CPS::Math::singlePhaseParameterToThreePhase(
                                    SwitchConstants::openResistance),
                                CPS::Math::singlePhaseParameterToThreePhase(
                                    SwitchConstants::closedResistance),
                                true // initially CLOSED (strong grid)
  );
  infeedSwStrong->connect({node1, node1s});

  auto infeedZStrong = EMT::Ph3::PiLine::make("infeed_impedance_strong");
  infeedZStrong->setParameters(
      CPS::Math::singlePhaseParameterToThreePhase(Rstrong),
      CPS::Math::singlePhaseParameterToThreePhase(Lstrong),
      CPS::Math::singlePhaseParameterToThreePhase(0.0));
  infeedZStrong->connect({node1s, node4});

  auto infeedSwWeak = EMT::Ph3::Switch::make("infeed_sw_weak");
  infeedSwWeak->setParameters(CPS::Math::singlePhaseParameterToThreePhase(
                                  SwitchConstants::openResistance),
                              CPS::Math::singlePhaseParameterToThreePhase(
                                  SwitchConstants::closedResistance),
                              false // initially OPEN (weak grid disconnected)
  );
  infeedSwWeak->connect({node1, node1w});

  auto infeedZWeak = EMT::Ph3::PiLine::make("infeed_impedance_weak");
  infeedZWeak->setParameters(CPS::Math::singlePhaseParameterToThreePhase(Rweak),
                             CPS::Math::singlePhaseParameterToThreePhase(Lweak),
                             CPS::Math::singlePhaseParameterToThreePhase(0.0));
  infeedZWeak->connect({node1w, node4});

  logger->logAttribute(VariableNames::vInfeed,
                       node4->attribute(AttributeNames::v));
  logger->logAttribute(VariableNames::iInfeedStrong,
                       infeedZStrong->attribute(AttributeNames::i));
  logger->logAttribute(VariableNames::iInfeedWeak,
                       infeedZWeak->attribute(AttributeNames::i));
  logger->logAttribute(VariableNames::fInfeed,
                       infeedSource->attribute(AttributeNames::f));

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

  // load-bus fault as shunt switch node5 -> GND
  auto loadBusFault = EMT::Ph3::Switch::make("load_bus_fault");
  loadBusFault->setParameters(
      CPS::Math::singlePhaseParameterToThreePhase(
          SwitchConstants::openResistance),
      CPS::Math::singlePhaseParameterToThreePhase(simParams.faultResistance),
      false // initially open (no fault)
  );
  loadBusFault->connect({node5, EMT::SimNode::GND});
  logger->logAttribute(VariableNames::iFault,
                       loadBusFault->attribute(AttributeNames::i));

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

  // ---------------- Topology ----------------
  auto systemNodeList = SystemNodeList{node1, node1s, node1w, node2, node3,
                                       node4, node5,  node6,  node7};

  auto componentList = SystemComponentList{
      infeedSource, infeedSwStrong, infeedZStrong,  infeedSwWeak, infeedZWeak,
      line1,        line2,          circuitBreaker, loadBusFault, load1,
      load1Switch,  load2,          load2Switch};

  auto systemTopology =
      SystemTopology(psParams.frequency, systemNodeList, componentList);

  const bool doStartupRamp = simParams.enableStartupRampPQ;
  auto conv1 = createEMTConverter(logger, psParams, systemTopology, node2, 1,
                                  doStartupRamp);
  auto conv2 = createEMTConverter(logger, psParams, systemTopology, node3, 2,
                                  doStartupRamp);

  // Optional converter3 at same PCC as converter1 (node2)
  EMTConverterHandle conv3{nullptr, 0.0, 0.0, 0.0, 0.0};
  const bool hasConv3 = simParams.enableConverter3;
  if (hasConv3) {
    conv3 = createEMTConverter(logger, psParams, systemTopology, node2, 3, doStartupRamp);
  }

  // ---------------- Simulation ----------------
  systemTopology.initWithPowerflow(systemTopologyPF, Domain::EMT);
  auto sim =
      setupSimulation(simName, simParams, systemTopology, logger, Domain::EMT);

  // ---------------- Events ----------------
  switch (psEvent) {
  case PowerSystemEventType::LoadStep: {
    auto disconnectLoad1 =
        DPsim::SwitchEvent3Ph::make(simParams.eventTime, load1Switch, false);
    auto connectLoad2 =
        DPsim::SwitchEvent3Ph::make(simParams.eventTime, load2Switch, true);
    sim.addEvent(disconnectLoad1);
    sim.addEvent(connectLoad2);
    break;
  }
  case PowerSystemEventType::LoadBusFault: {
    const double tOn = simParams.eventTime;
    const double tOff =
        simParams.eventTime + std::max(0.0, simParams.faultDuration);
    auto faultOn = DPsim::SwitchEvent3Ph::make(tOn, loadBusFault, true);
    auto faultOff = DPsim::SwitchEvent3Ph::make(tOff, loadBusFault, false);
    sim.addEvent(faultOn);
    sim.addEvent(faultOff);
    break;
  }
  case PowerSystemEventType::InfeedSCRStep: {
    // Open strong branch and close weak branch at eventTime
    const double t = simParams.eventTime;
    auto openStrong = DPsim::SwitchEvent3Ph::make(t, infeedSwStrong, false);
    auto closeWeak = DPsim::SwitchEvent3Ph::make(t, infeedSwWeak, true);
    sim.addEvent(openStrong);
    sim.addEvent(closeWeak);
    break;
  }
  case PowerSystemEventType::InfeedFrequencyRamp: {
    infeedSource->setParameters(
        CPS::Math::singlePhaseVariableToThreePhase(psParams.voltageLineToLine),
        psParams.frequency, simParams.rocof, simParams.eventTime,
        simParams.frequencyRampDuration, false);
    break;
  }
  case PowerSystemEventType::InfeedFrequencyStep: {
    infeedSource->setParameters(
        CPS::Math::singlePhaseVariableToThreePhase(psParams.voltageLineToLine),
        psParams.frequency, simParams.frequencyStepDelta / simParams.timeStep,
        simParams.eventTime, simParams.timeStep, false);
    break;
  }
  case PowerSystemEventType::InfeedVoltageAngleStep:
  case PowerSystemEventType::Converter1PrefStep:
  case PowerSystemEventType::None:
  default:
    break;
  }

  // ---- Hook: startup ramp + optional Converter1 Pref step + optional infeed angle step ----
  const double holdT = std::max(0.0, simParams.startupPQZeroHoldTime);
  const double rampDur = std::max(0.0, simParams.startupRampDuration);
  const double rampEndT = holdT + rampDur;

  bool rampDone = !doStartupRamp;
  bool printedHold = false;
  bool printedDone = false;

  const double prefStepTime = simParams.eventTime;
  const double newPref = psParams.converter1P * simParams.prefStepFactor;
  bool prefStepApplied = false;

  // Infeed voltage angle step state
  const double angleStepTime = simParams.eventTime;
  const double deltaRad = simParams.infeedVoltageAngleStepDeg * M_PI / 180.0;
  bool angleStepApplied = false;

  runStepped(sim, [&](DPsim::Simulation &s) {
    const double t = s.time();

    if (doStartupRamp && !rampDone) {
      double alpha = 0.0;

      if (t < holdT) {
        alpha = 0.0;
        if (!printedHold && t >= (0.0 + 0.5 * simParams.timeStep)) {
          printedHold = true;
          std::cout << "[HOOK][EMT] holding Pref/Qref at 0 for " << holdT
                    << " s\n";
        }
      } else if (rampDur <= 0.0) {
        alpha = 1.0;
      } else {
        alpha = clamp01((t - holdT) / rampDur);
      }

      conv1.conv->setParameters(conv1.sysOmega, conv1.sysVoltNom,
                                alpha * conv1.pFinal, alpha * conv1.qFinal);
      conv2.conv->setParameters(conv2.sysOmega, conv2.sysVoltNom,
                                alpha * conv2.pFinal, alpha * conv2.qFinal);
      if (hasConv3) {
        conv3.conv->setParameters(conv3.sysOmega, conv3.sysVoltNom,
                                  alpha * conv3.pFinal, alpha * conv3.qFinal);
      }

      if (!printedDone && t >= (rampEndT - 0.5 * simParams.timeStep)) {
        conv1.conv->setParameters(conv1.sysOmega, conv1.sysVoltNom,
                                  conv1.pFinal, conv1.qFinal);
        conv2.conv->setParameters(conv2.sysOmega, conv2.sysVoltNom,
                                  conv2.pFinal, conv2.qFinal);
        if (hasConv3) {
          conv3.conv->setParameters(conv3.sysOmega, conv3.sysVoltNom,
                                    conv3.pFinal, conv3.qFinal);
        }
        rampDone = true;
        printedDone = true;
        std::cout << "[HOOK][EMT] t=" << t
                  << " finished startup hold+ramp: hold=" << holdT
                  << "s, ramp=" << rampDur << "s\n";
      }
    }

    if (psEvent == PowerSystemEventType::Converter1PrefStep &&
        !prefStepApplied && t >= (prefStepTime - 0.5 * simParams.timeStep)) {

      conv1.conv->setParameters(conv1.sysOmega, conv1.sysVoltNom, newPref,
                                conv1.qFinal);

      prefStepApplied = true;
      std::cout << "[HOOK][EMT] t=" << t << " set Converter1 Pref=" << newPref
                << " W\n";
    }

    if (psEvent == PowerSystemEventType::InfeedVoltageAngleStep &&
        !angleStepApplied && t >= (angleStepTime - 0.5 * simParams.timeStep)) {

      const Complex Vref = std::polar(psParams.voltageLineToLine, deltaRad);

      // Use the existing 6-arg overload (same as frequency ramp/step) with rocof=0.
      // The complex Vref carries the desired phase angle.
      infeedSource->setParameters(
          CPS::Math::singlePhaseVariableToThreePhase(Vref), psParams.frequency,
          0.0,                // rocof = 0
          t,                  // "start time" (irrelevant when rocof=0)
          simParams.timeStep, // small duration
          false);

      angleStepApplied = true;
      std::cout << "[HOOK][EMT] t=" << t << " infeed voltage angle step = "
                << simParams.infeedVoltageAngleStepDeg << " deg\n";
    }
  });
}

void simulateDP(const SimulationParameters &simParams,
                const PowerSystemParameters &psParams,
                const SystemTopology &systemTopologyPF,
                const PowerSystemEventType &psEvent) {
  String simName = "DP_simulation";
  Logger::setLogDir("logs/" + simName);
  auto logger = DataLogger::make(simName);

  // ---------------- Nodes ----------------
  auto node1 = DP::SimNode::make("node1", PhaseType::Single);
  auto node2 = DP::SimNode::make("node2", PhaseType::Single);
  auto node3 = DP::SimNode::make("node3", PhaseType::Single);
  auto node4 = DP::SimNode::make("node4", PhaseType::Single);
  auto node5 = DP::SimNode::make("node5", PhaseType::Single);
  auto node6 = DP::SimNode::make("node6", PhaseType::Single);
  auto node7 = DP::SimNode::make("node7", PhaseType::Single);

  auto node1s = DP::SimNode::make("node1_strong", PhaseType::Single);
  auto node1w = DP::SimNode::make("node1_weak", PhaseType::Single);

  // ---------------- Components ----------------
  auto infeedSource = DP::Ph1::NetworkInjection::make("infeed_source");
  infeedSource->connect({node1});

  const double kZ = std::max(1e-9, simParams.infeedImpedanceStepFactor);
  const double Rstrong = psParams.infeedResistance;
  const double Lstrong = psParams.infeedInductance;
  const double Rweak = psParams.infeedResistance * kZ;
  const double Lweak = psParams.infeedInductance * kZ;

  auto infeedSwStrong = DP::Ph1::Switch::make("infeed_sw_strong");
  infeedSwStrong->setParameters(SwitchConstants::openResistance,
                                SwitchConstants::closedResistance,
                                true); // initially CLOSED
  infeedSwStrong->connect({node1, node1s});

  auto infeedZStrong = DP::Ph1::PiLine::make("infeed_impedance_strong");
  infeedZStrong->setParameters(Rstrong, Lstrong, 0.0);
  infeedZStrong->connect({node1s, node4});

  auto infeedSwWeak = DP::Ph1::Switch::make("infeed_sw_weak");
  infeedSwWeak->setParameters(SwitchConstants::openResistance,
                              SwitchConstants::closedResistance,
                              false); // initially OPEN
  infeedSwWeak->connect({node1, node1w});

  auto infeedZWeak = DP::Ph1::PiLine::make("infeed_impedance_weak");
  infeedZWeak->setParameters(Rweak, Lweak, 0.0);
  infeedZWeak->connect({node1w, node4});

  logger->logAttribute(VariableNames::vInfeed,
                       node4->attribute(AttributeNames::v));
  logger->logAttribute(VariableNames::iInfeedStrong,
                       infeedZStrong->attribute(AttributeNames::i));
  logger->logAttribute(VariableNames::iInfeedWeak,
                       infeedZWeak->attribute(AttributeNames::i));
  logger->logAttribute(VariableNames::fInfeed,
                       infeedSource->attribute(AttributeNames::f));

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

  // load-bus fault as shunt switch node5 -> GND
  auto loadBusFault = DP::Ph1::Switch::make("load_bus_fault");
  loadBusFault->setParameters(SwitchConstants::openResistance,
                              simParams.faultResistance, false);
  loadBusFault->connect({node5, DP::SimNode::GND});
  logger->logAttribute(VariableNames::iFault,
                       loadBusFault->attribute(AttributeNames::i));

  auto load1 = DP::Ph1::Resistor::make("load1");
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

  // ---------------- Topology ----------------
  auto systemNodeList = SystemNodeList{node1, node1s, node1w, node2, node3,
                                       node4, node5,  node6,  node7};

  auto componentList = SystemComponentList{
      infeedSource, infeedSwStrong, infeedZStrong,  infeedSwWeak, infeedZWeak,
      line1,        line2,          circuitBreaker, loadBusFault, load1,
      load1Switch,  load2,          load2Switch};

  auto systemTopology =
      SystemTopology(psParams.frequency, systemNodeList, componentList);

  const bool doStartupRamp = simParams.enableStartupRampPQ;
  auto conv1 = createDPConverter(logger, psParams, systemTopology, node2, 1,
                                 doStartupRamp);
  auto conv2 = createDPConverter(logger, psParams, systemTopology, node3, 2,
                                 doStartupRamp);

  // Optional converter3 at same PCC as converter1 (node2)
  DPConverterHandle conv3{nullptr, 0.0, 0.0, 0.0, 0.0};
  const bool hasConv3 = simParams.enableConverter3;
  if (hasConv3) {
    conv3 = createDPConverter(logger, psParams, systemTopology, node2, 3, doStartupRamp);
  }

  // ---------------- Simulation ----------------
  systemTopology.initWithPowerflow(systemTopologyPF, Domain::DP);
  auto sim =
      setupSimulation(simName, simParams, systemTopology, logger, Domain::DP);

  // ---------------- Events ----------------
  switch (psEvent) {
  case PowerSystemEventType::LoadStep: {
    auto disconnectLoad1 =
        DPsim::SwitchEvent::make(simParams.eventTime, load1Switch, false);
    auto connectLoad2 =
        DPsim::SwitchEvent::make(simParams.eventTime, load2Switch, true);
    sim.addEvent(disconnectLoad1);
    sim.addEvent(connectLoad2);
    break;
  }
  case PowerSystemEventType::LoadBusFault: {
    const double tOn = simParams.eventTime;
    const double tOff =
        simParams.eventTime + std::max(0.0, simParams.faultDuration);
    auto faultOn = DPsim::SwitchEvent::make(tOn, loadBusFault, true);
    auto faultOff = DPsim::SwitchEvent::make(tOff, loadBusFault, false);
    sim.addEvent(faultOn);
    sim.addEvent(faultOff);
    break;
  }
  case PowerSystemEventType::InfeedSCRStep: {
    const double t = simParams.eventTime;
    auto openStrong = DPsim::SwitchEvent::make(t, infeedSwStrong, false);
    auto closeWeak = DPsim::SwitchEvent::make(t, infeedSwWeak, true);
    sim.addEvent(openStrong);
    sim.addEvent(closeWeak);
    break;
  }
  case PowerSystemEventType::InfeedFrequencyRamp: {
    infeedSource->setParameters(Complex(psParams.voltageLineToLine, 0), 0.0,
                                simParams.rocof, simParams.eventTime,
                                simParams.frequencyRampDuration, false);
    break;
  }
  case PowerSystemEventType::InfeedFrequencyStep: {
    infeedSource->setParameters(Complex(psParams.voltageLineToLine, 0), 0.0,
                                simParams.frequencyStepDelta /
                                    simParams.timeStep,
                                simParams.eventTime, simParams.timeStep, false);
    break;
  }
  case PowerSystemEventType::InfeedVoltageAngleStep:
  case PowerSystemEventType::Converter1PrefStep:
  case PowerSystemEventType::None:
  default:
    break;
  }

  // ---- Hook: startup ramp + optional Converter1 Pref step + optional infeed angle step ----
  const double holdT = std::max(0.0, simParams.startupPQZeroHoldTime);
  const double rampDur = std::max(0.0, simParams.startupRampDuration);
  const double rampEndT = holdT + rampDur;

  bool rampDone = !doStartupRamp;
  bool printedHold = false;
  bool printedDone = false;

  const double prefStepTime = simParams.eventTime;
  const double newPref = psParams.converter1P * simParams.prefStepFactor;
  bool prefStepApplied = false;

  // Infeed voltage angle step state
  const double angleStepTime = simParams.eventTime;
  const double deltaRad = simParams.infeedVoltageAngleStepDeg * M_PI / 180.0;
  bool angleStepApplied = false;

  runStepped(sim, [&](DPsim::Simulation &s) {
    const double t = s.time();

    if (doStartupRamp && !rampDone) {
      double alpha = 0.0;

      if (t < holdT) {
        alpha = 0.0;
        if (!printedHold && t >= (0.0 + 0.5 * simParams.timeStep)) {
          printedHold = true;
          std::cout << "[HOOK][DP] holding Pref/Qref at 0 for " << holdT
                    << " s\n";
        }
      } else if (rampDur <= 0.0) {
        alpha = 1.0;
      } else {
        alpha = clamp01((t - holdT) / rampDur);
      }

      conv1.conv->setParameters(conv1.sysOmega, conv1.sysVoltNom,
                                alpha * conv1.pFinal, alpha * conv1.qFinal);
      conv2.conv->setParameters(conv2.sysOmega, conv2.sysVoltNom,
                                alpha * conv2.pFinal, alpha * conv2.qFinal);
      if (hasConv3) {
        conv3.conv->setParameters(conv3.sysOmega, conv3.sysVoltNom,
                                  alpha * conv3.pFinal, alpha * conv3.qFinal);
      }

      if (!printedDone && t >= (rampEndT - 0.5 * simParams.timeStep)) {
        conv1.conv->setParameters(conv1.sysOmega, conv1.sysVoltNom,
                                  conv1.pFinal, conv1.qFinal);
        conv2.conv->setParameters(conv2.sysOmega, conv2.sysVoltNom,
                                  conv2.pFinal, conv2.qFinal);
        if (hasConv3) {
          conv3.conv->setParameters(conv3.sysOmega, conv3.sysVoltNom,
                                    conv3.pFinal, conv3.qFinal);
        }
        rampDone = true;
        printedDone = true;
        std::cout << "[HOOK][DP] t=" << t << " finished startup hold+ramp\n";
      }
    }

    if (psEvent == PowerSystemEventType::Converter1PrefStep &&
        !prefStepApplied && t >= (prefStepTime - 0.5 * simParams.timeStep)) {

      conv1.conv->setParameters(conv1.sysOmega, conv1.sysVoltNom, newPref,
                                conv1.qFinal);

      prefStepApplied = true;
      std::cout << "[HOOK][DP] t=" << t << " set Converter1 Pref=" << newPref
                << " W\n";
    }

    if (psEvent == PowerSystemEventType::InfeedVoltageAngleStep &&
        !angleStepApplied && t >= (angleStepTime - 0.5 * simParams.timeStep)) {

      const Complex Vref = std::polar(psParams.voltageLineToLine, deltaRad);

      infeedSource->setParameters(Vref,
                                  0.0, // 
                                  0.0, // rocof = 0
                                  t,   // start time (irrelevant when rocof=0)
                                  simParams.timeStep, false);

      angleStepApplied = true;
      std::cout << "[HOOK][DP] t=" << t << " infeed voltage angle step = "
                << simParams.infeedVoltageAngleStepDeg << " deg\n";
    }
  });
}

void simulateSP(const SimulationParameters &simParams,
                const PowerSystemParameters &psParams,
                const SystemTopology &systemTopologyPF,
                const PowerSystemEventType &psEvent) {
  String simName = "SP_simulation";
  Logger::setLogDir("logs/" + simName);
  auto logger = DataLogger::make(simName);

  // ---------------- Nodes ----------------
  auto node1 = SP::SimNode::make("node1", PhaseType::Single);
  auto node2 = SP::SimNode::make("node2", PhaseType::Single);
  auto node3 = SP::SimNode::make("node3", PhaseType::Single);
  auto node4 = SP::SimNode::make("node4", PhaseType::Single);
  auto node5 = SP::SimNode::make("node5", PhaseType::Single);
  auto node6 = SP::SimNode::make("node6", PhaseType::Single);
  auto node7 = SP::SimNode::make("node7", PhaseType::Single);

  auto node1s = SP::SimNode::make("node1_strong", PhaseType::Single);
  auto node1w = SP::SimNode::make("node1_weak", PhaseType::Single);

  // ---------------- Components ----------------
  auto infeedSource = SP::Ph1::NetworkInjection::make("infeed_source");
  infeedSource->connect({node1});

  const double kZ = std::max(1e-9, simParams.infeedImpedanceStepFactor);
  const double Rstrong = psParams.infeedResistance;
  const double Lstrong = psParams.infeedInductance;
  const double Rweak = psParams.infeedResistance * kZ;
  const double Lweak = psParams.infeedInductance * kZ;

  auto infeedSwStrong = SP::Ph1::Switch::make("infeed_sw_strong");
  infeedSwStrong->setParameters(SwitchConstants::openResistance,
                                SwitchConstants::closedResistance, true);
  infeedSwStrong->connect({node1, node1s});

  auto infeedZStrong = SP::Ph1::PiLine::make("infeed_impedance_strong");
  infeedZStrong->setParameters(Rstrong, Lstrong, 0.0);
  infeedZStrong->connect({node1s, node4});

  auto infeedSwWeak = SP::Ph1::Switch::make("infeed_sw_weak");
  infeedSwWeak->setParameters(SwitchConstants::openResistance,
                              SwitchConstants::closedResistance, false);
  infeedSwWeak->connect({node1, node1w});

  auto infeedZWeak = SP::Ph1::PiLine::make("infeed_impedance_weak");
  infeedZWeak->setParameters(Rweak, Lweak, 0.0);
  infeedZWeak->connect({node1w, node4});

  logger->logAttribute(VariableNames::vInfeed,
                       node4->attribute(AttributeNames::v));
  logger->logAttribute(VariableNames::iInfeedStrong,
                       infeedZStrong->attribute(AttributeNames::i));
  logger->logAttribute(VariableNames::iInfeedWeak,
                       infeedZWeak->attribute(AttributeNames::i));
  logger->logAttribute(VariableNames::fInfeed,
                       infeedSource->attribute(AttributeNames::f));

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

  // load-bus fault as shunt switch node5 -> GND
  auto loadBusFault = SP::Ph1::Switch::make("load_bus_fault");
  loadBusFault->setParameters(SwitchConstants::openResistance,
                              simParams.faultResistance, false);
  loadBusFault->connect({node5, SP::SimNode::GND});
  logger->logAttribute(VariableNames::iFault,
                       loadBusFault->attribute(AttributeNames::i));

  auto load1 = SP::Ph1::Resistor::make("load1");
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

  // ---------------- Topology ----------------
  auto systemNodeList = SystemNodeList{node1, node1s, node1w, node2, node3,
                                       node4, node5,  node6,  node7};

  auto componentList = SystemComponentList{
      infeedSource, infeedSwStrong, infeedZStrong,  infeedSwWeak, infeedZWeak,
      line1,        line2,          circuitBreaker, loadBusFault, load1,
      load1Switch,  load2,          load2Switch};

  auto systemTopology =
      SystemTopology(psParams.frequency, systemNodeList, componentList);

  const bool doStartupRamp = simParams.enableStartupRampPQ;
  auto conv1 = createSPConverter(logger, psParams, systemTopology, node2, 1,
                                 doStartupRamp);
  auto conv2 = createSPConverter(logger, psParams, systemTopology, node3, 2,
                                 doStartupRamp);

  // Optional converter3 at same PCC as converter1 (node2)
  SPConverterHandle conv3{nullptr, 0.0, 0.0, 0.0, 0.0};
  const bool hasConv3 = simParams.enableConverter3;
  if (hasConv3) {
    conv3 = createSPConverter(logger, psParams, systemTopology, node2, 3, doStartupRamp);
  }

  // ---------------- Simulation ----------------
  systemTopology.initWithPowerflow(systemTopologyPF, Domain::SP);
  auto sim =
      setupSimulation(simName, simParams, systemTopology, logger, Domain::SP);

  // ---------------- Events ----------------
  switch (psEvent) {
  case PowerSystemEventType::LoadStep: {
    auto disconnectLoad1 =
        DPsim::SwitchEvent::make(simParams.eventTime, load1Switch, false);
    auto connectLoad2 =
        DPsim::SwitchEvent::make(simParams.eventTime, load2Switch, true);
    sim.addEvent(disconnectLoad1);
    sim.addEvent(connectLoad2);
    break;
  }
  case PowerSystemEventType::LoadBusFault: {
    const double tOn = simParams.eventTime;
    const double tOff =
        simParams.eventTime + std::max(0.0, simParams.faultDuration);
    auto faultOn = DPsim::SwitchEvent::make(tOn, loadBusFault, true);
    auto faultOff = DPsim::SwitchEvent::make(tOff, loadBusFault, false);
    sim.addEvent(faultOn);
    sim.addEvent(faultOff);
    break;
  }
  case PowerSystemEventType::InfeedSCRStep: {
    const double t = simParams.eventTime;
    auto openStrong = DPsim::SwitchEvent::make(t, infeedSwStrong, false);
    auto closeWeak = DPsim::SwitchEvent::make(t, infeedSwWeak, true);
    sim.addEvent(openStrong);
    sim.addEvent(closeWeak);
    break;
  }
  case PowerSystemEventType::InfeedFrequencyRamp: {
    infeedSource->setParameters(Complex(psParams.voltageLineToLine, 0), 0.0,
                                simParams.rocof, simParams.eventTime,
                                simParams.frequencyRampDuration, false);
    break;
  }
  case PowerSystemEventType::InfeedFrequencyStep: {
    infeedSource->setParameters(Complex(psParams.voltageLineToLine, 0), 0.0,
                                simParams.frequencyStepDelta /
                                    simParams.timeStep,
                                simParams.eventTime, simParams.timeStep, false);
    break;
  }
  case PowerSystemEventType::InfeedVoltageAngleStep:
  case PowerSystemEventType::Converter1PrefStep:
  case PowerSystemEventType::None:
  default:
    break;
  }

  // ---- Hook: startup ramp + optional Converter1 Pref step + optional infeed angle step ----
  const double holdT = std::max(0.0, simParams.startupPQZeroHoldTime);
  const double rampDur = std::max(0.0, simParams.startupRampDuration);
  const double rampEndT = holdT + rampDur;

  bool rampDone = !doStartupRamp;
  bool printedHold = false;
  bool printedDone = false;

  const double prefStepTime = simParams.eventTime;
  const double newPref = psParams.converter1P * simParams.prefStepFactor;
  bool prefStepApplied = false;

  // Infeed voltage angle step state
  const double angleStepTime = simParams.eventTime;
  const double deltaRad = simParams.infeedVoltageAngleStepDeg * M_PI / 180.0;
  bool angleStepApplied = false;

  runStepped(sim, [&](DPsim::Simulation &s) {
    const double t = s.time();

    if (doStartupRamp && !rampDone) {
      double alpha = 0.0;

      if (t < holdT) {
        alpha = 0.0;
        if (!printedHold && t >= (0.0 + 0.5 * simParams.timeStep)) {
          printedHold = true;
          std::cout << "[HOOK][SP] holding Pref/Qref at 0 for " << holdT
                    << " s\n";
        }
      } else if (rampDur <= 0.0) {
        alpha = 1.0;
      } else {
        alpha = clamp01((t - holdT) / rampDur);
      }

      conv1.conv->setParameters(conv1.sysOmega, conv1.sysVoltNom,
                                alpha * conv1.pFinal, alpha * conv1.qFinal);
      conv2.conv->setParameters(conv2.sysOmega, conv2.sysVoltNom,
                                alpha * conv2.pFinal, alpha * conv2.qFinal);
      if (hasConv3) {
        conv3.conv->setParameters(conv3.sysOmega, conv3.sysVoltNom,
                                  alpha * conv3.pFinal, alpha * conv3.qFinal);
      }

      if (!printedDone && t >= (rampEndT - 0.5 * simParams.timeStep)) {
        conv1.conv->setParameters(conv1.sysOmega, conv1.sysVoltNom,
                                  conv1.pFinal, conv1.qFinal);
        conv2.conv->setParameters(conv2.sysOmega, conv2.sysVoltNom,
                                  conv2.pFinal, conv2.qFinal);
        if (hasConv3) {
          conv3.conv->setParameters(conv3.sysOmega, conv3.sysVoltNom,
                                    conv3.pFinal, conv3.qFinal);
        }
        rampDone = true;
        printedDone = true;
        std::cout << "[HOOK][SP] t=" << t << " finished startup hold+ramp\n";
      }
    }

    if (psEvent == PowerSystemEventType::Converter1PrefStep &&
        !prefStepApplied && t >= (prefStepTime - 0.5 * simParams.timeStep)) {

      conv1.conv->setParameters(conv1.sysOmega, conv1.sysVoltNom, newPref,
                                conv1.qFinal);

      prefStepApplied = true;
      std::cout << "[HOOK][SP] t=" << t << " set Converter1 Pref=" << newPref
                << " W\n";
    }

    if (psEvent == PowerSystemEventType::InfeedVoltageAngleStep &&
        !angleStepApplied && t >= (angleStepTime - 0.5 * simParams.timeStep)) {

      const Complex Vref = std::polar(psParams.voltageLineToLine, deltaRad);
      infeedSource->setParameters(Vref,
                                  0.0, // 
                                  0.0, // rocof = 0
                                  t,   // start time (irrelevant when rocof=0)
                                  simParams.timeStep, false);

      angleStepApplied = true;
      std::cout << "[HOOK][SP] t=" << t << " infeed voltage angle step = "
                << simParams.infeedVoltageAngleStepDeg << " deg\n";
    }
  });
}

// --------- PF calculation (updated to include SCR-step branches + fault branch) ---------

SystemTopology calculatePF(const SimulationParameters &simParams,
                           const PowerSystemParameters &psParams) {
  String simName = "PF_calculation";
  Logger::setLogDir("logs/" + simName);
  auto logger = DataLogger::make(simName);

  // ---------------- Nodes (SP PF) ----------------
  auto node1 = SP::SimNode::make("node1", PhaseType::Single);
  auto node2 = SP::SimNode::make("node2", PhaseType::Single);
  auto node3 = SP::SimNode::make("node3", PhaseType::Single);
  auto node4 = SP::SimNode::make("node4", PhaseType::Single);
  auto node5 = SP::SimNode::make("node5", PhaseType::Single);
  auto node6 = SP::SimNode::make("node6", PhaseType::Single);
  auto node7 = SP::SimNode::make("node7", PhaseType::Single);

  // extra nodes for SCR-step branches
  auto node1s = SP::SimNode::make("node1_strong", PhaseType::Single);
  auto node1w = SP::SimNode::make("node1_weak", PhaseType::Single);

  // ---------------- Components ----------------
  auto infeedSource =
      SP::Ph1::NetworkInjection::make("infeed_source", Logger::Level::debug);
  infeedSource->setParameters(psParams.voltageLineToLine);
  infeedSource->setBaseVoltage(psParams.voltageLineToLine);
  infeedSource->modifyPowerFlowBusType(PowerflowBusType::VD);
  infeedSource->connect({node1});

  const double kZ = std::max(1e-9, simParams.infeedImpedanceStepFactor);
  const double Rstrong = psParams.infeedResistance;
  const double Lstrong = psParams.infeedInductance;
  const double Rweak = psParams.infeedResistance * kZ;
  const double Lweak = psParams.infeedInductance * kZ;

  // In PF we approximate "switches" as purely resistive PiLines.
  // Start in STRONG configuration: strong-switch closed, weak-switch open.
  auto infeedSwStrongPF =
      SP::Ph1::PiLine::make("infeed_sw_strong", Logger::Level::debug);
  infeedSwStrongPF->setParameters(SwitchConstants::closedResistance, 0.0);
  infeedSwStrongPF->setBaseVoltage(psParams.voltageLineToLine);
  infeedSwStrongPF->connect({node1, node1s});

  auto infeedZStrongPF =
      SP::Ph1::PiLine::make("infeed_impedance_strong", Logger::Level::debug);
  infeedZStrongPF->setParameters(Rstrong, Lstrong, 0.0);
  infeedZStrongPF->setBaseVoltage(psParams.voltageLineToLine);
  infeedZStrongPF->connect({node1s, node4});

  auto infeedSwWeakPF =
      SP::Ph1::PiLine::make("infeed_sw_weak", Logger::Level::debug);
  infeedSwWeakPF->setParameters(SwitchConstants::openResistance, 0.0);
  infeedSwWeakPF->setBaseVoltage(psParams.voltageLineToLine);
  infeedSwWeakPF->connect({node1, node1w});

  auto infeedZWeakPF =
      SP::Ph1::PiLine::make("infeed_impedance_weak", Logger::Level::debug);
  infeedZWeakPF->setParameters(Rweak, Lweak, 0.0);
  infeedZWeakPF->setBaseVoltage(psParams.voltageLineToLine);
  infeedZWeakPF->connect({node1w, node4});

  auto line1 = SP::Ph1::PiLine::make("line1", Logger::Level::debug);
  line1->setParameters(psParams.line1Resistance, psParams.line1Inductance,
                       psParams.line1Capacitance);
  line1->setBaseVoltage(psParams.voltageLineToLine);
  line1->connect({node2, node4});

  auto line2 = SP::Ph1::PiLine::make("line2", Logger::Level::debug);
  line2->setParameters(psParams.line2Resistance, psParams.line2Inductance,
                       psParams.line2Capacitance);
  line2->setBaseVoltage(psParams.voltageLineToLine);
  line2->connect({node3, node4});

  auto circuitBreaker =
      SP::Ph1::PiLine::make("circuit_breaker", Logger::Level::debug);
  circuitBreaker->setParameters(SwitchConstants::closedResistance, 0.0);
  circuitBreaker->setBaseVoltage(psParams.voltageLineToLine);
  circuitBreaker->connect({node4, node5});

  // include fault branch in PF as "open" (very large shunt resistance)
  auto loadBusFaultPF =
      SP::Ph1::PiLine::make("load_bus_fault", Logger::Level::debug);
  loadBusFaultPF->setParameters(SwitchConstants::openResistance, 0.0);
  loadBusFaultPF->setBaseVoltage(psParams.voltageLineToLine);
  loadBusFaultPF->connect({node5, SP::SimNode::GND});

  auto load1 = SP::Ph1::PiLine::make("load1", Logger::Level::debug);
  load1->setParameters(psParams.loadResistance1, 0.0);
  load1->setBaseVoltage(psParams.voltageLineToLine);
  load1->connect({node6, SP::SimNode::GND});

  auto load1Switch =
      SP::Ph1::PiLine::make("load1_switch", Logger::Level::debug);
  load1Switch->setParameters(SwitchConstants::closedResistance, 0.0);
  load1Switch->setBaseVoltage(psParams.voltageLineToLine);
  load1Switch->connect({node5, node6});

  auto load2 = SP::Ph1::PiLine::make("load2", Logger::Level::debug);
  load2->setParameters(psParams.loadResistance2, 0.0);
  load2->setBaseVoltage(psParams.voltageLineToLine);
  load2->connect({node7, SP::SimNode::GND});

  auto load2Switch =
      SP::Ph1::PiLine::make("load2_switch", Logger::Level::debug);
  load2Switch->setParameters(SwitchConstants::openResistance, 0.0);
  load2Switch->setBaseVoltage(psParams.voltageLineToLine);
  load2Switch->connect({node5, node7});

  auto converter1 = SP::Ph1::Load::make("Converter1", Logger::Level::debug);
  converter1->setParameters(-psParams.converter1P, -psParams.converter1Q,
                            psParams.voltageLineToLine);
  converter1->modifyPowerFlowBusType(PowerflowBusType::PQ);
  converter1->connect({node2});

  // Optional converter3 in PF (same as converter1, same PCC node2)
  std::shared_ptr<SP::Ph1::Load> converter3 = nullptr;
  if (simParams.enableConverter3) {
    converter3 = SP::Ph1::Load::make("Converter3", Logger::Level::debug);
    converter3->setParameters(-psParams.converter1P, -psParams.converter1Q,
                              psParams.voltageLineToLine);
    converter3->modifyPowerFlowBusType(PowerflowBusType::PQ);
    converter3->connect({node2});
  }

  auto converter2 = SP::Ph1::Load::make("Converter2", Logger::Level::debug);
  converter2->setParameters(-psParams.converter2P, -psParams.converter2Q,
                            psParams.voltageLineToLine);
  converter2->modifyPowerFlowBusType(PowerflowBusType::PQ);
  converter2->connect({node3});

  // ---------------- Topology ----------------
  auto systemNodeList = SystemNodeList{node1, node1s, node1w, node2, node3,
                                       node4, node5,  node6,  node7};

  // Build component list (so Converter3 can be optional without inserting nullptrs)
  SystemComponentList componentList;
  componentList.push_back(infeedSource);
  componentList.push_back(infeedSwStrongPF);
  componentList.push_back(infeedZStrongPF);
  componentList.push_back(infeedSwWeakPF);
  componentList.push_back(infeedZWeakPF);

  componentList.push_back(converter1);
  if (converter3) {
    componentList.push_back(converter3);
  }
  componentList.push_back(line1);

  componentList.push_back(converter2);
  componentList.push_back(line2);

  componentList.push_back(circuitBreaker);
  componentList.push_back(loadBusFaultPF);
  componentList.push_back(load1);
  componentList.push_back(load1Switch);
  componentList.push_back(load2);
  componentList.push_back(load2Switch);

  auto systemTopology =
      SystemTopology(psParams.frequency, systemNodeList, componentList);

  // logging (optional)
  logger->logAttribute(VariableNames::vInfeed,
                       node1->attribute(AttributeNames::v));
  logger->logAttribute("vConverter1", node2->attribute(AttributeNames::v));
  if (simParams.enableConverter3) {
    logger->logAttribute("vConverter3", node2->attribute(AttributeNames::v));
  }
  logger->logAttribute("vLoadBus", node5->attribute(AttributeNames::v));

  // ---------------- Simulation ----------------
  Simulation sim(simName, Logger::Level::debug);
  sim.setSystem(systemTopology);
  sim.setTimeStep(simParams.finalTime);
  sim.setFinalTime(2 * simParams.finalTime);
  sim.setDomain(Domain::SP);
  sim.setSolverType(Solver::Type::NRP);
  sim.setSolverAndComponentBehaviour(Solver::Behaviour::Initialization);
  sim.doInitFromNodesAndTerminals(false);
  sim.addLogger(logger);
  sim.run();

  return systemTopology;
}

} // namespace HVDCWise

int main() {
  HVDCWise::SimulationParameters simParams;
  HVDCWise::PowerSystemInputParameters psInputParams;

  HVDCWise::PowerSystemParameters psParams =
      HVDCWise::calculatePowerSystemParameters(psInputParams);

  // Choose one:
  // auto psEvent = HVDCWise::PowerSystemEventType::LoadBusFault;
  auto psEvent = HVDCWise::PowerSystemEventType::Converter1PrefStep;
  // auto psEvent = HVDCWise::PowerSystemEventType::InfeedSCRStep;
  // auto psEvent = HVDCWise::PowerSystemEventType::InfeedVoltageAngleStep;

  auto systemTopologyPF = HVDCWise::calculatePF(simParams, psParams);

  HVDCWise::simulateEMT(simParams, psParams, systemTopologyPF, psEvent);
  HVDCWise::simulateDP(simParams, psParams, systemTopologyPF, psEvent);
  HVDCWise::simulateSP(simParams, psParams, systemTopologyPF, psEvent);

  return 0;
}