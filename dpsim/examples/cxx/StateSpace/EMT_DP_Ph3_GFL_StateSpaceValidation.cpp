// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include <DPsim.h>

#include <dpsim-models/DP/DP_Ph3_AvVoltSourceInverterStateSpace.h>
#include <dpsim-models/DP/DP_Ph3_SSN_GFL_Split.h>
#include <dpsim-models/EMT/EMT_Ph3_SSN_GFL.h>
#include <dpsim-models/EMT/EMT_Ph3_SSN_GFL_Split.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace CPS;
using namespace DPsim;

namespace {

constexpr UInt ControllerStateCount = 8;
constexpr UInt PhysicalStateCount = 9;
constexpr UInt NoDelayStateCount = ControllerStateCount + PhysicalStateCount;
constexpr UInt DelayStateCount = 3;
constexpr UInt SplitReferenceStateCount = NoDelayStateCount + DelayStateCount;

const Real K23 = std::sqrt(2.0 / 3.0);

enum class ModelKind { EmtVariable, EmtSplit, DpVariable, DpSplit };

struct Parameters {
  Real frequency = 50.0;
  Real omega = 2.0 * PI * frequency;
  Real voltageRmsLineToLine = 400.0;
  Real gridResistance = 0.3;
  Real gridInductance = 1.0e-3;
  Real filterInductance = 2.0e-3;
  Real filterCapacitance = 10.0e-6;
  Real filterResistance = 0.2;
  Real couplingResistance = 0.2;
  Real kpPll = 0.25;
  Real kiPll = 0.2;
  Real powerCutoff = omega;
  Real activePowerReference = 10000.0;
  Real reactivePowerReference = 5000.0;
  Real kpPowerBase = 0.05;
  Real kiPowerBase = 0.2;
  Real kpCurrentBase = 0.25;
  Real kiCurrentBase = 1.0;
  Real stableGainScale = 3.7;
  Real unstableGainScale = 3.8;
  Real perturbationTime = 1.0 / frequency;
  Real perturbationDuration = 0.5e-3;
  Real gridVoltagePulseRelative = 1.0e-3;
  Real gainCaseFinalTime = perturbationTime + 0.2;
  Real largeStepFinalTime = perturbationTime + 0.5;
};

struct OperatingPoint {
  Complex source;
  Complex intermediate;
  Complex terminal;
  Complex capacitorVoltage;
  Complex gridCurrent;
  Complex filterCurrent;
  Complex bridgeVoltage;
};

struct ModalResult {
  VectorComp discreteEigenvalues;
  VectorComp continuousEigenvalues;
  MatrixComp participationFactors;
  std::vector<String> stateNames;
  Real extractionTime = 0.0;
  UInt stateCount = 0;
};

struct ReferenceResult {
  Matrix noDelayA;
  Matrix noDelayAd;
  Matrix splitAd;
};

struct EigenvalueRecord {
  String study;
  String model;
  String method;
  Real gainScale;
  Real timeStep;
  VectorComp discrete;
  VectorComp continuous;
};

struct SummaryRecord {
  String study;
  String model;
  Real gainScale;
  Real timeStep;
  UInt stateCount;
  Real extractionTime;
  Real maxRealPart;
  Real spectralRadius;
  Real discreteReferenceError;
  Real referenceError;
  Real noDelayDeviation;
};

struct ParticipationRecord {
  String model;
  UInt mode;
  UInt state;
  String stateName;
  Complex value;
};

struct TimeDomainRecord {
  String model;
  String stabilityCase;
  Real gainScale;
  Real timeStep;
  Real finalTime;
  Real perturbationTime;
  Real perturbationDuration;
  Real perturbationRelative;
  Real prePerturbationActivePowerError;
  Real prePerturbationReactivePowerError;
  String logPath;
};

struct TimeDomainResult {
  String logPath;
  Real prePerturbationActivePowerError;
  Real prePerturbationReactivePowerError;
};

String modelName(ModelKind model) {
  switch (model) {
  case ModelKind::EmtVariable:
    return "EMT variable";
  case ModelKind::EmtSplit:
    return "EMT split";
  case ModelKind::DpVariable:
    return "DP Ph3 variable";
  case ModelKind::DpSplit:
    return "DP Ph3 split";
  }
  throw std::logic_error("Unsupported model kind.");
}

String fileToken(ModelKind model) {
  switch (model) {
  case ModelKind::EmtVariable:
    return "emt_variable";
  case ModelKind::EmtSplit:
    return "emt_split";
  case ModelKind::DpVariable:
    return "dp_ph3_variable";
  case ModelKind::DpSplit:
    return "dp_ph3_split";
  }
  throw std::logic_error("Unsupported model kind.");
}

Bool isEmt(ModelKind model) {
  return model == ModelKind::EmtVariable || model == ModelKind::EmtSplit;
}

Bool isSplit(ModelKind model) {
  return model == ModelKind::EmtSplit || model == ModelKind::DpSplit;
}

String timeStepToken(Real timeStep) {
  const Real microseconds = 1e6 * timeStep;
  if (microseconds < 1000.0)
    return std::to_string(static_cast<Int>(std::llround(microseconds))) + "us";
  return std::to_string(static_cast<Int>(std::llround(1e3 * timeStep))) +
         "ms";
}

String gainToken(Real gainScale) {
  String result = std::to_string(gainScale);
  while (!result.empty() && result.back() == '0')
    result.pop_back();
  if (!result.empty() && result.back() == '.')
    result.pop_back();
  std::replace(result.begin(), result.end(), '.', '_');
  return result;
}

MatrixComp balancedEnvelope(const Complex &dqValue) {
  const Complex phaseA = K23 * dqValue;
  MatrixComp result(3, 1);
  result << phaseA, phaseA * SHIFT_TO_PHASE_B, phaseA * SHIFT_TO_PHASE_C;
  return result;
}

Matrix trapezoidalStateMatrix(const Matrix &a, Real timeStep) {
  const Matrix identity = Matrix::Identity(a.rows(), a.cols());
  return (identity - 0.5 * timeStep * a)
      .partialPivLu()
      .solve(identity + 0.5 * timeStep * a);
}

void trapezoidalMatrices(const Matrix &a, const Matrix &b, Real timeStep,
                         Matrix &ad, Matrix &bd) {
  const Matrix identity = Matrix::Identity(a.rows(), a.cols());
  const Matrix lhs = identity - 0.5 * timeStep * a;
  ad = lhs.partialPivLu().solve(identity + 0.5 * timeStep * a);
  bd = lhs.partialPivLu().solve(0.5 * timeStep * b);
}

VectorComp eigenvalues(const Matrix &matrix) {
  if (!matrix.allFinite())
    throw std::runtime_error(
        "Eigenvalue computation received a non-finite matrix.");
  Eigen::EigenSolver<Matrix> solver(matrix);
  if (solver.info() != Eigen::Success)
    throw std::runtime_error("Eigenvalue computation failed.");
  return solver.eigenvalues();
}

VectorComp bilinearContinuousEigenvalues(const VectorComp &values,
                                         Real timeStep) {
  VectorComp result(values.rows());
  for (Eigen::Index idx = 0; idx < values.rows(); ++idx)
    result(idx) = (2.0 / timeStep) * (values(idx) - Complex(1.0, 0.0)) /
                  (values(idx) + Complex(1.0, 0.0));
  return result;
}

VectorComp logarithmicContinuousEigenvalues(const VectorComp &values,
                                            Real timeStep) {
  VectorComp result(values.rows());
  for (Eigen::Index idx = 0; idx < values.rows(); ++idx)
    result(idx) = std::log(values(idx)) / timeStep;
  return result;
}

Real directedEigenvalueDistance(const VectorComp &reference,
                                const VectorComp &candidate) {
  Real maximum = 0.0;
  for (Eigen::Index refIdx = 0; refIdx < reference.rows(); ++refIdx) {
    Real nearest = std::numeric_limits<Real>::max();
    for (Eigen::Index valueIdx = 0; valueIdx < candidate.rows(); ++valueIdx)
      nearest = std::min(
          nearest, std::abs(reference(refIdx) - candidate(valueIdx)));
    maximum = std::max(maximum, nearest);
  }
  return maximum;
}

Bool isFinite(const Complex &value) {
  return std::isfinite(value.real()) && std::isfinite(value.imag());
}

Real directedFiniteEigenvalueDistance(const VectorComp &reference,
                                      const VectorComp &candidate) {
  Real maximum = 0.0;
  Bool hasFiniteReference = false;
  for (Eigen::Index refIdx = 0; refIdx < reference.rows(); ++refIdx) {
    if (!isFinite(reference(refIdx)))
      continue;
    hasFiniteReference = true;
    Real nearest = std::numeric_limits<Real>::infinity();
    for (Eigen::Index valueIdx = 0; valueIdx < candidate.rows(); ++valueIdx) {
      if (!isFinite(candidate(valueIdx)))
        continue;
      nearest = std::min(
          nearest, std::abs(reference(refIdx) - candidate(valueIdx)));
    }
    maximum = std::max(maximum, nearest);
  }
  return hasFiniteReference ? maximum
                            : std::numeric_limits<Real>::quiet_NaN();
}

Real maximumRealPart(const VectorComp &values) {
  Real result = -std::numeric_limits<Real>::infinity();
  for (Eigen::Index idx = 0; idx < values.rows(); ++idx) {
    if (std::isfinite(values(idx).real()))
      result = std::max(result, values(idx).real());
  }
  return result;
}

Real spectralRadius(const VectorComp &values) {
  Real result = 0.0;
  for (Eigen::Index idx = 0; idx < values.rows(); ++idx) {
    if (isFinite(values(idx)))
      result = std::max(result, std::abs(values(idx)));
  }
  return result;
}

void appendParticipationFromMatrix(
    const String &model, const Matrix &matrix,
    const std::vector<String> &stateNames,
    std::vector<ParticipationRecord> &records) {
  Eigen::EigenSolver<Matrix> solver(matrix, true);
  if (solver.info() != Eigen::Success)
    throw std::runtime_error("Reference modal decomposition failed.");
  const MatrixComp right = solver.eigenvectors();
  const MatrixComp left = right.fullPivLu().inverse();
  const MatrixComp participation =
      Math::elementwiseProduct(right, left.transpose());
  for (UInt mode = 0; mode < static_cast<UInt>(matrix.rows()); ++mode) {
    for (UInt state = 0; state < static_cast<UInt>(matrix.rows()); ++state) {
      records.push_back({model, mode, state, stateNames[state],
                         participation(state, mode)});
    }
  }
}

Matrix numericalJacobian(const std::function<Matrix(const Matrix &)> &function,
                         const Matrix &operatingPoint) {
  const Matrix value = function(operatingPoint);
  Matrix jacobian(value.rows(), operatingPoint.rows());
  for (Eigen::Index column = 0; column < operatingPoint.rows(); ++column) {
    const Real step =
        1e-6 * std::max(1.0, std::abs(operatingPoint(column, 0)));
    Matrix plus = operatingPoint;
    Matrix minus = operatingPoint;
    plus(column, 0) += step;
    minus(column, 0) -= step;
    jacobian.col(column) = (function(plus) - function(minus)) / (2.0 * step);
  }
  return jacobian;
}

Matrix numericalInputJacobian(
    const std::function<Matrix(const Matrix &, const Matrix &)> &function,
    const Matrix &state, const Matrix &input) {
  const Matrix value = function(state, input);
  Matrix jacobian(value.rows(), input.rows());
  for (Eigen::Index column = 0; column < input.rows(); ++column) {
    const Real step = 1e-6 * std::max(1.0, std::abs(input(column, 0)));
    Matrix plus = input;
    Matrix minus = input;
    plus(column, 0) += step;
    minus(column, 0) -= step;
    jacobian.col(column) =
        (function(state, plus) - function(state, minus)) / (2.0 * step);
  }
  return jacobian;
}

void writeParameters(const std::filesystem::path &path, const Parameters &p) {
  std::ofstream stream(path);
  stream << std::setprecision(std::numeric_limits<Real>::max_digits10)
         << "parameter,value,unit\n"
         << "frequency," << p.frequency << ",Hz\n"
         << "voltage_rms_ll," << p.voltageRmsLineToLine << ",V\n"
         << "grid_resistance," << p.gridResistance << ",ohm\n"
         << "grid_inductance," << p.gridInductance << ",H\n"
         << "filter_inductance," << p.filterInductance << ",H\n"
         << "filter_capacitance," << p.filterCapacitance << ",F\n"
         << "filter_resistance," << p.filterResistance << ",ohm\n"
         << "coupling_resistance," << p.couplingResistance << ",ohm\n"
         << "kp_pll," << p.kpPll << ",-\n"
         << "ki_pll," << p.kiPll << ",-\n"
         << "power_cutoff," << p.powerCutoff << ",rad_per_s\n"
         << "active_power_reference," << p.activePowerReference << ",W\n"
         << "reactive_power_reference," << p.reactivePowerReference
         << ",var\n"
         << "kp_power_base," << p.kpPowerBase << ",-\n"
         << "ki_power_base," << p.kiPowerBase << ",-\n"
         << "kp_current_base," << p.kpCurrentBase << ",-\n"
         << "ki_current_base," << p.kiCurrentBase << ",-\n"
         << "stable_gain_scale," << p.stableGainScale << ",-\n"
         << "unstable_gain_scale," << p.unstableGainScale << ",-\n";
  const Real operatingApparentPower =
      std::hypot(p.activePowerReference, p.reactivePowerReference);
  const Real gridImpedance =
      std::hypot(p.gridResistance, p.omega * p.gridInductance);
  const Real resonanceFrequency =
      std::sqrt((p.filterInductance + p.gridInductance) /
                (p.filterInductance * p.gridInductance *
                 p.filterCapacitance)) /
      (2.0 * PI);
  stream << "operating_apparent_power," << operatingApparentPower << ",VA\n"
         << "grid_x_over_r,"
         << p.omega * p.gridInductance / p.gridResistance << ",-\n"
         << "operating_point_scr,"
         << p.voltageRmsLineToLine * p.voltageRmsLineToLine /
                (gridImpedance * operatingApparentPower)
         << ",-\n"
         << "undamped_lcl_resonance_frequency," << resonanceFrequency
         << ",Hz\n"
         << "perturbation_time," << p.perturbationTime << ",s\n"
         << "perturbation_duration," << p.perturbationDuration << ",s\n"
         << "grid_voltage_pulse_relative," << p.gridVoltagePulseRelative
         << ",-\n"
         << "gain_case_final_time," << p.gainCaseFinalTime << ",s\n"
         << "large_step_final_time," << p.largeStepFinalTime << ",s\n";
}

} // namespace

class EMTDPPh3GFLStateSpaceValidation {
public:
  explicit EMTDPPh3GFLStateSpaceValidation(Bool runTimeDomain)
      : mRunTimeDomain(runTimeDomain),
        mOutputDirectory(std::filesystem::path("logs") /
                         "EMT_DP_Ph3_GFL_StateSpaceValidation") {}

  void run();

private:
  struct SystemHandles {
    SystemTopology system;
    Attribute<MatrixComp>::Ptr sourceVoltageReference;
    Attribute<Real>::Ptr activePower;
    Attribute<Real>::Ptr reactivePower;
  };

  OperatingPoint operatingPoint() const;
  SystemTopology runPowerFlow(const String &name) const;
  SystemHandles buildEmtSystem(ModelKind model, const SystemTopology &powerFlow,
                               Real gainScale,
                               const std::shared_ptr<DataLogger> &logger) const;
  SystemHandles buildDpSystem(ModelKind model, const OperatingPoint &op,
                              Real gainScale,
                              const std::shared_ptr<DataLogger> &logger) const;
  SystemHandles buildSystem(ModelKind model, const SystemTopology &powerFlow,
                            const OperatingPoint &op, Real gainScale,
                            const std::shared_ptr<DataLogger> &logger) const;

  ModalResult extractOneStep(ModelKind model, const SystemTopology &powerFlow,
                             const OperatingPoint &op, Real gainScale,
                             Real timeStep, const String &study) const;
  VectorComp calculateMonodromy(ModelKind model,
                                const SystemTopology &powerFlow,
                                const OperatingPoint &op, Real gainScale,
                                Real timeStep) const;
  TimeDomainResult
  runTimeDomainCase(ModelKind model, const SystemTopology &powerFlow,
                    const OperatingPoint &op, Real gainScale, Real timeStep,
                    Real finalTime, const String &stabilityCase) const;

  Matrix controllerDerivative(const Matrix &controllerState,
                              const Matrix &measurement,
                              Real gainScale) const;
  Matrix controllerOutput(const Matrix &controllerState,
                          const Matrix &measurement, Real gainScale) const;
  Matrix physicalDerivative(const Matrix &physicalState,
                            const Matrix &bridgeVoltage) const;
  ReferenceResult buildReferences(const OperatingPoint &op, Real gainScale,
                                  Real timeStep) const;

  void writeEigenvalues(const std::vector<EigenvalueRecord> &records) const;
  void writeSummary(const std::vector<SummaryRecord> &records) const;
  void writeParticipation(
      const std::vector<ParticipationRecord> &records) const;
  void writeTimeDomainManifest(
      const std::vector<TimeDomainRecord> &records) const;

  Parameters mParameters;
  Bool mRunTimeDomain;
  std::filesystem::path mOutputDirectory;
};

OperatingPoint EMTDPPh3GFLStateSpaceValidation::operatingPoint() const {
  // The analytical reference is expressed in the power-invariant global dq0
  // frame. Its d-axis voltage magnitude therefore equals the line-to-line RMS
  // voltage of the balanced system.
  const Complex source(mParameters.voltageRmsLineToLine, 0.0);
  const Complex power(mParameters.activePowerReference,
                      mParameters.reactivePowerReference);
  const Complex impedance(
      mParameters.couplingResistance + mParameters.gridResistance,
      mParameters.omega * mParameters.gridInductance);

  Complex capacitorVoltage = source;
  Complex gridCurrent(0.0, 0.0);
  for (UInt iteration = 0; iteration < 100; ++iteration) {
    gridCurrent = std::conj(power / capacitorVoltage);
    const Complex next = source + impedance * gridCurrent;
    if (std::abs(next - capacitorVoltage) < 1e-13) {
      capacitorVoltage = next;
      break;
    }
    capacitorVoltage = next;
  }

  gridCurrent = std::conj(power / capacitorVoltage);
  const Complex terminal =
      capacitorVoltage - mParameters.couplingResistance * gridCurrent;
  const Complex intermediate =
      source + Complex(0.0, mParameters.omega * mParameters.gridInductance) *
                   gridCurrent;
  const Complex filterCurrent =
      gridCurrent + Complex(0.0, mParameters.omega *
                                     mParameters.filterCapacitance) *
                        capacitorVoltage;
  const Complex bridgeVoltage =
      capacitorVoltage +
      Complex(mParameters.filterResistance,
              mParameters.omega * mParameters.filterInductance) *
          filterCurrent;

  return {source, intermediate, terminal, capacitorVoltage, gridCurrent,
          filterCurrent, bridgeVoltage};
}

SystemTopology
EMTDPPh3GFLStateSpaceValidation::runPowerFlow(const String &name) const {
  auto grid = SimNode<Complex>::make("nGrid", PhaseType::Single);
  auto middle = SimNode<Complex>::make("nMiddle", PhaseType::Single);
  auto pcc = SimNode<Complex>::make("nPcc", PhaseType::Single);

  auto source = SP::Ph1::NetworkInjection::make("GridSourcePF");
  source->setParameters(mParameters.voltageRmsLineToLine);
  source->setBaseVoltage(mParameters.voltageRmsLineToLine);
  source->modifyPowerFlowBusType(PowerflowBusType::VD);

  auto resistance = SP::Ph1::PiLine::make("GridResistancePF");
  resistance->setParameters(mParameters.gridResistance, 0.0, 0.0, 0.0);
  resistance->setBaseVoltage(mParameters.voltageRmsLineToLine);
  auto inductance = SP::Ph1::PiLine::make("GridInductancePF");
  inductance->setParameters(0.0, mParameters.gridInductance, 0.0, 0.0);
  inductance->setBaseVoltage(mParameters.voltageRmsLineToLine);

  // P/Q are controlled at the filter-capacitor node, which is separated from
  // the external PCC by Rc inside the inverter model. Initialize the external
  // power-flow network with the corresponding PCC power, including the Rc
  // loss, rather than applying the capacitor-side reference directly at PCC.
  const OperatingPoint op = operatingPoint();
  const Complex pccPower = op.terminal * std::conj(op.gridCurrent);
  auto inverter = SP::Ph1::Load::make("InverterPF");
  inverter->setParameters(-pccPower.real(), -pccPower.imag(),
                          mParameters.voltageRmsLineToLine);
  inverter->modifyPowerFlowBusType(PowerflowBusType::PQ);

  source->connect({grid});
  inductance->connect({grid, middle});
  resistance->connect({middle, pcc});
  inverter->connect({pcc});

  SystemTopology system(
      mParameters.frequency, SystemNodeList{grid, middle, pcc},
      SystemComponentList{source, inductance, resistance, inverter});
  Simulation simulation(name, Logger::Level::warn);
  simulation.setSystem(system);
  simulation.setDomain(Domain::SP);
  simulation.setSolverType(Solver::Type::NRP);
  simulation.setSolverAndComponentBehaviour(Solver::Behaviour::Initialization);
  simulation.setTimeStep(1.0);
  simulation.setFinalTime(2.0);
  simulation.doInitFromNodesAndTerminals(false);
  simulation.run();
  return system;
}

EMTDPPh3GFLStateSpaceValidation::SystemHandles
EMTDPPh3GFLStateSpaceValidation::buildEmtSystem(
    ModelKind model, const SystemTopology &powerFlow, Real gainScale,
    const std::shared_ptr<DataLogger> &logger) const {
  auto grid = SimNode<Real>::make("nGrid", PhaseType::ABC);
  auto middle = SimNode<Real>::make("nMiddle", PhaseType::ABC);
  auto pcc = SimNode<Real>::make("nPcc", PhaseType::ABC);
  auto source = EMT::Ph3::NetworkInjection::make("GridSource");
  auto inductance = EMT::Ph3::Inductor::make("GridInductor");
  inductance->setParameters(Math::singlePhaseParameterToThreePhase(
      mParameters.gridInductance));
  auto resistance = EMT::Ph3::Resistor::make("GridResistor");
  resistance->setParameters(Math::singlePhaseParameterToThreePhase(
      mParameters.gridResistance));

  Attribute<Real>::Ptr activePower;
  Attribute<Real>::Ptr reactivePower;
  SystemComponentList components{source, inductance, resistance};

  const Real kpPower = gainScale * mParameters.kpPowerBase;
  const Real kiPower = gainScale * mParameters.kiPowerBase;
  const Real kpCurrent = gainScale * mParameters.kpCurrentBase;
  const Real kiCurrent = gainScale * mParameters.kiCurrentBase;
  if (model == ModelKind::EmtVariable) {
    auto inverter = EMT::Ph3::SSN_GFL::make("Inverter");
    inverter->setParameters(
        mParameters.filterInductance, mParameters.filterCapacitance,
        mParameters.filterResistance, mParameters.couplingResistance,
        mParameters.omega, mParameters.kpPll, mParameters.kiPll,
        mParameters.powerCutoff, mParameters.activePowerReference,
        mParameters.reactivePowerReference, kpPower, kiPower, kpCurrent,
        kiCurrent);
    inverter->connect({EMT::SimNode::GND, pcc});
    activePower = inverter->attributeTyped<Real>("p_inst");
    reactivePower = inverter->attributeTyped<Real>("q_inst");
    components.push_back(inverter);
  } else {
    auto inverter = EMT::Ph3::SSN_GFL_Split::make("Inverter");
    inverter->setParameters(
        mParameters.filterInductance, mParameters.filterCapacitance,
        mParameters.filterResistance, mParameters.couplingResistance,
        mParameters.omega, mParameters.kpPll, mParameters.kiPll,
        mParameters.powerCutoff, mParameters.activePowerReference,
        mParameters.reactivePowerReference, kpPower, kiPower, kpCurrent,
        kiCurrent);
    inverter->connect({EMT::SimNode::GND, pcc});
    activePower = inverter->attributeTyped<Real>("p_inst");
    reactivePower = inverter->attributeTyped<Real>("q_inst");
    components.push_back(inverter);
  }

  source->connect({grid});
  inductance->connect({grid, middle});
  resistance->connect({middle, pcc});
  SystemTopology system(mParameters.frequency,
                        SystemNodeList{grid, middle, pcc}, components);
  system.initWithPowerflow(powerFlow, Domain::EMT);
  if (logger) {
    logger->logAttribute("p_inst", activePower);
    logger->logAttribute("q_inst", reactivePower);
  }
  return {system, source->attributeTyped<MatrixComp>("V_ref"), activePower,
          reactivePower};
}

EMTDPPh3GFLStateSpaceValidation::SystemHandles
EMTDPPh3GFLStateSpaceValidation::buildDpSystem(
    ModelKind model, const OperatingPoint &op, Real gainScale,
    const std::shared_ptr<DataLogger> &logger) const {
  auto pcc = SimNode<Complex>::make("nPcc", PhaseType::ABC);
  auto middle = SimNode<Complex>::make("nMiddle", PhaseType::ABC);
  auto grid = SimNode<Complex>::make("nGrid", PhaseType::ABC);
  // op is expressed in the power-invariant global dq0 convention, whose
  // balanced voltage magnitude equals the line-to-line RMS node convention.
  pcc->setInitialVoltage(op.terminal);
  middle->setInitialVoltage(op.intermediate);
  grid->setInitialVoltage(op.source);

  auto source = DP::Ph3::VoltageSource::make("GridSource");
  source->setParameters(balancedEnvelope(op.source), 0.0);
  auto inductance = DP::Ph3::Inductor::make("GridInductor");
  inductance->setParameters(Math::singlePhaseParameterToThreePhase(
      mParameters.gridInductance));
  auto resistance = DP::Ph3::Resistor::make("GridResistor");
  resistance->setParameters(Math::singlePhaseParameterToThreePhase(
      mParameters.gridResistance));

  Attribute<Real>::Ptr activePower;
  Attribute<Real>::Ptr reactivePower;
  SystemComponentList components;
  const Real kpPower = gainScale * mParameters.kpPowerBase;
  const Real kiPower = gainScale * mParameters.kiPowerBase;
  const Real kpCurrent = gainScale * mParameters.kpCurrentBase;
  const Real kiCurrent = gainScale * mParameters.kiCurrentBase;
  if (model == ModelKind::DpVariable) {
    auto inverter = DP::Ph3::AvVoltSourceInverterStateSpace::make("Inverter");
    inverter->setParameters(
        mParameters.filterInductance, mParameters.filterCapacitance,
        mParameters.filterResistance, mParameters.couplingResistance,
        mParameters.omega, mParameters.kpPll, mParameters.kiPll,
        mParameters.powerCutoff, mParameters.activePowerReference,
        mParameters.reactivePowerReference, kpPower, kiPower, kpCurrent,
        kiCurrent);
    inverter->connect({DP::SimNode::GND, pcc});
    activePower = inverter->attributeTyped<Real>("p_inst");
    reactivePower = inverter->attributeTyped<Real>("q_inst");
    components.push_back(inverter);
  } else {
    auto inverter = DP::Ph3::SSN_GFL_Split::make("Inverter");
    inverter->setParameters(
        mParameters.filterInductance, mParameters.filterCapacitance,
        mParameters.filterResistance, mParameters.couplingResistance,
        mParameters.omega, mParameters.kpPll, mParameters.kiPll,
        mParameters.powerCutoff, mParameters.activePowerReference,
        mParameters.reactivePowerReference, kpPower, kiPower, kpCurrent,
        kiCurrent);
    inverter->connect({DP::SimNode::GND, pcc});
    activePower = inverter->attributeTyped<Real>("p_inst");
    reactivePower = inverter->attributeTyped<Real>("q_inst");
    components.push_back(inverter);
  }

  resistance->connect({pcc, middle});
  inductance->connect({middle, grid});
  source->connect({DP::SimNode::GND, grid});
  components.push_back(resistance);
  components.push_back(inductance);
  components.push_back(source);
  if (logger) {
    logger->logAttribute("p_inst", activePower);
    logger->logAttribute("q_inst", reactivePower);
  }
  return {SystemTopology(mParameters.frequency,
                         SystemNodeList{pcc, middle, grid}, components),
          source->attributeTyped<MatrixComp>("V_ref"), activePower,
          reactivePower};
}

EMTDPPh3GFLStateSpaceValidation::SystemHandles
EMTDPPh3GFLStateSpaceValidation::buildSystem(
    ModelKind model, const SystemTopology &powerFlow, const OperatingPoint &op,
    Real gainScale, const std::shared_ptr<DataLogger> &logger) const {
  if (isEmt(model))
    return buildEmtSystem(model, powerFlow, gainScale, logger);
  return buildDpSystem(model, op, gainScale, logger);
}

ModalResult EMTDPPh3GFLStateSpaceValidation::extractOneStep(
    ModelKind model, const SystemTopology &powerFlow, const OperatingPoint &op,
    Real gainScale, Real timeStep, const String &study) const {
  const String name = "GFLValidation_" + study + "_" + fileToken(model) +
                      "_g" + gainToken(gainScale) + "_dt_" +
                      timeStepToken(timeStep);
  SystemHandles handles =
      buildSystem(model, powerFlow, op, gainScale, nullptr);
  Simulation simulation(name, Logger::Level::warn);
  simulation.setSystem(handles.system);
  simulation.setDomain(isEmt(model) ? Domain::EMT : Domain::DP);
  simulation.setSolverType(Solver::Type::MNA);
  simulation.doSystemMatrixRecomputation(true);
  simulation.doInitFromNodesAndTerminals(true);
  simulation.doStateSpaceExtraction(true);
  simulation.setTimeStep(timeStep);
  simulation.setFinalTime(timeStep);
  simulation.initialize();
  simulation.start();
  simulation.step();

  const auto &extractor = simulation.getStateSpaceExtractor();
  StateSpaceModalAnalysis modal(extractor);
  if (isEmt(model)) {
    modal.setAnalysisFrame(StateSpaceAnalysisFrame::GlobalDQ0);
    modal.setGlobalDq0Frame(mParameters.omega);
  }
  if (isSplit(model))
    modal.setPoleMapping(StateSpacePoleMapping::Logarithmic);
  modal.update();

  ModalResult result{modal.getDiscreteEigenvalues(),
                     modal.getContinuousEigenvalues(),
                     modal.getParticipationFactors(),
                     modal.getStateNames(),
                     extractor.getLastExtractionTime(),
                     extractor.getStateCount()};
  const UInt expectedStateCount =
      model == ModelKind::EmtVariable
          ? 17
          : model == ModelKind::EmtSplit
                ? 20
                : model == ModelKind::DpVariable ? 26 : 32;
  if (result.stateCount != expectedStateCount)
    throw std::runtime_error("Unexpected state count for " + modelName(model) +
                             ".");
  simulation.stop();
  return result;
}

VectorComp EMTDPPh3GFLStateSpaceValidation::calculateMonodromy(
    ModelKind model, const SystemTopology &powerFlow, const OperatingPoint &op,
    Real gainScale, Real timeStep) const {
  if (!isEmt(model))
    throw std::invalid_argument("Monodromy is only required for EMT models.");
  const UInt steps =
      static_cast<UInt>(std::llround(1.0 / (mParameters.frequency * timeStep)));
  const String name = "GFLValidation_monodromy_" + fileToken(model) + "_dt_" +
                      timeStepToken(timeStep);
  SystemHandles handles =
      buildSystem(model, powerFlow, op, gainScale, nullptr);
  Simulation simulation(name, Logger::Level::warn);
  simulation.setSystem(handles.system);
  simulation.setDomain(Domain::EMT);
  simulation.setSolverType(Solver::Type::MNA);
  simulation.doSystemMatrixRecomputation(true);
  simulation.doInitFromNodesAndTerminals(true);
  simulation.doStateSpaceExtraction(true);
  simulation.setTimeStep(timeStep);
  simulation.setFinalTime(steps * timeStep);
  simulation.initialize();
  simulation.start();

  Matrix transition;
  for (UInt step = 0; step < steps; ++step) {
    simulation.step();
    const Matrix &ad =
        simulation.getStateSpaceExtractor().getDiscreteStateMatrix();
    if (step == 0)
      transition = Matrix::Identity(ad.rows(), ad.cols());
    transition = ad * transition;
  }
  simulation.stop();
  return eigenvalues(transition);
}

TimeDomainResult EMTDPPh3GFLStateSpaceValidation::runTimeDomainCase(
    ModelKind model, const SystemTopology &powerFlow, const OperatingPoint &op,
    Real gainScale, Real timeStep, Real finalTime,
    const String &stabilityCase) const {
  const String name = "EMT_DP_Ph3_GFL_StateSpaceValidation_" +
                      fileToken(model) + "_" + stabilityCase + "_dt_" +
                      timeStepToken(timeStep);
  Logger::setLogDir((mOutputDirectory / "time_domain" / name).string());
  auto logger = DataLogger::make(name);
  SystemHandles handles = buildSystem(model, powerFlow, op, gainScale, logger);
  Simulation simulation(name, Logger::Level::warn);
  simulation.setSystem(handles.system);
  simulation.addLogger(logger);
  simulation.setDomain(isEmt(model) ? Domain::EMT : Domain::DP);
  simulation.setSolverType(Solver::Type::MNA);
  simulation.doSystemMatrixRecomputation(true);
  simulation.doInitFromNodesAndTerminals(true);
  simulation.setTimeStep(timeStep);
  simulation.setFinalTime(finalTime);
  simulation.initialize();

  const MatrixComp nominalVoltage = handles.sourceVoltageReference->get();
  const MatrixComp perturbedVoltage =
      (1.0 + mParameters.gridVoltagePulseRelative) * nominalVoltage;
  Bool perturbationApplied = false;
  Bool perturbationCleared = false;
  Real maximumActivePowerError = 0.0;
  Real maximumReactivePowerError = 0.0;
  simulation.start();
  while (simulation.time() < simulation.finalTime() - DOUBLE_EPSILON) {
    if (!perturbationApplied &&
        simulation.time() >= mParameters.perturbationTime - DOUBLE_EPSILON) {
      **handles.sourceVoltageReference = perturbedVoltage;
      perturbationApplied = true;
    }
    if (perturbationApplied && !perturbationCleared &&
        simulation.time() >= mParameters.perturbationTime +
                                 mParameters.perturbationDuration -
                                 DOUBLE_EPSILON) {
      **handles.sourceVoltageReference = nominalVoltage;
      perturbationCleared = true;
    }
    simulation.step();
    if (simulation.time() <=
        mParameters.perturbationTime + DOUBLE_EPSILON) {
      maximumActivePowerError =
          std::max(maximumActivePowerError,
                   std::abs(**handles.activePower -
                            mParameters.activePowerReference) /
                       std::max(1.0,
                                std::abs(mParameters.activePowerReference)));
      maximumReactivePowerError =
          std::max(maximumReactivePowerError,
                   std::abs(**handles.reactivePower -
                            mParameters.reactivePowerReference) /
                       std::max(1.0,
                                std::abs(mParameters.reactivePowerReference)));
    }
  }
  simulation.stop();
  return {(mOutputDirectory / "time_domain" / name / (name + ".csv"))
              .string(),
          maximumActivePowerError, maximumReactivePowerError};
}

Matrix EMTDPPh3GFLStateSpaceValidation::controllerDerivative(
    const Matrix &x, const Matrix &measurement, Real gainScale) const {
  const Complex voltage(measurement(0, 0), measurement(1, 0));
  const Complex current(measurement(3, 0), measurement(4, 0));
  const Complex rotation = std::exp(Complex(0.0, -x(0, 0)));
  const Complex voltageLocal = voltage * rotation;
  const Complex currentLocal = current * rotation;
  const Complex power = voltageLocal * std::conj(currentLocal);
  const Real kpPower = gainScale * mParameters.kpPowerBase;
  const Real kiPower = gainScale * mParameters.kiPowerBase;
  const Complex currentReference(
      kpPower * (mParameters.activePowerReference - x(2, 0)) +
          kiPower * x(4, 0),
      kpPower * (x(3, 0) - mParameters.reactivePowerReference) +
          kiPower * x(5, 0));

  Matrix derivative = Matrix::Zero(ControllerStateCount, 1);
  derivative(0, 0) = mParameters.kpPll * voltageLocal.imag() +
                     mParameters.kiPll * x(1, 0);
  derivative(1, 0) = voltageLocal.imag();
  derivative(2, 0) =
      mParameters.powerCutoff * (power.real() - x(2, 0));
  derivative(3, 0) =
      mParameters.powerCutoff * (power.imag() - x(3, 0));
  derivative(4, 0) = mParameters.activePowerReference - x(2, 0);
  derivative(5, 0) = x(3, 0) - mParameters.reactivePowerReference;
  derivative(6, 0) = currentReference.real() - currentLocal.real();
  derivative(7, 0) = currentReference.imag() - currentLocal.imag();
  return derivative;
}

Matrix EMTDPPh3GFLStateSpaceValidation::controllerOutput(
    const Matrix &x, const Matrix &measurement, Real gainScale) const {
  const Complex current(measurement(3, 0), measurement(4, 0));
  const Complex rotation = std::exp(Complex(0.0, -x(0, 0)));
  const Complex currentLocal = current * rotation;
  const Real kpPower = gainScale * mParameters.kpPowerBase;
  const Real kiPower = gainScale * mParameters.kiPowerBase;
  const Real kpCurrent = gainScale * mParameters.kpCurrentBase;
  const Real kiCurrent = gainScale * mParameters.kiCurrentBase;
  const Complex currentReference(
      kpPower * (mParameters.activePowerReference - x(2, 0)) +
          kiPower * x(4, 0),
      kpPower * (x(3, 0) - mParameters.reactivePowerReference) +
          kiPower * x(5, 0));
  const Complex voltageLocal =
      kpCurrent * (currentReference - currentLocal) +
      kiCurrent * Complex(x(6, 0), x(7, 0));
  const Complex voltageGlobal =
      voltageLocal * std::exp(Complex(0.0, x(0, 0)));
  Matrix output = Matrix::Zero(3, 1);
  output(0, 0) = voltageGlobal.real();
  output(1, 0) = voltageGlobal.imag();
  return output;
}

Matrix EMTDPPh3GFLStateSpaceValidation::physicalDerivative(
    const Matrix &x, const Matrix &bridgeVoltage) const {
  Matrix derivative = Matrix::Zero(PhysicalStateCount, 1);
  const auto addCrossCoupling = [this, &x](Matrix &value, UInt offset) {
    value(offset, 0) += mParameters.omega * x(offset + 1, 0);
    value(offset + 1, 0) -= mParameters.omega * x(offset, 0);
  };

  derivative.block(0, 0, 3, 1) =
      (x.block(3, 0, 3, 1) - x.block(6, 0, 3, 1)) /
      mParameters.filterCapacitance;
  derivative.block(3, 0, 3, 1) =
      (bridgeVoltage - x.block(0, 0, 3, 1) -
       mParameters.filterResistance * x.block(3, 0, 3, 1)) /
      mParameters.filterInductance;
  Matrix source = Matrix::Zero(3, 1);
  source(0, 0) = mParameters.voltageRmsLineToLine;
  derivative.block(6, 0, 3, 1) =
      (x.block(0, 0, 3, 1) - source -
       (mParameters.couplingResistance + mParameters.gridResistance) *
           x.block(6, 0, 3, 1)) /
      mParameters.gridInductance;
  addCrossCoupling(derivative, 0);
  addCrossCoupling(derivative, 3);
  addCrossCoupling(derivative, 6);
  return derivative;
}

ReferenceResult EMTDPPh3GFLStateSpaceValidation::buildReferences(
    const OperatingPoint &op, Real gainScale, Real timeStep) const {
  const Real angle = std::arg(op.capacitorVoltage);
  const Complex rotation = std::exp(Complex(0.0, -angle));
  const Complex current = op.gridCurrent * rotation;
  const Complex bridge = op.bridgeVoltage * rotation;
  const Real kiPower = gainScale * mParameters.kiPowerBase;
  const Real kiCurrent = gainScale * mParameters.kiCurrentBase;

  Matrix controller = Matrix::Zero(ControllerStateCount, 1);
  controller(0, 0) = angle;
  controller(2, 0) = mParameters.activePowerReference;
  controller(3, 0) = mParameters.reactivePowerReference;
  controller(4, 0) = current.real() / kiPower;
  controller(5, 0) = current.imag() / kiPower;
  controller(6, 0) = bridge.real() / kiCurrent;
  controller(7, 0) = bridge.imag() / kiCurrent;

  Matrix measurement = Matrix::Zero(6, 1);
  measurement(0, 0) = op.capacitorVoltage.real();
  measurement(1, 0) = op.capacitorVoltage.imag();
  measurement(3, 0) = op.gridCurrent.real();
  measurement(4, 0) = op.gridCurrent.imag();
  Matrix physical = Matrix::Zero(PhysicalStateCount, 1);
  physical(0, 0) = op.capacitorVoltage.real();
  physical(1, 0) = op.capacitorVoltage.imag();
  physical(3, 0) = op.filterCurrent.real();
  physical(4, 0) = op.filterCurrent.imag();
  physical(6, 0) = op.gridCurrent.real();
  physical(7, 0) = op.gridCurrent.imag();
  Matrix bridgeInput = Matrix::Zero(3, 1);
  bridgeInput(0, 0) = op.bridgeVoltage.real();
  bridgeInput(1, 0) = op.bridgeVoltage.imag();

  const auto controllerFunction = [this, gainScale](const Matrix &x,
                                                     const Matrix &u) {
    return controllerDerivative(x, u, gainScale);
  };
  const auto outputFunction = [this, gainScale](const Matrix &x,
                                                 const Matrix &u) {
    return controllerOutput(x, u, gainScale);
  };
  const auto physicalFunction = [this](const Matrix &x, const Matrix &u) {
    return physicalDerivative(x, u);
  };
  const Real equilibriumResidual =
      std::max(controllerFunction(controller, measurement).norm(),
               physicalFunction(physical, bridgeInput).norm());
  if (equilibriumResidual > 1e-4)
    throw std::runtime_error(
        "Analytical reference is not at the intended equilibrium; residual "
        "norm = " +
        std::to_string(equilibriumResidual) + ".");
  const Matrix ac = numericalJacobian(
      [&controllerFunction, &measurement](const Matrix &x) {
        return controllerFunction(x, measurement);
      },
      controller);
  const Matrix bc =
      numericalInputJacobian(controllerFunction, controller, measurement);
  const Matrix cc = numericalJacobian(
      [&outputFunction, &measurement](const Matrix &x) {
        return outputFunction(x, measurement);
      },
      controller);
  const Matrix dc =
      numericalInputJacobian(outputFunction, controller, measurement);
  const Matrix ap = numericalJacobian(
      [&physicalFunction, &bridgeInput](const Matrix &x) {
        return physicalFunction(x, bridgeInput);
      },
      physical);
  const Matrix bp =
      numericalInputJacobian(physicalFunction, physical, bridgeInput);

  Matrix measurementFromPhysical = Matrix::Zero(6, PhysicalStateCount);
  measurementFromPhysical.block(0, 0, 3, 3).setIdentity();
  measurementFromPhysical.block(3, 6, 3, 3).setIdentity();

  Matrix noDelay = Matrix::Zero(NoDelayStateCount, NoDelayStateCount);
  noDelay.block(0, 0, ControllerStateCount, ControllerStateCount) = ac;
  noDelay.block(0, ControllerStateCount, ControllerStateCount,
                PhysicalStateCount) = bc * measurementFromPhysical;
  noDelay.block(ControllerStateCount, 0, PhysicalStateCount,
                ControllerStateCount) = bp * cc;
  noDelay.block(ControllerStateCount, ControllerStateCount,
                PhysicalStateCount, PhysicalStateCount) =
      ap + bp * dc * measurementFromPhysical;

  Matrix controllerAd, controllerBd;
  Matrix physicalAd, physicalBd;
  trapezoidalMatrices(ac, bc, timeStep, controllerAd, controllerBd);
  trapezoidalMatrices(ap, bp, timeStep, physicalAd, physicalBd);
  Matrix split = Matrix::Zero(SplitReferenceStateCount,
                              SplitReferenceStateCount);
  const UInt physicalOffset = ControllerStateCount;
  const UInt delayOffset = NoDelayStateCount;
  split.block(physicalOffset, physicalOffset, PhysicalStateCount,
              PhysicalStateCount) = physicalAd;
  split.block(physicalOffset, delayOffset, PhysicalStateCount,
              DelayStateCount) = 2.0 * physicalBd;
  split.block(0, 0, ControllerStateCount, ControllerStateCount) = controllerAd;
  split.block(0, physicalOffset, ControllerStateCount, PhysicalStateCount) =
      controllerBd * measurementFromPhysical *
      (Matrix::Identity(PhysicalStateCount, PhysicalStateCount) + physicalAd);
  split.block(0, delayOffset, ControllerStateCount, DelayStateCount) =
      2.0 * controllerBd * measurementFromPhysical * physicalBd;
  split.block(delayOffset, 0, DelayStateCount, ControllerStateCount) =
      cc * controllerAd;
  split.block(delayOffset, physicalOffset, DelayStateCount,
              PhysicalStateCount) =
      cc * controllerBd * measurementFromPhysical *
          (Matrix::Identity(PhysicalStateCount, PhysicalStateCount) +
           physicalAd) +
      dc * measurementFromPhysical * physicalAd;
  split.block(delayOffset, delayOffset, DelayStateCount, DelayStateCount) =
      2.0 * (cc * controllerBd * measurementFromPhysical +
             dc * measurementFromPhysical) *
      physicalBd;

  return {noDelay, trapezoidalStateMatrix(noDelay, timeStep), split};
}

void EMTDPPh3GFLStateSpaceValidation::writeEigenvalues(
    const std::vector<EigenvalueRecord> &records) const {
  std::ofstream stream(mOutputDirectory / "eigenvalues.csv");
  stream << std::setprecision(std::numeric_limits<Real>::max_digits10)
         << "study,model,method,gain_scale,time_step_s,time_step_us,index,"
            "z_real,z_imag,lambda_real,lambda_imag\n";
  for (const auto &record : records) {
    for (Eigen::Index idx = 0; idx < record.continuous.rows(); ++idx) {
      const Complex z = idx < record.discrete.rows()
                            ? record.discrete(idx)
                            : Complex(std::numeric_limits<Real>::quiet_NaN(),
                                      std::numeric_limits<Real>::quiet_NaN());
      stream << record.study << ',' << record.model << ',' << record.method
             << ',' << record.gainScale << ',' << record.timeStep << ','
             << 1e6 * record.timeStep << ',' << idx << ',' << z.real() << ','
             << z.imag() << ',' << record.continuous(idx).real() << ','
             << record.continuous(idx).imag() << '\n';
    }
  }
}

void EMTDPPh3GFLStateSpaceValidation::writeSummary(
    const std::vector<SummaryRecord> &records) const {
  std::ofstream stream(mOutputDirectory / "summary.csv");
  stream << std::setprecision(std::numeric_limits<Real>::max_digits10)
         << "study,model,gain_scale,time_step_s,time_step_us,state_count,"
            "extraction_sample_time_s,max_real_lambda,spectral_radius,"
            "reference_error,"
            "split_to_no_delay_deviation,reference_error_z\n";
  for (const auto &record : records)
    stream << record.study << ',' << record.model << ',' << record.gainScale
           << ',' << record.timeStep << ',' << 1e6 * record.timeStep << ','
           << record.stateCount << ',' << record.extractionTime << ','
           << record.maxRealPart << ',' << record.spectralRadius << ','
           << record.referenceError << ',' << record.noDelayDeviation << ','
           << record.discreteReferenceError << '\n';
}

void EMTDPPh3GFLStateSpaceValidation::writeParticipation(
    const std::vector<ParticipationRecord> &records) const {
  std::ofstream stream(mOutputDirectory / "participation.csv");
  stream << std::setprecision(std::numeric_limits<Real>::max_digits10)
         << "model,mode_index,state_index,state_name,p_real,p_imag,p_abs\n";
  for (const auto &record : records)
    stream << record.model << ',' << record.mode << ',' << record.state << ','
           << record.stateName << ',' << record.value.real() << ','
           << record.value.imag() << ',' << std::abs(record.value) << '\n';
}

void EMTDPPh3GFLStateSpaceValidation::writeTimeDomainManifest(
    const std::vector<TimeDomainRecord> &records) const {
  std::ofstream stream(mOutputDirectory / "time_domain_manifest.csv");
  stream << std::setprecision(std::numeric_limits<Real>::max_digits10)
         << "model,stability_case,gain_scale,time_step_s,time_step_us,"
            "final_time_s,perturbation_time_s,perturbation_duration_s,"
            "perturbation_relative,pre_perturbation_p_error,"
            "pre_perturbation_q_error,log_path\n";
  for (const auto &record : records)
    stream << record.model << ',' << record.stabilityCase << ','
           << record.gainScale << ',' << record.timeStep << ','
           << 1e6 * record.timeStep << ',' << record.finalTime << ','
           << record.perturbationTime << ',' << record.perturbationDuration
           << ',' << record.perturbationRelative << ','
           << record.prePerturbationActivePowerError << ','
           << record.prePerturbationReactivePowerError << ','
           << record.logPath << '\n';
}

void EMTDPPh3GFLStateSpaceValidation::run() {
  std::filesystem::create_directories(mOutputDirectory);
  writeParameters(mOutputDirectory / "parameters.csv", mParameters);
  const OperatingPoint op = operatingPoint();
  const SystemTopology powerFlow = runPowerFlow("GFLValidation_PowerFlow");
  const std::array<ModelKind, 4> models = {
      ModelKind::EmtVariable, ModelKind::EmtSplit, ModelKind::DpVariable,
      ModelKind::DpSplit};
  const std::array<Real, 8> timeSteps = {1e-6,   5e-6,   10e-6,  50e-6,
                                         100e-6, 250e-6, 500e-6, 1e-3};
  std::vector<EigenvalueRecord> eigenvalueRecords;
  std::vector<SummaryRecord> summaryRecords;
  std::vector<ParticipationRecord> participationRecords;

  std::cout << "\n============================================================\n"
            << "EMT/DP Ph3 GFL state-space validation\n"
            << "============================================================\n"
            << "Topology: ideal source -> R/L grid -> averaged GFL inverter\n"
            << "Models: EMT/DP Ph3 variable and split SSN\n"
            << "Stable gain scale: " << mParameters.stableGainScale << "\n"
            << "Time-domain logging: "
            << (mRunTimeDomain ? "enabled" : "disabled (use --time-domain)")
            << "\n";

  for (const Real timeStep : timeSteps) {
    const ReferenceResult reference =
        buildReferences(op, mParameters.stableGainScale, timeStep);
    const VectorComp noDelayZ = eigenvalues(reference.noDelayAd);
    const VectorComp noDelayLambda =
        bilinearContinuousEigenvalues(noDelayZ, timeStep);
    const VectorComp splitZ = eigenvalues(reference.splitAd);
    const VectorComp splitLambda =
        logarithmicContinuousEigenvalues(splitZ, timeStep);
    eigenvalueRecords.push_back({"time_step_sweep", "analytical no delay",
                                 "trapezoidal", mParameters.stableGainScale,
                                 timeStep, noDelayZ, noDelayLambda});
    eigenvalueRecords.push_back({"time_step_sweep", "analytical split delay",
                                 "partitioned_trapezoidal",
                                 mParameters.stableGainScale, timeStep, splitZ,
                                 splitLambda});
    std::cout << "\n  dt = " << timeStepToken(timeStep) << "\n";

    if (timeStep == timeSteps.front()) {
      const std::vector<String> noDelayNames = {
          "psi",       "phi_pll", "p_filtered", "q_filtered",
          "phi_d",     "phi_q",   "gamma_d",    "gamma_q",
          "vc_d",      "vc_q",    "vc_0",       "if_d",
          "if_q",      "if_0",    "i_line_d",   "i_line_q",
          "i_line_0"};
      std::vector<String> splitNames = noDelayNames;
      splitNames.push_back("v_inv_delay_d");
      splitNames.push_back("v_inv_delay_q");
      splitNames.push_back("v_inv_delay_0");
      appendParticipationFromMatrix("analytical no delay", reference.noDelayAd,
                                    noDelayNames, participationRecords);
      appendParticipationFromMatrix("analytical split delay", reference.splitAd,
                                    splitNames, participationRecords);
    }

    for (const ModelKind model : models) {
      const ModalResult result = extractOneStep(
          model, powerFlow, op, mParameters.stableGainScale, timeStep,
          "time_step_sweep");
      const VectorComp referenceLambda =
          isSplit(model) ? splitLambda : noDelayLambda;
      const VectorComp referenceZ = isSplit(model) ? splitZ : noDelayZ;
      const Real referenceErrorZ = directedEigenvalueDistance(
          referenceZ, result.discreteEigenvalues);
      const Real referenceError = directedFiniteEigenvalueDistance(
          referenceLambda, result.continuousEigenvalues);
      const Real noDelayDeviation =
          isSplit(model)
              ? directedFiniteEigenvalueDistance(
                    noDelayLambda, result.continuousEigenvalues)
              : referenceError;
      eigenvalueRecords.push_back(
          {"time_step_sweep", modelName(model), "one_step",
           mParameters.stableGainScale, timeStep, result.discreteEigenvalues,
           result.continuousEigenvalues});
      summaryRecords.push_back(
          {"time_step_sweep", modelName(model), mParameters.stableGainScale,
           timeStep, result.stateCount, result.extractionTime,
           maximumRealPart(result.continuousEigenvalues),
           spectralRadius(result.discreteEigenvalues), referenceErrorZ,
           referenceError, noDelayDeviation});
      std::cout << "    " << modelName(model)
                << ": reference error |Delta z| = " << referenceErrorZ
                << ", finite |Delta lambda| = " << referenceError
                << " 1/s, max Re(lambda) = "
                << maximumRealPart(result.continuousEigenvalues) << " 1/s\n";

      if (timeStep == timeSteps.front()) {
        for (UInt mode = 0; mode < result.stateCount; ++mode) {
          for (UInt state = 0; state < result.stateCount; ++state) {
            const String stateName =
                state < result.stateNames.size()
                    ? result.stateNames[state]
                    : "x" + std::to_string(state);
            participationRecords.push_back(
                {modelName(model), mode, state, stateName,
                 result.participationFactors(state, mode)});
          }
        }
      }

      if (isEmt(model)) {
        try {
          const VectorComp multipliers = calculateMonodromy(
              model, powerFlow, op, mParameters.stableGainScale, timeStep);
          VectorComp floquetLambda(multipliers.rows());
          const Real period = 1.0 / mParameters.frequency;
          for (Eigen::Index idx = 0; idx < multipliers.rows(); ++idx)
            floquetLambda(idx) = std::log(multipliers(idx)) / period;
          eigenvalueRecords.push_back(
              {"time_step_sweep", modelName(model), "monodromy",
               mParameters.stableGainScale, timeStep, multipliers,
               floquetLambda});
        } catch (const std::exception &error) {
          std::cerr << "    " << modelName(model)
                    << " monodromy skipped: " << error.what() << '\n';
        }
      }
    }

    // Preserve all completed cases even if a later, deliberately unstable
    // case cannot be evaluated numerically.
    writeEigenvalues(eigenvalueRecords);
    writeSummary(summaryRecords);
    writeParticipation(participationRecords);
  }

  // A compact gain sweep shows that the two time-domain gain cases lie on
  // opposite sides of a continuously moving filter/grid mode rather than
  // being isolated, hand-picked operating points.
  constexpr Real gainSweepTimeStep = 1e-6;
  const std::array<Real, 6> gainScales = {3.5, 3.6, 3.7, 3.75, 3.8, 3.9};
  for (const Real gainScale : gainScales) {
    const ReferenceResult reference =
        buildReferences(op, gainScale, gainSweepTimeStep);
    const VectorComp noDelayZ = eigenvalues(reference.noDelayAd);
    const VectorComp noDelayLambda =
        bilinearContinuousEigenvalues(noDelayZ, gainSweepTimeStep);
    const VectorComp splitZ = eigenvalues(reference.splitAd);
    const VectorComp splitLambda =
        logarithmicContinuousEigenvalues(splitZ, gainSweepTimeStep);
    eigenvalueRecords.push_back(
        {"gain_sweep", "analytical no delay", "trapezoidal", gainScale,
         gainSweepTimeStep, noDelayZ, noDelayLambda});
    eigenvalueRecords.push_back(
        {"gain_sweep", "analytical split delay", "partitioned_trapezoidal",
         gainScale, gainSweepTimeStep, splitZ, splitLambda});
    for (const ModelKind model : models) {
      const ModalResult result = extractOneStep(
          model, powerFlow, op, gainScale, gainSweepTimeStep, "gain_sweep");
      const VectorComp &referenceZ = isSplit(model) ? splitZ : noDelayZ;
      const VectorComp &referenceLambda =
          isSplit(model) ? splitLambda : noDelayLambda;
      eigenvalueRecords.push_back(
          {"gain_sweep", modelName(model), "one_step", gainScale,
           gainSweepTimeStep, result.discreteEigenvalues,
           result.continuousEigenvalues});
      summaryRecords.push_back(
          {"gain_sweep", modelName(model), gainScale, gainSweepTimeStep,
           result.stateCount, result.extractionTime,
           maximumRealPart(result.continuousEigenvalues),
           spectralRadius(result.discreteEigenvalues),
           directedEigenvalueDistance(referenceZ,
                                      result.discreteEigenvalues),
           directedFiniteEigenvalueDistance(referenceLambda,
                                            result.continuousEigenvalues),
           isSplit(model)
               ? directedFiniteEigenvalueDistance(
                     noDelayLambda, result.continuousEigenvalues)
               : 0.0});
    }
  }
  writeEigenvalues(eigenvalueRecords);
  writeSummary(summaryRecords);

  for (const Real gainScale : {mParameters.stableGainScale,
                               mParameters.unstableGainScale}) {
    const String stabilityCase = gainScale == mParameters.stableGainScale
                                     ? "stable"
                                     : "unstable";
    for (const Real timeStep : {1e-6, 1e-3}) {
      const ReferenceResult reference = buildReferences(op, gainScale, timeStep);
      const VectorComp noDelayLambda = bilinearContinuousEigenvalues(
          eigenvalues(reference.noDelayAd), timeStep);
      const VectorComp splitLambda = logarithmicContinuousEigenvalues(
          eigenvalues(reference.splitAd), timeStep);
      for (const ModelKind model : models) {
        const ModalResult result = extractOneStep(
            model, powerFlow, op, gainScale, timeStep,
            "stability_" + stabilityCase);
        const VectorComp referenceLambda =
            isSplit(model) ? splitLambda : noDelayLambda;
        const VectorComp referenceZ =
            isSplit(model) ? eigenvalues(reference.splitAd)
                           : eigenvalues(reference.noDelayAd);
        eigenvalueRecords.push_back(
            {"stability_" + stabilityCase, modelName(model), "one_step",
             gainScale, timeStep, result.discreteEigenvalues,
             result.continuousEigenvalues});
        summaryRecords.push_back(
            {"stability_" + stabilityCase, modelName(model), gainScale,
             timeStep, result.stateCount, result.extractionTime,
             maximumRealPart(result.continuousEigenvalues),
             spectralRadius(result.discreteEigenvalues),
             directedEigenvalueDistance(referenceZ,
                                        result.discreteEigenvalues),
             directedFiniteEigenvalueDistance(referenceLambda,
                                              result.continuousEigenvalues),
             isSplit(model)
                 ? directedFiniteEigenvalueDistance(
                       noDelayLambda, result.continuousEigenvalues)
                 : 0.0});
      }
      writeEigenvalues(eigenvalueRecords);
      writeSummary(summaryRecords);
    }
  }

  std::vector<TimeDomainRecord> timeDomainRecords;
  if (mRunTimeDomain) {
    struct TimeDomainCase {
      String name;
      Real gainScale;
      Real timeStep;
      Real finalTime;
    };
    // Every case first runs undisturbed for one fundamental period. A short,
    // balanced grid-voltage pulse then excites the free response without using
    // model-specific state coordinates. The 1 us cases validate the
    // gain-induced modal stability change. The 50 and 100 us cases bracket the
    // split-delay stability crossing while resolving its approximately 436 Hz
    // critical pair. A 1 ms waveform would undersample the 1.9 kHz
    // gain-sensitive mode and can overflow; it remains a modal diagnostic.
    const std::array<TimeDomainCase, 4> timeDomainCases = {
        TimeDomainCase{"stable", mParameters.stableGainScale, 1e-6,
                       mParameters.gainCaseFinalTime},
        TimeDomainCase{"unstable", mParameters.unstableGainScale, 1e-6,
                       mParameters.gainCaseFinalTime},
        TimeDomainCase{"delay_stable", mParameters.stableGainScale, 50e-6,
                       mParameters.largeStepFinalTime},
        TimeDomainCase{"large_step", mParameters.stableGainScale, 100e-6,
                       mParameters.largeStepFinalTime}};
    for (const auto &timeDomainCase : timeDomainCases) {
      for (const ModelKind model : models) {
        const TimeDomainResult result = runTimeDomainCase(
            model, powerFlow, op, timeDomainCase.gainScale,
            timeDomainCase.timeStep, timeDomainCase.finalTime,
            timeDomainCase.name);
        timeDomainRecords.push_back(
            {modelName(model), timeDomainCase.name,
             timeDomainCase.gainScale, timeDomainCase.timeStep,
             timeDomainCase.finalTime, mParameters.perturbationTime,
             mParameters.perturbationDuration,
             mParameters.gridVoltagePulseRelative,
             result.prePerturbationActivePowerError,
             result.prePerturbationReactivePowerError, result.logPath});
        writeTimeDomainManifest(timeDomainRecords);
      }
    }
  }

  writeEigenvalues(eigenvalueRecords);
  writeSummary(summaryRecords);
  writeParticipation(participationRecords);
  writeTimeDomainManifest(timeDomainRecords);
  std::cout << "Results written to " << mOutputDirectory.string() << "\n";
  if (!mRunTimeDomain)
    std::cout
        << "Use --time-domain to additionally generate Study 3 traces.\n";
}

int main(int argc, char **argv) {
  Bool runTimeDomain = false;
  for (Int idx = 1; idx < argc; ++idx) {
    if (String(argv[idx]) == "--time-domain" || String(argv[idx]) == "--all")
      runTimeDomain = true;
  }
  EMTDPPh3GFLStateSpaceValidation example(runTimeDomain);
  example.run();
  return 0;
}
