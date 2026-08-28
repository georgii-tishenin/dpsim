// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include <DPsim.h>

#include <dpsim-models/DP/DP_Ph3_AvVoltSourceInverterStateSpace.h>
#include <dpsim-models/DP/DP_Ph3_SSN_GFL_Split.h>
#include <dpsim-models/EMT/EMT_Ph3_AvVoltSourceInverterStateSpace.h>
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
#include <utility>
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

enum class ModelKind {
  EmtVariable,
  EmtLegacyVariable,
  EmtSplit,
  DpVariable,
  DpSplit
};

struct Parameters {
  // Benchmark definition. All impedance-like physical parameters and all PI
  // gains are specified on this power-invariant dq per-unit base. The derived
  // SI quantities below are the only values passed to the component models.
  Real frequency = 50.0;
  Real voltageBaseRmsLineToLine = 3.3e3;
  Real apparentPowerBase = 12.0e6;
  Real activePowerReferencePu = 0.8;
  Real reactivePowerReferencePu = 0.1;

  // The external Thevenin impedance is referred to the converter terminals
  // and represents the aggregate transformer and upstream network impedance.
  Real shortCircuitRatio = 2.5;
  Real gridXOverR = 10.0;
  // Filter values retain the normalized physical design of the original
  // benchmark. This separates the change of rating/grid strength from a
  // controller or filter redesign.
  Real filterInductiveReactancePu =
      2.0 * PI * 50.0 * 2.0e-3 /
      (400.0 * 400.0 / std::hypot(10000.0, 5000.0));
  Real filterCapacitiveSusceptancePu =
      2.0 * PI * 50.0 * 10.0e-6 *
      (400.0 * 400.0 / std::hypot(10000.0, 5000.0));
  Real filterResistancePu =
      0.2 / (400.0 * 400.0 / std::hypot(10000.0, 5000.0));
  Real couplingResistancePu = filterResistancePu;

  // These are the original stable-case controller gains expressed on the
  // power-invariant dq base. No controller equations or tuning relationships
  // are changed when the converter voltage and power rating are rescaled.
  Real kpPllPuPerSecond = 0.25 * 400.0;
  Real kiPllPuPerSecondSquared = 0.2 * 400.0;
  Real powerMeasurementBandwidthHz = 50.0;
  Real kpPowerPu = 3.7 * 0.05 * 400.0;
  Real kiPowerPuPerSecond = 3.7 * 0.2 * 400.0;
  Real kpCurrentPu =
      3.7 * 0.25 * std::hypot(10000.0, 5000.0) / (400.0 * 400.0);
  Real kiCurrentPuPerSecond =
      3.7 * 1.0 * std::hypot(10000.0, 5000.0) / (400.0 * 400.0);

  // Derived bases and SI parameters.
  Real omega = 0.0;
  Real currentBase = 0.0;
  Real impedanceBase = 0.0;
  Real inductanceBase = 0.0;
  Real capacitanceBase = 0.0;
  Real voltageRmsLineToLine = 0.0;
  Real gridResistancePu = 0.0;
  Real gridReactancePu = 0.0;
  Real gridResistance = 0.0;
  Real gridInductance = 0.0;
  Real filterInductance = 0.0;
  Real filterCapacitance = 0.0;
  Real filterResistance = 0.0;
  Real couplingResistance = 0.0;
  Real kpPll = 0.0;
  Real kiPll = 0.0;
  Real powerCutoff = 0.0;
  Real activePowerReference = 0.0;
  Real reactivePowerReference = 0.0;
  Real kpPowerBase = 0.0;
  Real kiPowerBase = 0.0;
  Real kpCurrentBase = 0.0;
  Real kiCurrentBase = 0.0;

  // Study 3 varies only the outer active/reactive-power proportional gain.
  // All integral, current-controller, and PLL gains remain fixed.
  Real stableGainScale = 0.4;
  Real unstableGainScale = 0.45;
  Real perturbationTime = 0.0;
  Real perturbationDuration = 0.5e-3;
  Real gridVoltagePulseRelative = 1.0e-3;
  Real gainCaseFinalTime = perturbationTime + 0.2;
  Real largeStepFinalTime = perturbationTime + 0.15;

  Parameters() { updateDerived(); }

  void updateDerived() {
    omega = 2.0 * PI * frequency;
    voltageRmsLineToLine = voltageBaseRmsLineToLine;
    currentBase = apparentPowerBase / voltageBaseRmsLineToLine;
    impedanceBase = voltageBaseRmsLineToLine * voltageBaseRmsLineToLine /
                    apparentPowerBase;
    inductanceBase = impedanceBase / omega;
    capacitanceBase = 1.0 / (omega * impedanceBase);

    const Real gridImpedancePu = 1.0 / shortCircuitRatio;
    gridResistancePu =
        gridImpedancePu / std::sqrt(1.0 + gridXOverR * gridXOverR);
    gridReactancePu = gridXOverR * gridResistancePu;
    gridResistance = gridResistancePu * impedanceBase;
    gridInductance = gridReactancePu * inductanceBase;
    filterInductance = filterInductiveReactancePu * inductanceBase;
    filterCapacitance =
        filterCapacitiveSusceptancePu * capacitanceBase;
    filterResistance = filterResistancePu * impedanceBase;
    couplingResistance = couplingResistancePu * impedanceBase;

    activePowerReference = activePowerReferencePu * apparentPowerBase;
    reactivePowerReference = reactivePowerReferencePu * apparentPowerBase;
    powerCutoff = 2.0 * PI * powerMeasurementBandwidthHz;

    kpPll = kpPllPuPerSecond / voltageBaseRmsLineToLine;
    kiPll = kiPllPuPerSecondSquared / voltageBaseRmsLineToLine;

    kpPowerBase = kpPowerPu / voltageBaseRmsLineToLine;
    kiPowerBase = kiPowerPuPerSecond / voltageBaseRmsLineToLine;
    kpCurrentBase = kpCurrentPu * impedanceBase;
    kiCurrentBase = kiCurrentPuPerSecond * impedanceBase;
  }
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
  String study;
  String model;
  Real timeStep;
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
  case ModelKind::EmtLegacyVariable:
    return "EMT legacy variable";
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
  case ModelKind::EmtLegacyVariable:
    return "emt_legacy_variable";
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
  return model == ModelKind::EmtVariable ||
         model == ModelKind::EmtLegacyVariable ||
         model == ModelKind::EmtSplit;
}

Bool isSplit(ModelKind model) {
  return model == ModelKind::EmtSplit || model == ModelKind::DpSplit;
}

String timeStepToken(Real timeStep) {
  const Real microseconds = 1e6 * timeStep;
  if (microseconds < 1.0)
    return std::to_string(static_cast<Int>(std::llround(1e9 * timeStep))) +
           "ns";
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

Matrix parkTransformDQ0(Real theta) {
  Matrix transform(3, 3);
  const Real k = std::sqrt(2.0 / 3.0);
  const Real k0 = 1.0 / std::sqrt(3.0);
  transform.row(0) << k * std::cos(theta),
      k * std::cos(theta - 2.0 * PI / 3.0),
      k * std::cos(theta + 2.0 * PI / 3.0);
  transform.row(1) << -k * std::sin(theta),
      -k * std::sin(theta - 2.0 * PI / 3.0),
      -k * std::sin(theta + 2.0 * PI / 3.0);
  transform.row(2) << k0, k0, k0;
  return transform;
}

Matrix parkTransformDerivativeDQ0(Real theta, Real omega) {
  Matrix derivative = Matrix::Zero(3, 3);
  const Real k = std::sqrt(2.0 / 3.0);
  derivative.row(0) << -omega * k * std::sin(theta),
      -omega * k * std::sin(theta - 2.0 * PI / 3.0),
      -omega * k * std::sin(theta + 2.0 * PI / 3.0);
  derivative.row(1) << -omega * k * std::cos(theta),
      -omega * k * std::cos(theta - 2.0 * PI / 3.0),
      -omega * k * std::cos(theta + 2.0 * PI / 3.0);
  return derivative;
}

Matrix globalDq0Transform(Real theta) {
  Matrix transform = Matrix::Identity(NoDelayStateCount, NoDelayStateCount);
  const Matrix park = parkTransformDQ0(theta);
  for (const UInt offset : {8u, 11u, 14u})
    transform.block(offset, offset, 3, 3) = park;
  return transform;
}

Matrix globalDq0TransformDerivative(Real theta, Real omega) {
  Matrix derivative = Matrix::Zero(NoDelayStateCount, NoDelayStateCount);
  const Matrix parkDerivative = parkTransformDerivativeDQ0(theta, omega);
  for (const UInt offset : {8u, 11u, 14u})
    derivative.block(offset, offset, 3, 3) = parkDerivative;
  return derivative;
}

Matrix frozenAbcTrapezoidalStep(const Matrix &globalDqA, Real omega,
                                Real time, Real timeStep) {
  const Matrix transformNow = globalDq0Transform(omega * time);
  const Matrix transformNext = globalDq0Transform(omega * (time + timeStep));
  const Matrix transformDerivative =
      globalDq0TransformDerivative(omega * time, omega);
  const Matrix nativeA = transformNow.transpose() * globalDqA * transformNow -
                         transformNow.transpose() * transformDerivative;
  const Matrix identity =
      Matrix::Identity(NoDelayStateCount, NoDelayStateCount);
  const Matrix nativeAd =
      (identity - 0.5 * timeStep * nativeA)
          .partialPivLu()
          .solve(identity + 0.5 * timeStep * nativeA);
  return transformNext * nativeAd * transformNow.transpose();
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
    const String &study, const String &model, Real timeStep,
    const Matrix &matrix,
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
      records.push_back({study, model, timeStep, mode, state, stateNames[state],
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
  const Real phaseCurrentRmsBase =
      p.apparentPowerBase /
      (std::sqrt(3.0) * p.voltageRmsLineToLine);
  std::ofstream stream(path);
  stream << std::setprecision(std::numeric_limits<Real>::max_digits10)
         << "parameter,value,unit\n"
         << "frequency," << p.frequency << ",Hz\n"
         << "apparent_power_base," << p.apparentPowerBase << ",VA\n"
         << "voltage_rms_ll," << p.voltageRmsLineToLine << ",V\n"
         << "current_base_dq," << p.currentBase << ",A\n"
         << "impedance_base," << p.impedanceBase << ",ohm\n"
         << "inductance_base," << p.inductanceBase << ",H\n"
         << "capacitance_base," << p.capacitanceBase << ",F\n"
         << "active_power_reference_pu," << p.activePowerReferencePu
         << ",pu\n"
         << "reactive_power_reference_pu," << p.reactivePowerReferencePu
         << ",pu\n"
         << "short_circuit_ratio," << p.shortCircuitRatio << ",-\n"
         << "grid_x_over_r," << p.gridXOverR << ",-\n"
         << "grid_resistance_pu," << p.gridResistancePu << ",pu\n"
         << "grid_reactance_pu," << p.gridReactancePu << ",pu\n"
         << "filter_inductive_reactance_pu,"
         << p.filterInductiveReactancePu << ",pu\n"
         << "filter_capacitive_susceptance_pu,"
         << p.filterCapacitiveSusceptancePu << ",pu\n"
         << "filter_resistance_pu," << p.filterResistancePu << ",pu\n"
         << "coupling_resistance_pu," << p.couplingResistancePu
         << ",pu\n"
         << "grid_resistance," << p.gridResistance << ",ohm\n"
         << "grid_inductance," << p.gridInductance << ",H\n"
         << "filter_inductance," << p.filterInductance << ",H\n"
         << "filter_capacitance," << p.filterCapacitance << ",F\n"
         << "filter_resistance," << p.filterResistance << ",ohm\n"
         << "coupling_resistance," << p.couplingResistance << ",ohm\n"
         << "current_cross_coupling_coefficient,"
         << p.omega * p.filterInductance << ",ohm\n"
         << "current_cross_coupling_default_enabled,0,-\n"
         << "kp_pll," << p.kpPll << ",rad_per_V_s\n"
         << "ki_pll," << p.kiPll << ",rad_per_V_s2\n"
         << "kp_pll_pu," << p.kpPllPuPerSecond << ",pu_per_s\n"
         << "ki_pll_pu," << p.kiPllPuPerSecondSquared
         << ",pu_per_s2\n"
         << "power_measurement_bandwidth,"
         << p.powerMeasurementBandwidthHz << ",Hz\n"
         << "power_cutoff," << p.powerCutoff << ",rad_per_s\n"
         << "active_power_reference," << p.activePowerReference << ",W\n"
         << "reactive_power_reference," << p.reactivePowerReference
         << ",var\n"
         << "kp_power_base," << p.kpPowerBase << ",A_per_W\n"
         << "ki_power_base," << p.kiPowerBase << ",A_per_W_s\n"
         << "kp_current_base," << p.kpCurrentBase << ",V_per_A\n"
         << "ki_current_base," << p.kiCurrentBase << ",V_per_A_s\n"
         << "kp_power_pu," << p.kpPowerPu << ",pu_A_per_pu_P\n"
         << "ki_power_pu," << p.kiPowerPuPerSecond
         << ",pu_A_per_pu_P_s\n"
         << "kp_current_pu," << p.kpCurrentPu << ",pu_V_per_pu_A\n"
         << "ki_current_pu," << p.kiCurrentPuPerSecond
         << ",pu_V_per_pu_A_s\n"
         << "stable_kp_power_multiplier," << p.stableGainScale << ",-\n"
         << "unstable_kp_power_multiplier," << p.unstableGainScale
         << ",-\n";
  const Real resonanceFrequency =
      std::sqrt((p.filterInductance + p.gridInductance) /
                (p.filterInductance * p.gridInductance *
                 p.filterCapacitance)) /
      (2.0 * PI);
  stream << "phase_current_rms_base," << phaseCurrentRmsBase << ",A\n"
         << "dq_current_base," << p.currentBase << ",A\n"
         << "kp_power_base_normalized,"
         << p.kpPowerBase * p.apparentPowerBase / p.currentBase
         << ",pu_A_per_pu_P\n"
         << "ki_power_base_normalized,"
         << p.kiPowerBase * p.apparentPowerBase / p.currentBase
         << ",pu_A_per_pu_P_s\n"
         << "kp_current_base_normalized,"
         << p.kpCurrentBase * p.currentBase / p.voltageRmsLineToLine
         << ",pu_V_per_pu_A\n"
         << "ki_current_base_normalized,"
         << p.kiCurrentBase * p.currentBase / p.voltageRmsLineToLine
         << ",pu_V_per_pu_A_s\n"
         << "operating_apparent_power,"
         << std::hypot(p.activePowerReference, p.reactivePowerReference)
         << ",VA\n"
         << "operating_point_scr,"
         << p.voltageRmsLineToLine * p.voltageRmsLineToLine /
                (std::hypot(p.gridResistance,
                            p.omega * p.gridInductance) *
                 p.apparentPowerBase)
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

namespace {

class EMTDPPh3GFLStateSpaceValidation {
public:
  explicit EMTDPPh3GFLStateSpaceValidation(
      Bool runTimeDomain, Bool enableCurrentCrossCoupling = false,
      const String &resultDirectory =
          "EMT_DP_Ph3_GFL_StateSpaceValidation")
      : mRunTimeDomain(runTimeDomain),
        mEnableCurrentCrossCoupling(enableCurrentCrossCoupling),
        mOutputDirectory(std::filesystem::path("logs") / resultDirectory) {
    if (mEnableCurrentCrossCoupling) {
      // The original 0.4/0.45 pair is already unstable with decoupling because
      // a separate electrical mode crosses first.  Use a lower pair that
      // brackets that cross-coupled stability boundary without retuning any
      // other controller coefficient.
      mParameters.shortCircuitRatio = 1.4;
      mParameters.stableGainScale = 0.2;
      mParameters.unstableGainScale = 0.26;
      mParameters.updateDerived();
    }
  }

  void run();
  void runBenchmarkSelectionScan();
  void runWeakGridStudy();
  void runCandidateCheck();
  void runEmtVariableDiagnostics();

private:
  struct SystemHandles {
    SystemTopology system;
    Attribute<MatrixComp>::Ptr sourceVoltageReference;
    Attribute<Real>::Ptr activePower;
    Attribute<Real>::Ptr reactivePower;
    std::shared_ptr<EMT::Ph3::SSN_GFL> emtVariableInverter;
    std::shared_ptr<EMT::VTypeVariableSSNComp> emtVariableBase;
    std::shared_ptr<EMT::Ph3::Inductor> emtGridInductor;
  };

  OperatingPoint operatingPoint() const;
  SystemTopology runPowerFlow(const String &name) const;
  SystemHandles buildEmtSystem(ModelKind model, const SystemTopology &powerFlow,
                               Real gainScale,
                               const std::shared_ptr<DataLogger> &logger,
                               Bool enableCurrentCrossCoupling = false) const;
  SystemHandles buildDpSystem(ModelKind model,
                              const SystemTopology &powerFlow, Real gainScale,
                              const std::shared_ptr<DataLogger> &logger,
                              Bool enableCurrentCrossCoupling = false) const;
  SystemHandles buildSystem(ModelKind model, const SystemTopology &powerFlow,
                            const OperatingPoint &op, Real gainScale,
                            const std::shared_ptr<DataLogger> &logger,
                            Bool enableCurrentCrossCoupling = false) const;

  ModalResult extractOneStep(ModelKind model, const SystemTopology &powerFlow,
                             const OperatingPoint &op, Real gainScale,
                             Real timeStep, const String &study,
                             Bool useCompactHistoryState = false,
                             Bool enableCurrentCrossCoupling = false) const;
  VectorComp calculateMonodromy(ModelKind model,
                                const SystemTopology &powerFlow,
                                const OperatingPoint &op, Real gainScale,
                                Real timeStep,
                                UInt warmupPeriods = 0,
                                Bool useCompactHistoryState = false,
                                Bool enableCurrentCrossCoupling = false) const;
  TimeDomainResult
  runTimeDomainCase(ModelKind model, const SystemTopology &powerFlow,
                    const OperatingPoint &op, Real gainScale, Real timeStep,
                    Real finalTime, const String &stabilityCase,
                    Real perturbationRelative =
                        std::numeric_limits<Real>::quiet_NaN()) const;

  Matrix controllerDerivative(const Matrix &controllerState,
                              const Matrix &measurement,
                              Real gainScale) const;
  Matrix controllerOutput(const Matrix &controllerState,
                          const Matrix &measurement, Real gainScale,
                          Bool enableCurrentCrossCoupling = false) const;
  Matrix physicalDerivative(const Matrix &physicalState,
                            const Matrix &bridgeVoltage) const;
  ReferenceResult buildReferences(const OperatingPoint &op, Real gainScale,
                                  Real timeStep,
                                  Bool enableCurrentCrossCoupling = false) const;

  void writeEigenvalues(const std::vector<EigenvalueRecord> &records) const;
  void writeSummary(const std::vector<SummaryRecord> &records) const;
  void writeParticipation(
      const std::vector<ParticipationRecord> &records) const;
  void writeTimeDomainManifest(
      const std::vector<TimeDomainRecord> &records) const;

  Parameters mParameters;
  Bool mRunTimeDomain;
  Bool mEnableCurrentCrossCoupling;
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
    const std::shared_ptr<DataLogger> &logger,
    Bool enableCurrentCrossCoupling) const {
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
  std::shared_ptr<EMT::Ph3::SSN_GFL> emtVariableInverter;
  std::shared_ptr<EMT::VTypeVariableSSNComp> emtVariableBase;
  SystemComponentList components{source, inductance, resistance};

  const Real kpPower = gainScale * mParameters.kpPowerBase;
  const Real kiPower = mParameters.kiPowerBase;
  const Real kpCurrent = mParameters.kpCurrentBase;
  const Real kiCurrent = mParameters.kiCurrentBase;
  if (model == ModelKind::EmtVariable) {
    auto inverter = EMT::Ph3::SSN_GFL::make("Inverter");
    emtVariableInverter = inverter;
    emtVariableBase = inverter;
    inverter->setParameters(
        mParameters.filterInductance, mParameters.filterCapacitance,
        mParameters.filterResistance, mParameters.couplingResistance,
        mParameters.omega, mParameters.kpPll, mParameters.kiPll,
        mParameters.powerCutoff, mParameters.activePowerReference,
        mParameters.reactivePowerReference, kpPower, kiPower, kpCurrent,
        kiCurrent);
    inverter->setEnableCurrentCrossCoupling(enableCurrentCrossCoupling);
    inverter->connect({EMT::SimNode::GND, pcc});
    activePower = inverter->attributeTyped<Real>("p_inst");
    reactivePower = inverter->attributeTyped<Real>("q_inst");
    components.push_back(inverter);
  } else if (model == ModelKind::EmtLegacyVariable) {
    auto inverter =
        EMT::Ph3::AvVoltSourceInverterStateSpace::make("Inverter");
    emtVariableBase = inverter;
    inverter->setParameters(
        mParameters.filterInductance, mParameters.filterCapacitance,
        mParameters.filterResistance, mParameters.couplingResistance,
        mParameters.omega, mParameters.kpPll, mParameters.kiPll,
        mParameters.powerCutoff, mParameters.activePowerReference,
        mParameters.reactivePowerReference, kpPower, kiPower, kpCurrent,
        kiCurrent);
    if (enableCurrentCrossCoupling)
      throw std::invalid_argument(
          "Cross-coupling comparison is not implemented for the legacy EMT "
          "diagnostic model.");
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
    inverter->setEnableCurrentCrossCoupling(enableCurrentCrossCoupling);
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
          reactivePower, emtVariableInverter, emtVariableBase, inductance};
}

EMTDPPh3GFLStateSpaceValidation::SystemHandles
EMTDPPh3GFLStateSpaceValidation::buildDpSystem(
    ModelKind model, const SystemTopology &powerFlow, Real gainScale,
    const std::shared_ptr<DataLogger> &logger,
    Bool enableCurrentCrossCoupling) const {
  auto pcc = SimNode<Complex>::make("nPcc", PhaseType::ABC);
  auto middle = SimNode<Complex>::make("nMiddle", PhaseType::ABC);
  auto grid = SimNode<Complex>::make("nGrid", PhaseType::ABC);

  auto source = DP::Ph3::VoltageSource::make("GridSource");
  source->setParameters(
      balancedEnvelope(Complex(mParameters.voltageRmsLineToLine, 0.0)), 0.0);
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
  const Real kiPower = mParameters.kiPowerBase;
  const Real kpCurrent = mParameters.kpCurrentBase;
  const Real kiCurrent = mParameters.kiCurrentBase;
  if (model == ModelKind::DpVariable) {
    auto inverter = DP::Ph3::AvVoltSourceInverterStateSpace::make("Inverter");
    inverter->setParameters(
        mParameters.filterInductance, mParameters.filterCapacitance,
        mParameters.filterResistance, mParameters.couplingResistance,
        mParameters.omega, mParameters.kpPll, mParameters.kiPll,
        mParameters.powerCutoff, mParameters.activePowerReference,
        mParameters.reactivePowerReference, kpPower, kiPower, kpCurrent,
        kiCurrent);
    inverter->setEnableCurrentCrossCoupling(enableCurrentCrossCoupling);
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
    inverter->setEnableCurrentCrossCoupling(enableCurrentCrossCoupling);
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
  SystemTopology system(mParameters.frequency,
                        SystemNodeList{pcc, middle, grid}, components);
  // Use the solved SP node voltages for DP initialization, just as the EMT
  // branch does. The analytical operating point is reserved for references.
  system.initWithPowerflow(powerFlow, Domain::DP);
  return {system, source->attributeTyped<MatrixComp>("V_ref"), activePower,
          reactivePower, nullptr, nullptr, nullptr};
}

EMTDPPh3GFLStateSpaceValidation::SystemHandles
EMTDPPh3GFLStateSpaceValidation::buildSystem(
    ModelKind model, const SystemTopology &powerFlow, const OperatingPoint &op,
    Real gainScale, const std::shared_ptr<DataLogger> &logger,
    Bool enableCurrentCrossCoupling) const {
  if (isEmt(model))
    return buildEmtSystem(model, powerFlow, gainScale, logger,
                          enableCurrentCrossCoupling);
  return buildDpSystem(model, powerFlow, gainScale, logger,
                       enableCurrentCrossCoupling);
}

ModalResult EMTDPPh3GFLStateSpaceValidation::extractOneStep(
    ModelKind model, const SystemTopology &powerFlow, const OperatingPoint &op,
    Real gainScale, Real timeStep, const String &study,
    Bool useCompactHistoryState, Bool enableCurrentCrossCoupling) const {
  const String name = "GFLValidation_" + study + "_" + fileToken(model) +
                      "_g" + gainToken(gainScale) + "_dt_" +
                      timeStepToken(timeStep);
  SystemHandles handles =
      buildSystem(model, powerFlow, op, gainScale, nullptr,
                  enableCurrentCrossCoupling);
  if (handles.emtVariableBase)
    handles.emtVariableBase->useAugmentedPhysicalStateExtraction(
        !useCompactHistoryState);
  else if (useCompactHistoryState)
    throw std::invalid_argument(
        "Compact history-state extraction is only available for an EMT "
        "variable SSN component.");
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
  if ((model == ModelKind::EmtVariable ||
       model == ModelKind::EmtLegacyVariable) &&
      !useCompactHistoryState)
    modal.setReduceAuxiliaryStates(true);
  modal.update();

  ModalResult result{modal.getDiscreteEigenvalues(),
                     modal.getContinuousEigenvalues(),
                     modal.getParticipationFactors(),
                     modal.getStateNames(),
                     extractor.getLastExtractionTime(),
                     static_cast<UInt>(modal.getStateNames().size())};
  UInt expectedStateCount = 32;
  if (model == ModelKind::EmtVariable ||
      model == ModelKind::EmtLegacyVariable)
    expectedStateCount = 17;
  else if (model == ModelKind::EmtSplit)
    expectedStateCount = 20;
  else if (model == ModelKind::DpVariable)
    expectedStateCount = 26;
  if (result.stateCount != expectedStateCount)
    throw std::runtime_error("Unexpected state count for " + modelName(model) +
                             ".");
  simulation.stop();
  return result;
}

VectorComp EMTDPPh3GFLStateSpaceValidation::calculateMonodromy(
    ModelKind model, const SystemTopology &powerFlow, const OperatingPoint &op,
    Real gainScale, Real timeStep, UInt warmupPeriods,
    Bool useCompactHistoryState, Bool enableCurrentCrossCoupling) const {
  if (!isEmt(model))
    throw std::invalid_argument("Monodromy is only required for EMT models.");
  const UInt steps =
      static_cast<UInt>(std::llround(1.0 / (mParameters.frequency * timeStep)));
  const String name = "GFLValidation_monodromy_" + fileToken(model) + "_dt_" +
                      timeStepToken(timeStep);
  SystemHandles handles =
      buildSystem(model, powerFlow, op, gainScale, nullptr,
                  enableCurrentCrossCoupling);
  if (handles.emtVariableBase)
    handles.emtVariableBase->useAugmentedPhysicalStateExtraction(
        !useCompactHistoryState);
  else if (useCompactHistoryState)
    throw std::invalid_argument(
        "Compact history-state monodromy is only available for an EMT "
        "variable SSN component.");
  Simulation simulation(name, Logger::Level::warn);
  simulation.setSystem(handles.system);
  simulation.setDomain(Domain::EMT);
  simulation.setSolverType(Solver::Type::MNA);
  simulation.doSystemMatrixRecomputation(true);
  simulation.doInitFromNodesAndTerminals(true);
  simulation.doStateSpaceExtraction(true);
  simulation.setTimeStep(timeStep);
  simulation.setFinalTime((warmupPeriods + 1) * steps * timeStep);
  simulation.initialize();
  simulation.start();

  for (UInt step = 0; step < warmupPeriods * steps; ++step)
    simulation.step();

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
    const String &stabilityCase, Real perturbationRelative) const {
  const String name = "EMT_DP_Ph3_GFL_StateSpaceValidation_" +
                      fileToken(model) + "_" + stabilityCase + "_dt_" +
                      timeStepToken(timeStep);
  Logger::setLogDir((mOutputDirectory / "time_domain" / name).string());
  auto logger = DataLogger::make(name);
  SystemHandles handles = buildSystem(model, powerFlow, op, gainScale, logger,
                                      mEnableCurrentCrossCoupling);
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
  const Real appliedPerturbation =
      std::isfinite(perturbationRelative)
          ? perturbationRelative
          : mParameters.gridVoltagePulseRelative;
  const MatrixComp perturbedVoltage =
      (1.0 + appliedPerturbation) * nominalVoltage;
  // Study 3 starts from the common power-flow initialization and applies the
  // scheduled pulse before the first MNA step. This avoids mixing an arbitrary
  // pre-disturbance transient with the free-response plotting interval.
  **handles.sourceVoltageReference = perturbedVoltage;
  Bool perturbationApplied = true;
  Bool perturbationCleared = false;
  const Bool hasPrePerturbationInterval = mParameters.perturbationTime > 0.0;
  Real maximumActivePowerError =
      hasPrePerturbationInterval
          ? 0.0
          : std::numeric_limits<Real>::quiet_NaN();
  Real maximumReactivePowerError =
      hasPrePerturbationInterval
          ? 0.0
          : std::numeric_limits<Real>::quiet_NaN();
  simulation.start();
  while (simulation.time() < simulation.finalTime() - DOUBLE_EPSILON) {
    if (perturbationApplied && !perturbationCleared &&
        simulation.time() >= mParameters.perturbationTime +
                                 mParameters.perturbationDuration -
                                 DOUBLE_EPSILON) {
      **handles.sourceVoltageReference = nominalVoltage;
      perturbationCleared = true;
    }
    simulation.step();
    if (hasPrePerturbationInterval &&
        simulation.time() <= mParameters.perturbationTime + DOUBLE_EPSILON) {
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
  const Real kiPower = mParameters.kiPowerBase;
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
    const Matrix &x, const Matrix &measurement, Real gainScale,
    Bool enableCurrentCrossCoupling) const {
  const Complex current(measurement(3, 0), measurement(4, 0));
  const Complex filterCurrent(measurement(6, 0), measurement(7, 0));
  const Complex rotation = std::exp(Complex(0.0, -x(0, 0)));
  const Complex currentLocal = current * rotation;
  const Real kpPower = gainScale * mParameters.kpPowerBase;
  const Real kiPower = mParameters.kiPowerBase;
  const Real kpCurrent = mParameters.kpCurrentBase;
  const Real kiCurrent = mParameters.kiCurrentBase;
  const Complex currentReference(
      kpPower * (mParameters.activePowerReference - x(2, 0)) +
          kiPower * x(4, 0),
      kpPower * (x(3, 0) - mParameters.reactivePowerReference) +
          kiPower * x(5, 0));
  const Complex voltageLocal =
      kpCurrent * (currentReference - currentLocal) +
      kiCurrent * Complex(x(6, 0), x(7, 0)) +
      (enableCurrentCrossCoupling
           ? Complex(0.0, mParameters.omega *
                              mParameters.filterInductance) *
                 filterCurrent * rotation
           : Complex(0.0, 0.0));
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
    const OperatingPoint &op, Real gainScale, Real timeStep,
    Bool enableCurrentCrossCoupling) const {
  const Real angle = std::arg(op.capacitorVoltage);
  const Complex rotation = std::exp(Complex(0.0, -angle));
  const Complex current = op.gridCurrent * rotation;
  const Complex bridge = op.bridgeVoltage * rotation;
  const Complex filterCurrent = op.filterCurrent * rotation;
  const Real kiPower = mParameters.kiPowerBase;
  const Real kiCurrent = mParameters.kiCurrentBase;

  Matrix controller = Matrix::Zero(ControllerStateCount, 1);
  controller(0, 0) = angle;
  controller(2, 0) = mParameters.activePowerReference;
  controller(3, 0) = mParameters.reactivePowerReference;
  controller(4, 0) = current.real() / kiPower;
  controller(5, 0) = current.imag() / kiPower;
  const Complex crossVoltage =
      enableCurrentCrossCoupling
          ? Complex(0.0, mParameters.omega *
                             mParameters.filterInductance) *
                filterCurrent
          : Complex(0.0, 0.0);
  controller(6, 0) = (bridge.real() - crossVoltage.real()) / kiCurrent;
  controller(7, 0) = (bridge.imag() - crossVoltage.imag()) / kiCurrent;

  Matrix measurement = Matrix::Zero(9, 1);
  measurement(0, 0) = op.capacitorVoltage.real();
  measurement(1, 0) = op.capacitorVoltage.imag();
  measurement(3, 0) = op.gridCurrent.real();
  measurement(4, 0) = op.gridCurrent.imag();
  measurement(6, 0) = op.filterCurrent.real();
  measurement(7, 0) = op.filterCurrent.imag();
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
  const auto outputFunction = [this, gainScale,
                               enableCurrentCrossCoupling](const Matrix &x,
                                                           const Matrix &u) {
    return controllerOutput(x, u, gainScale, enableCurrentCrossCoupling);
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

  Matrix measurementFromPhysical = Matrix::Zero(9, PhysicalStateCount);
  measurementFromPhysical.block(0, 0, 3, 3).setIdentity();
  measurementFromPhysical.block(3, 6, 3, 3).setIdentity();
  measurementFromPhysical.block(6, 3, 3, 3).setIdentity();

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
  trapezoidalMatrices(ac, bc, timeStep, controllerAd, controllerBd);

  // Assemble the split reference component by component. At finite time steps,
  // discretizing the already closed nine-state physical subsystem is not
  // equivalent to the simulator's separate inverter/grid trapezoidal maps and
  // algebraic Rc closure.
  constexpr UInt inverterStateCount = 6;
  constexpr UInt gridStateCount = 3;
  Matrix inverterA = Matrix::Zero(inverterStateCount, inverterStateCount);
  Matrix inverterTerminalB = Matrix::Zero(inverterStateCount, 3);
  Matrix inverterBridgeB = Matrix::Zero(inverterStateCount, 3);
  inverterA.block(0, 0, 3, 3) =
      -1.0 / (mParameters.filterCapacitance *
              mParameters.couplingResistance) *
      Matrix::Identity(3, 3);
  inverterA.block(0, 3, 3, 3) =
      1.0 / mParameters.filterCapacitance * Matrix::Identity(3, 3);
  inverterA.block(3, 0, 3, 3) =
      -1.0 / mParameters.filterInductance * Matrix::Identity(3, 3);
  inverterA.block(3, 3, 3, 3) =
      -mParameters.filterResistance / mParameters.filterInductance *
      Matrix::Identity(3, 3);
  inverterTerminalB.block(0, 0, 3, 3) =
      1.0 / (mParameters.filterCapacitance *
             mParameters.couplingResistance) *
      Matrix::Identity(3, 3);
  inverterBridgeB.block(3, 0, 3, 3) =
      1.0 / mParameters.filterInductance * Matrix::Identity(3, 3);
  for (const UInt offset : {0U, 3U}) {
    inverterA(offset, offset + 1) += mParameters.omega;
    inverterA(offset + 1, offset) -= mParameters.omega;
  }

  Matrix gridA = -mParameters.gridResistance / mParameters.gridInductance *
                 Matrix::Identity(gridStateCount, gridStateCount);
  gridA(0, 1) += mParameters.omega;
  gridA(1, 0) -= mParameters.omega;
  const Matrix gridTerminalB =
      1.0 / mParameters.gridInductance * Matrix::Identity(3, 3);

  Matrix inverterAd, inverterBd;
  Matrix gridAd, gridBd;
  trapezoidalMatrices(inverterA, inverterTerminalB, timeStep, inverterAd,
                      inverterBd);
  trapezoidalMatrices(gridA, gridTerminalB, timeStep, gridAd, gridBd);
  const Matrix inverterLhs =
      Matrix::Identity(inverterStateCount, inverterStateCount) -
      0.5 * timeStep * inverterA;
  const Matrix inverterBridgeBd = inverterLhs.fullPivLu().solve(
      0.5 * timeStep * inverterBridgeB);

  Matrix terminalFromPhysical = Matrix::Zero(3, PhysicalStateCount);
  terminalFromPhysical.block(0, 0, 3, 3).setIdentity();
  terminalFromPhysical.block(0, inverterStateCount, 3, 3) =
      -mParameters.couplingResistance * Matrix::Identity(3, 3);
  Matrix terminalInput = Matrix::Zero(PhysicalStateCount, 3);
  terminalInput.block(0, 0, inverterStateCount, 3) = inverterBd;
  terminalInput.block(inverterStateCount, 0, gridStateCount, 3) = gridBd;
  Matrix delayInput = Matrix::Zero(PhysicalStateCount, DelayStateCount);
  delayInput.block(0, 0, inverterStateCount, DelayStateCount) =
      2.0 * inverterBridgeBd;

  const Matrix algebraicInverse =
      (Matrix::Identity(PhysicalStateCount, PhysicalStateCount) -
       terminalInput * terminalFromPhysical)
          .fullPivLu()
          .solve(Matrix::Identity(PhysicalStateCount, PhysicalStateCount));
  const Matrix physicalFromDelay = algebraicInverse * delayInput;

  Matrix historyAd = Matrix::Zero(PhysicalStateCount, PhysicalStateCount);
  historyAd.block(0, 0, inverterStateCount, inverterStateCount) = inverterAd;
  historyAd.block(inverterStateCount, inverterStateCount, gridStateCount,
                  gridStateCount) = gridAd;
  Matrix historyTerminal = Matrix::Zero(PhysicalStateCount, 3);
  historyTerminal.block(0, 0, inverterStateCount, 3) =
      (inverterAd + Matrix::Identity(inverterStateCount, inverterStateCount)) *
      inverterBd;
  historyTerminal.block(inverterStateCount, 0, gridStateCount, 3) =
      (gridAd + Matrix::Identity(gridStateCount, gridStateCount)) * gridBd;
  Matrix historyDelay =
      Matrix::Zero(PhysicalStateCount, DelayStateCount);
  historyDelay.block(0, 0, inverterStateCount, DelayStateCount) =
      2.0 * inverterAd * inverterBridgeBd;

  const Matrix terminalFromHistory = terminalFromPhysical * algebraicInverse;
  const Matrix terminalFromDelay =
      terminalFromPhysical * physicalFromDelay;
  const Matrix physicalAd =
      historyAd + historyTerminal * terminalFromHistory;
  const Matrix physicalDelay =
      historyDelay + historyTerminal * terminalFromDelay;
  const Matrix measurementFromHistory =
      measurementFromPhysical * algebraicInverse;
  const Matrix measurementFromDelay =
      measurementFromPhysical * physicalFromDelay;
  const Matrix controllerInputUpdate =
      controllerAd * controllerBd + controllerBd;
  const Matrix delayMeasurementGain = cc * controllerBd + dc;

  Matrix split = Matrix::Zero(SplitReferenceStateCount,
                              SplitReferenceStateCount);
  const UInt physicalOffset = ControllerStateCount;
  const UInt delayOffset = NoDelayStateCount;
  split.block(physicalOffset, physicalOffset, PhysicalStateCount,
              PhysicalStateCount) = physicalAd;
  split.block(physicalOffset, delayOffset, PhysicalStateCount,
              DelayStateCount) = physicalDelay;
  split.block(0, 0, ControllerStateCount, ControllerStateCount) = controllerAd;
  split.block(0, physicalOffset, ControllerStateCount, PhysicalStateCount) =
      controllerInputUpdate * measurementFromHistory;
  split.block(0, delayOffset, ControllerStateCount, DelayStateCount) =
      controllerInputUpdate * measurementFromDelay;
  split.block(delayOffset, 0, DelayStateCount, ControllerStateCount) =
      cc;
  split.block(delayOffset, physicalOffset, DelayStateCount,
              PhysicalStateCount) =
      delayMeasurementGain * measurementFromHistory;
  split.block(delayOffset, delayOffset, DelayStateCount, DelayStateCount) =
      delayMeasurementGain * measurementFromDelay;

  return {noDelay, trapezoidalStateMatrix(noDelay, timeStep), split};
}

void EMTDPPh3GFLStateSpaceValidation::writeEigenvalues(
    const std::vector<EigenvalueRecord> &records) const {
  std::ofstream stream(mOutputDirectory / "eigenvalues.csv");
  stream << std::setprecision(std::numeric_limits<Real>::max_digits10)
         << "study,model,method,kp_power_multiplier,time_step_s,time_step_us,index,"
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
         << "study,model,kp_power_multiplier,time_step_s,time_step_us,state_count,"
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
         << "study,model,time_step_s,time_step_us,mode_index,state_index,"
            "state_name,p_real,p_imag,p_abs\n";
  for (const auto &record : records)
    stream << record.study << ',' << record.model << ',' << record.timeStep
           << ',' << 1e6 * record.timeStep << ',' << record.mode << ','
           << record.state << ',' << record.stateName << ','
           << record.value.real() << ',' << record.value.imag() << ','
           << std::abs(record.value) << '\n';
}

void EMTDPPh3GFLStateSpaceValidation::writeTimeDomainManifest(
    const std::vector<TimeDomainRecord> &records) const {
  std::ofstream stream(mOutputDirectory / "time_domain_manifest.csv");
  stream << std::setprecision(std::numeric_limits<Real>::max_digits10)
         << "model,stability_case,kp_power_multiplier,time_step_s,time_step_us,"
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
  // Use round engineering time steps, with a 5 us refinement around the
  // observed split-model stability boundary.  Floquet analysis below is
  // additionally restricted to steps that close exactly over one 50 Hz
  // period; one-step extraction and analytical references have no such
  // restriction.
  const std::vector<Real> timeSteps =
      mEnableCurrentCrossCoupling
          ? std::vector<Real>{1e-6,   5e-6,   10e-6, 20e-6, 25e-6, 30e-6,
                              35e-6,  40e-6,  45e-6, 50e-6, 55e-6, 60e-6,
                              65e-6,  70e-6,  75e-6, 80e-6, 85e-6, 90e-6,
                              95e-6,  100e-6, 250e-6, 500e-6, 1e-3}
          : std::vector<Real>{1e-6,   5e-6,   10e-6, 50e-6, 55e-6, 60e-6,
                              65e-6,  70e-6,  75e-6, 80e-6, 85e-6, 90e-6,
                              95e-6,  100e-6, 250e-6, 500e-6, 1e-3};
  std::vector<EigenvalueRecord> eigenvalueRecords;
  std::vector<SummaryRecord> summaryRecords;
  std::vector<ParticipationRecord> participationRecords;

  std::cout << "\n============================================================\n"
            << "EMT/DP Ph3 GFL state-space validation\n"
            << "============================================================\n"
            << "Topology: ideal source -> R/L grid -> averaged GFL inverter\n"
            << "Models: EMT/DP Ph3 variable and split SSN\n"
            << "Stable Kp,P multiplier: " << mParameters.stableGainScale
            << "\nCurrent cross-coupling: "
            << (mEnableCurrentCrossCoupling ? "enabled" : "disabled")
            << "\n"
            << "Time-domain logging: "
            << (mRunTimeDomain ? "enabled" : "disabled (use --time-domain)")
            << "\n";

  for (const Real timeStep : timeSteps) {
    const ReferenceResult reference = buildReferences(
        op, mParameters.stableGainScale, timeStep,
        mEnableCurrentCrossCoupling);
    const VectorComp noDelayZ = eigenvalues(reference.noDelayAd);
    const VectorComp noDelayLambda =
        bilinearContinuousEigenvalues(noDelayZ, timeStep);
    const VectorComp splitZ = eigenvalues(reference.splitAd);
    const VectorComp splitLambda =
        logarithmicContinuousEigenvalues(splitZ, timeStep);
    eigenvalueRecords.push_back({"time_step_sweep", "analytical no delay",
                                 "trapezoidal", mParameters.stableGainScale,
                                 timeStep, noDelayZ, noDelayLambda});
    eigenvalueRecords.push_back({"time_step_sweep", "exact split companion",
                                 "partitioned_trapezoidal",
                                 mParameters.stableGainScale, timeStep, splitZ,
                                 splitLambda});
    std::cout << "\n  dt = " << timeStepToken(timeStep) << "\n";

    const Bool storeParticipation =
        timeStep == timeSteps.front() ||
        (mEnableCurrentCrossCoupling
             ? (std::abs(timeStep - 40e-6) < 1e-12 ||
                std::abs(timeStep - 45e-6) < 1e-12)
             : (std::abs(timeStep - 60e-6) < 1e-12 ||
                std::abs(timeStep - 80e-6) < 1e-12));
    if (storeParticipation) {
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
      appendParticipationFromMatrix("time_step_sweep", "analytical no delay",
                                    timeStep, reference.noDelayAd, noDelayNames,
                                    participationRecords);
      appendParticipationFromMatrix("time_step_sweep",
                                    "exact split companion", timeStep,
                                    reference.splitAd, splitNames,
                                    participationRecords);
    }

    for (const ModelKind model : models) {
      const ModalResult result = extractOneStep(
          model, powerFlow, op, mParameters.stableGainScale, timeStep,
          "time_step_sweep", false, mEnableCurrentCrossCoupling);
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
      if (storeParticipation) {
        for (UInt mode = 0; mode < result.stateCount; ++mode) {
          for (UInt state = 0; state < result.stateCount; ++state) {
            const String stateName =
                state < result.stateNames.size()
                    ? result.stateNames[state]
                    : "x" + std::to_string(state);
            participationRecords.push_back(
                {"time_step_sweep", modelName(model), timeStep, mode, state,
                 stateName,
                 result.participationFactors(state, mode)});
          }
        }
      }

      const Real period = 1.0 / mParameters.frequency;
      const Real stepsPerPeriod = period / timeStep;
      const Bool closesFundamentalPeriod =
          std::abs(stepsPerPeriod - std::round(stepsPerPeriod)) < 1e-9;
      if (isEmt(model) && closesFundamentalPeriod) {
        try {
          const VectorComp multipliers = calculateMonodromy(
              model, powerFlow, op, mParameters.stableGainScale, timeStep, 0,
              false, mEnableCurrentCrossCoupling);
          VectorComp floquetLambda(multipliers.rows());
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

  // Only the outer power-loop proportional gain is varied. The compact sweep
  // shows that the time-domain cases bracket a continuously moving
  // converter-grid interaction mode rather than isolated operating points.
  constexpr Real gainSweepTimeStep = 1e-6;
  const std::vector<Real> gainScales =
      mEnableCurrentCrossCoupling
          ? std::vector<Real>{0.1, 0.15, 0.19, 0.2, 0.225, 0.25, 0.255,
                              0.26, 0.275, 0.3, 0.35}
          : std::vector<Real>{0.35, 0.375, 0.4, 0.4125, 0.425, 0.4375, 0.45,
                              0.475};
  for (const Real gainScale : gainScales) {
    const ReferenceResult reference = buildReferences(
        op, gainScale, gainSweepTimeStep, mEnableCurrentCrossCoupling);
    const VectorComp noDelayZ = eigenvalues(reference.noDelayAd);
    const VectorComp noDelayLambda =
        bilinearContinuousEigenvalues(noDelayZ, gainSweepTimeStep);
    const VectorComp splitZ = eigenvalues(reference.splitAd);
    const VectorComp splitLambda =
        logarithmicContinuousEigenvalues(splitZ, gainSweepTimeStep);
    eigenvalueRecords.push_back(
        {"kp_power_sweep", "analytical no delay", "trapezoidal", gainScale,
         gainSweepTimeStep, noDelayZ, noDelayLambda});
    eigenvalueRecords.push_back(
        {"kp_power_sweep", "exact split companion", "partitioned_trapezoidal",
         gainScale, gainSweepTimeStep, splitZ, splitLambda});
    for (const ModelKind model : models) {
      const ModalResult result = extractOneStep(
          model, powerFlow, op, gainScale, gainSweepTimeStep,
          "kp_power_sweep", false, mEnableCurrentCrossCoupling);
      const VectorComp &referenceZ = isSplit(model) ? splitZ : noDelayZ;
      const VectorComp &referenceLambda =
          isSplit(model) ? splitLambda : noDelayLambda;
      eigenvalueRecords.push_back(
          {"kp_power_sweep", modelName(model), "one_step", gainScale,
           gainSweepTimeStep, result.discreteEigenvalues,
           result.continuousEigenvalues});
      summaryRecords.push_back(
          {"kp_power_sweep", modelName(model), gainScale, gainSweepTimeStep,
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

  // Additional controller study: repeat the same one-parameter gain sweep
  // with nominal-frequency filter-inductor decoupling enabled.  Keeping this
  // in a separate study preserves the original controller as the baseline and
  // isolates the effect of +j*omega_N*L_f*i_f,dq from timestep effects.
  if (!mEnableCurrentCrossCoupling) {
  for (const Real gainScale : gainScales) {
    const ReferenceResult reference = buildReferences(
        op, gainScale, gainSweepTimeStep, true);
    const VectorComp noDelayZ = eigenvalues(reference.noDelayAd);
    const VectorComp noDelayLambda =
        bilinearContinuousEigenvalues(noDelayZ, gainSweepTimeStep);
    const VectorComp splitZ = eigenvalues(reference.splitAd);
    const VectorComp splitLambda =
        logarithmicContinuousEigenvalues(splitZ, gainSweepTimeStep);
    eigenvalueRecords.push_back(
        {"current_cross_coupling_kp_sweep", "analytical no delay",
         "trapezoidal", gainScale, gainSweepTimeStep, noDelayZ,
         noDelayLambda});
    eigenvalueRecords.push_back(
        {"current_cross_coupling_kp_sweep", "exact split companion",
         "partitioned_trapezoidal", gainScale, gainSweepTimeStep, splitZ,
         splitLambda});

    const Bool storeCrossCouplingParticipation =
        std::abs(gainScale - mParameters.stableGainScale) < 1e-12;
    if (storeCrossCouplingParticipation) {
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
      appendParticipationFromMatrix(
          "current_cross_coupling_kp_sweep", "analytical no delay",
          gainSweepTimeStep, reference.noDelayAd, noDelayNames,
          participationRecords);
      appendParticipationFromMatrix(
          "current_cross_coupling_kp_sweep", "exact split companion",
          gainSweepTimeStep, reference.splitAd, splitNames,
          participationRecords);
    }

    for (const ModelKind model : models) {
      const ModalResult result = extractOneStep(
          model, powerFlow, op, gainScale, gainSweepTimeStep,
          "current_cross_coupling_kp_sweep", false, true);
      const VectorComp &referenceZ = isSplit(model) ? splitZ : noDelayZ;
      const VectorComp &referenceLambda =
          isSplit(model) ? splitLambda : noDelayLambda;
      eigenvalueRecords.push_back(
          {"current_cross_coupling_kp_sweep", modelName(model), "one_step",
           gainScale, gainSweepTimeStep, result.discreteEigenvalues,
           result.continuousEigenvalues});
      summaryRecords.push_back(
          {"current_cross_coupling_kp_sweep", modelName(model), gainScale,
           gainSweepTimeStep, result.stateCount, result.extractionTime,
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
      if (storeCrossCouplingParticipation) {
        for (UInt mode = 0; mode < result.stateCount; ++mode) {
          for (UInt state = 0; state < result.stateCount; ++state) {
            const String stateName =
                state < result.stateNames.size()
                    ? result.stateNames[state]
                    : "x" + std::to_string(state);
            participationRecords.push_back(
                {"current_cross_coupling_kp_sweep", modelName(model),
                 gainSweepTimeStep, mode, state, stateName,
                 result.participationFactors(state, mode)});
          }
        }
      }
    }
  }
  }
  writeEigenvalues(eigenvalueRecords);
  writeSummary(summaryRecords);
  writeParticipation(participationRecords);

  for (const Real gainScale : {mParameters.stableGainScale,
                               mParameters.unstableGainScale}) {
    const String stabilityCase = gainScale == mParameters.stableGainScale
                                     ? "stable"
                                     : "unstable";
    for (const Real timeStep : {1e-6, 1e-3}) {
      const ReferenceResult reference = buildReferences(
          op, gainScale, timeStep, mEnableCurrentCrossCoupling);
      const VectorComp noDelayLambda = bilinearContinuousEigenvalues(
          eigenvalues(reference.noDelayAd), timeStep);
      const VectorComp splitLambda = logarithmicContinuousEigenvalues(
          eigenvalues(reference.splitAd), timeStep);
      for (const ModelKind model : models) {
        const ModalResult result = extractOneStep(
            model, powerFlow, op, gainScale, timeStep,
            "stability_" + stabilityCase, false,
            mEnableCurrentCrossCoupling);
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
    // A short balanced grid-voltage pulse is present from t = 0 and excites
    // the response without using model-specific state coordinates. The 1 us
    // cases validate the Kp,P-induced modal stability change. The additional
    // cross-coupled 35 and 45 us cases straddle its refined split-delay
    // boundary; the non-cross-coupled benchmark keeps its 60 and 80 us pair.
    // A 1 ms waveform would undersample the approximately 1.1 kHz branch and
    // can overflow; it remains a modal diagnostic only.
    const std::array<TimeDomainCase, 6> timeDomainCases = {
        TimeDomainCase{"stable",
                       mEnableCurrentCrossCoupling
                           ? 0.19
                           : mParameters.stableGainScale,
                       1e-6,
                       mEnableCurrentCrossCoupling
                           ? 1.5
                           : mParameters.gainCaseFinalTime},
        TimeDomainCase{"unstable",
                       mEnableCurrentCrossCoupling
                           ? 0.255
                           : mParameters.unstableGainScale,
                       1e-6,
                       mEnableCurrentCrossCoupling
                           ? 1.5
                           : mParameters.gainCaseFinalTime},
        TimeDomainCase{"delay_stable", mParameters.stableGainScale,
                       mEnableCurrentCrossCoupling ? 30e-6 : 50e-6,
                       mParameters.largeStepFinalTime},
        TimeDomainCase{"near_transition_stable",
                       mParameters.stableGainScale,
                       mEnableCurrentCrossCoupling ? 35e-6 : 60e-6,
                       mEnableCurrentCrossCoupling
                           ? 1.5
                           : mParameters.largeStepFinalTime},
        TimeDomainCase{"near_transition_unstable",
                       mParameters.stableGainScale,
                       mEnableCurrentCrossCoupling ? 45e-6 : 80e-6,
                       mEnableCurrentCrossCoupling
                           ? 1.5
                           : mParameters.largeStepFinalTime},
        TimeDomainCase{"large_step", mParameters.stableGainScale,
                       mEnableCurrentCrossCoupling ? 80e-6 : 100e-6,
                       mParameters.largeStepFinalTime}};
    // The weak-grid cross-coupled base case requires a smaller perturbation
    // to remain in the local linear regime. Keep the established excitation
    // of the non-cross-coupled validation unchanged.
    const Real perturbationRelative =
        mEnableCurrentCrossCoupling ? 2e-5
                                    : mParameters.gridVoltagePulseRelative;
    for (const auto &timeDomainCase : timeDomainCases) {
      for (const ModelKind model : models) {
        const TimeDomainResult result = runTimeDomainCase(
            model, powerFlow, op, timeDomainCase.gainScale,
            timeDomainCase.timeStep, timeDomainCase.finalTime,
            timeDomainCase.name, perturbationRelative);
        timeDomainRecords.push_back(
            {modelName(model), timeDomainCase.name,
             timeDomainCase.gainScale, timeDomainCase.timeStep,
             timeDomainCase.finalTime, mParameters.perturbationTime,
             mParameters.perturbationDuration,
             perturbationRelative,
             result.prePerturbationActivePowerError,
             result.prePerturbationReactivePowerError, result.logPath});
        writeTimeDomainManifest(timeDomainRecords);
      }
    }
  }

  writeEigenvalues(eigenvalueRecords);
  writeSummary(summaryRecords);
  writeParticipation(participationRecords);
  if (mRunTimeDomain)
    writeTimeDomainManifest(timeDomainRecords);
  std::cout << "Results written to " << mOutputDirectory.string() << "\n";
  if (!mRunTimeDomain)
    std::cout
        << "Use --time-domain to additionally generate Study 3 traces.\n";
}

void EMTDPPh3GFLStateSpaceValidation::runBenchmarkSelectionScan() {
  std::filesystem::create_directories(mOutputDirectory);
  std::ofstream stream(mOutputDirectory / "benchmark_selection_scan.csv");
  std::ofstream summaryStream(mOutputDirectory /
                              "benchmark_selection_summary.csv");
  std::ofstream participationStream(
      mOutputDirectory / "benchmark_selection_participation.csv");
  stream << std::setprecision(std::numeric_limits<Real>::max_digits10)
         << "sweep_parameter,sweep_value,sweep_unit,mode_index,"
            "lambda_real,lambda_imag,frequency_hz,damping_ratio\n";
  summaryStream << std::setprecision(std::numeric_limits<Real>::max_digits10)
                << "sweep_parameter,sweep_value,sweep_unit,status,"
                   "max_real_lambda,critical_oscillatory_real,"
                   "critical_oscillatory_frequency_hz\n";
  participationStream
      << std::setprecision(std::numeric_limits<Real>::max_digits10)
      << "sweep_parameter,sweep_value,sweep_unit,mode_index,state_index,"
         "state_name,p_abs,p_normalized\n";

  const Parameters baseline = mParameters;
  const auto evaluate = [this, &stream, &summaryStream, &participationStream](
                            const String &parameter, Real value,
                            const String &unit, Parameters candidate) {
    candidate.updateDerived();
    mParameters = candidate;
    try {
      const OperatingPoint op = operatingPoint();
      const Real gainScale = mEnableCurrentCrossCoupling
                                 ? mParameters.stableGainScale
                                 : 1.0;
      const ReferenceResult reference = buildReferences(
          op, gainScale, 1e-6, mEnableCurrentCrossCoupling);
      Eigen::EigenSolver<Matrix> modalSolver(reference.noDelayA, true);
      if (modalSolver.info() != Eigen::Success)
        throw std::runtime_error("Candidate modal decomposition failed.");
      const VectorComp modes = modalSolver.eigenvalues();
      Real maxReal = -std::numeric_limits<Real>::infinity();
      Real criticalReal = -std::numeric_limits<Real>::infinity();
      Real criticalFrequency = std::numeric_limits<Real>::quiet_NaN();
      Eigen::Index criticalModeIndex = -1;
      for (Eigen::Index idx = 0; idx < modes.rows(); ++idx) {
        const Complex mode = modes(idx);
        if (!std::isfinite(mode.real()) || !std::isfinite(mode.imag()))
          continue;
        maxReal = std::max(maxReal, mode.real());
        if (mode.imag() > 2.0 * PI * 0.1 && mode.real() > criticalReal) {
          criticalReal = mode.real();
          criticalFrequency = mode.imag() / (2.0 * PI);
          criticalModeIndex = idx;
        }
        if (mode.imag() < -DOUBLE_EPSILON)
          continue;
        const Real magnitude = std::abs(mode);
        const Real dampingRatio =
            magnitude > 0.0 ? -mode.real() / magnitude
                            : std::numeric_limits<Real>::quiet_NaN();
        stream << parameter << ',' << value << ',' << unit << ',' << idx
               << ',' << mode.real() << ',' << mode.imag() << ','
               << std::abs(mode.imag()) / (2.0 * PI) << ',' << dampingRatio
               << '\n';
      }
      const String status = maxReal < 0.0 ? "stable" : "unstable";
      summaryStream << parameter << ',' << value << ',' << unit << ','
                    << status << ',' << maxReal << ',' << criticalReal << ','
                    << criticalFrequency << '\n';
      if (parameter == "short_circuit_ratio" && criticalModeIndex >= 0) {
        const std::vector<String> stateNames = {
            "psi",       "phi_pll", "p_filtered", "q_filtered",
            "phi_d",     "phi_q",   "gamma_d",    "gamma_q",
            "vc_d",      "vc_q",    "vc_0",       "if_d",
            "if_q",      "if_0",    "i_line_d",   "i_line_q",
            "i_line_0"};
        const MatrixComp right = modalSolver.eigenvectors();
        const MatrixComp left = right.fullPivLu().inverse();
        const MatrixComp factors =
            Math::elementwiseProduct(right, left.transpose());
        Real total = 0.0;
        for (Eigen::Index state = 0; state < factors.rows(); ++state)
          total += std::abs(factors(state, criticalModeIndex));
        for (Eigen::Index state = 0; state < factors.rows(); ++state) {
          const Real magnitude = std::abs(factors(state, criticalModeIndex));
          const String name = state < static_cast<Eigen::Index>(stateNames.size())
                                  ? stateNames[state]
                                  : "x" + std::to_string(state);
          participationStream << parameter << ',' << value << ',' << unit
                              << ',' << criticalModeIndex << ',' << state
                              << ',' << name << ',' << magnitude << ','
                              << (total > 0.0 ? magnitude / total : 0.0)
                              << '\n';
        }
      }
      std::cout << std::setw(26) << parameter << " = " << std::setw(8)
                << value << ' ' << std::setw(2) << unit << ": " << status
                << ", max Re(lambda) = " << maxReal
                << " 1/s, critical f = " << criticalFrequency << " Hz\n";
    } catch (const std::exception &error) {
      const Real nan = std::numeric_limits<Real>::quiet_NaN();
      summaryStream << parameter << ',' << value << ',' << unit
                    << ",invalid," << nan << ',' << nan << ',' << nan << '\n';
      std::cout << std::setw(26) << parameter << " = " << value << ' '
                << unit << ": invalid (" << error.what() << ")\n";
    }
  };

  for (const Real scr : {1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0,
                         5.0, 7.5, 10.0, 20.0, 33.0}) {
    Parameters candidate = baseline;
    candidate.shortCircuitRatio = scr;
    evaluate("short_circuit_ratio", scr, "-", candidate);
  }
  for (const Real activePowerPu : {0.2, 0.4, 0.6, 0.7, 0.8, 0.9}) {
    Parameters candidate = baseline;
    candidate.activePowerReferencePu = activePowerPu;
    evaluate("active_power_reference", activePowerPu, "pu", candidate);
  }
  for (const Real multiplier : {0.1, 0.25, 0.5, 0.75, 1.0, 1.25}) {
    Parameters candidate = baseline;
    candidate.kpPowerPu = multiplier * baseline.kpPowerPu;
    evaluate("kp_power_multiplier", multiplier, "-", candidate);
  }
  for (const Real multiplier : {0.1, 0.25, 0.5, 0.75, 1.0, 1.25}) {
    Parameters candidate = baseline;
    candidate.kiPowerPuPerSecond =
        multiplier * baseline.kiPowerPuPerSecond;
    evaluate("ki_power_multiplier", multiplier, "-", candidate);
  }
  for (const Real multiplier : {0.1, 0.25, 0.5, 0.75, 1.0, 1.25}) {
    Parameters candidate = baseline;
    candidate.kpCurrentPu = multiplier * baseline.kpCurrentPu;
    evaluate("kp_current_multiplier", multiplier, "-", candidate);
  }
  for (const Real multiplier : {0.1, 0.25, 0.5, 0.75, 1.0, 1.25}) {
    Parameters candidate = baseline;
    candidate.kiCurrentPuPerSecond =
        multiplier * baseline.kiCurrentPuPerSecond;
    evaluate("ki_current_multiplier", multiplier, "-", candidate);
  }
  for (const Real multiplier : {0.1, 0.25, 0.5, 0.75, 1.0, 1.25}) {
    Parameters candidate = baseline;
    candidate.kpPllPuPerSecond = multiplier * baseline.kpPllPuPerSecond;
    evaluate("kp_pll_multiplier", multiplier, "-", candidate);
  }
  for (const Real multiplier : {0.1, 0.25, 0.5, 0.75, 1.0, 1.25}) {
    Parameters candidate = baseline;
    candidate.kiPllPuPerSecondSquared =
        multiplier * baseline.kiPllPuPerSecondSquared;
    evaluate("ki_pll_multiplier", multiplier, "-", candidate);
  }
  for (const Real multiplier : {0.5, 1.0, 2.0, 3.0, 4.0, 6.0}) {
    Parameters candidate = baseline;
    candidate.filterInductiveReactancePu =
        multiplier * baseline.filterInductiveReactancePu;
    evaluate("filter_inductance_multiplier", multiplier, "-", candidate);
  }
  for (const Real multiplier : {0.5, 1.0, 2.0, 3.0, 4.0, 6.0}) {
    Parameters candidate = baseline;
    candidate.filterCapacitiveSusceptancePu =
        multiplier * baseline.filterCapacitiveSusceptancePu;
    evaluate("filter_capacitance_multiplier", multiplier, "-", candidate);
  }
  for (const Real multiplier : {1.0, 1.5, 2.0, 3.0, 4.0}) {
    Parameters candidate = baseline;
    candidate.filterInductiveReactancePu =
        multiplier * baseline.filterInductiveReactancePu;
    candidate.filterCapacitiveSusceptancePu =
        multiplier * baseline.filterCapacitiveSusceptancePu;
    evaluate("combined_filter_multiplier", multiplier, "-", candidate);
  }

  // Candidate fixed controller: only the outer power-loop proportional gain
  // is reduced. These sweeps determine whether grid strength can serve as the
  // single physically interpretable Study 3 continuation parameter.
  Parameters selectedController = baseline;
  selectedController.kpPowerPu = 0.25 * baseline.kpPowerPu;
  for (const Real scr : {1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.5, 10.0}) {
    Parameters candidate = selectedController;
    candidate.shortCircuitRatio = scr;
    evaluate("selected_controller_scr", scr, "-", candidate);
  }
  for (const Real activePowerPu : {0.2, 0.4, 0.6, 0.7, 0.8, 0.9}) {
    Parameters candidate = selectedController;
    candidate.activePowerReferencePu = activePowerPu;
    evaluate("selected_controller_active_power", activePowerPu, "pu",
             candidate);
  }
  for (const Real multiplier : {0.25, 0.30, 0.35, 0.40, 0.425, 0.45, 0.475,
                                0.50}) {
    Parameters candidate = baseline;
    candidate.kpPowerPu = multiplier * baseline.kpPowerPu;
    evaluate("selected_kp_power_multiplier", multiplier, "-", candidate);
  }

  mParameters = baseline;
  writeParameters(mOutputDirectory / "parameters_candidate.csv", mParameters);
  std::cout << "Benchmark-selection modes written to "
            << (mOutputDirectory / "benchmark_selection_scan.csv").string()
            << "\n";
}

void EMTDPPh3GFLStateSpaceValidation::runWeakGridStudy() {
  if (!mEnableCurrentCrossCoupling) {
    std::cout << "The weak-grid extension is defined for the cross-coupled "
                 "benchmark only.\n";
    return;
  }

  std::filesystem::create_directories(mOutputDirectory);
  std::ofstream modeStream(mOutputDirectory / "weak_grid_eigenvalues.csv");
  std::ofstream participationStream(mOutputDirectory /
                                    "weak_grid_participation.csv");
  struct WeakGridTimeDomainRecord {
    Real shortCircuitRatio;
    TimeDomainRecord record;
  };
  std::vector<WeakGridTimeDomainRecord> timeDomainRecords;
  modeStream << std::setprecision(std::numeric_limits<Real>::max_digits10)
             << "short_circuit_ratio,model,method,index,z_real,z_imag,"
                "lambda_real,lambda_imag,frequency_hz\n";
  participationStream
      << std::setprecision(std::numeric_limits<Real>::max_digits10)
      << "short_circuit_ratio,model,mode_index,state_index,state_name,p_abs,"
         "p_normalized\n";

  const Parameters baseline = mParameters;
  constexpr Real timeStep = 1e-6;
  const std::array<ModelKind, 4> models = {
      ModelKind::EmtVariable, ModelKind::EmtSplit, ModelKind::DpVariable,
      ModelKind::DpSplit};
  const std::vector<String> noDelayNames = {
      "psi",       "phi_pll", "p_filtered", "q_filtered",
      "phi_d",     "phi_q",   "gamma_d",    "gamma_q",
      "vc_d",      "vc_q",    "vc_0",       "if_d",
      "if_q",      "if_0",    "i_line_d",   "i_line_q",
      "i_line_0"};

  const auto criticalMode = [](const VectorComp &continuous) {
    Eigen::Index selected = -1;
    Real maximum = -std::numeric_limits<Real>::infinity();
    for (Eigen::Index idx = 0; idx < continuous.rows(); ++idx) {
      const Complex mode = continuous(idx);
      if (mode.imag() > 2.0 * PI * 0.1 && std::isfinite(mode.real()) &&
          mode.real() > maximum) {
        maximum = mode.real();
        selected = idx;
      }
    }
    return selected;
  };

  const auto appendModes = [&modeStream](Real scr, const String &model,
                                         const String &method,
                                         const VectorComp &discrete,
                                         const VectorComp &continuous) {
    for (Eigen::Index idx = 0; idx < continuous.rows(); ++idx) {
      const Complex z = idx < discrete.rows()
                            ? discrete(idx)
                            : Complex(std::numeric_limits<Real>::quiet_NaN(),
                                      std::numeric_limits<Real>::quiet_NaN());
      const Complex mode = continuous(idx);
      modeStream << scr << ',' << model << ',' << method << ',' << idx << ','
                 << z.real() << ',' << z.imag() << ',' << mode.real() << ','
                 << mode.imag() << ',' << std::abs(mode.imag()) / (2.0 * PI)
                 << '\n';
    }
  };

  const auto appendParticipation = [&participationStream](
                                       Real scr, const String &model,
                                       Eigen::Index mode,
                                       const MatrixComp &factors,
                                       const std::vector<String> &names) {
    if (mode < 0 || mode >= factors.cols())
      return;
    Real total = 0.0;
    for (Eigen::Index state = 0; state < factors.rows(); ++state)
      total += std::abs(factors(state, mode));
    for (Eigen::Index state = 0; state < factors.rows(); ++state) {
      const Real magnitude = std::abs(factors(state, mode));
      const String name = state < static_cast<Eigen::Index>(names.size())
                              ? names[state]
                              : "x" + std::to_string(state);
      participationStream << scr << ',' << model << ',' << mode << ','
                          << state << ',' << name << ',' << magnitude << ','
                          << (total > 0.0 ? magnitude / total : 0.0) << '\n';
    }
  };

  const auto appendReference = [&](Real scr, const String &model,
                                   const String &method, const Matrix &matrix,
                                   Bool logarithmic,
                                   const std::vector<String> &names) {
    Eigen::EigenSolver<Matrix> solver(matrix, true);
    if (solver.info() != Eigen::Success)
      throw std::runtime_error("Weak-grid reference decomposition failed.");
    const VectorComp discrete = solver.eigenvalues();
    const VectorComp continuous =
        logarithmic ? logarithmicContinuousEigenvalues(discrete, timeStep)
                    : bilinearContinuousEigenvalues(discrete, timeStep);
    appendModes(scr, model, method, discrete, continuous);
    const MatrixComp right = solver.eigenvectors();
    const MatrixComp left = right.fullPivLu().inverse();
    appendParticipation(
        scr, model, criticalMode(continuous),
        Math::elementwiseProduct(right, left.transpose()), names);
  };

  for (const Real scr : {1.3, 1.325, 1.35, 1.375, 1.4, 1.45, 1.5, 1.6,
                         1.8, 2.0}) {
    mParameters = baseline;
    mParameters.shortCircuitRatio = scr;
    mParameters.updateDerived();
    try {
      const OperatingPoint op = operatingPoint();
      const SystemTopology powerFlow =
          runPowerFlow("GFLValidation_WeakGrid_SCR_" +
                       std::to_string(static_cast<Int>(1000.0 * scr)));
      const ReferenceResult reference = buildReferences(
          op, mParameters.stableGainScale, timeStep, true);
      appendReference(scr, "analytical no delay", "trapezoidal",
                      reference.noDelayAd, false, noDelayNames);
      std::vector<String> splitNames = noDelayNames;
      splitNames.push_back("v_inv_delay_d");
      splitNames.push_back("v_inv_delay_q");
      splitNames.push_back("v_inv_delay_0");
      appendReference(scr, "exact split companion",
                      "partitioned_trapezoidal", reference.splitAd, true,
                      splitNames);

      for (const ModelKind model : models) {
        const ModalResult result = extractOneStep(
            model, powerFlow, op, mParameters.stableGainScale, timeStep,
            "weak_grid_scr_sweep", false, true);
        appendModes(scr, modelName(model), "one_step",
                    result.discreteEigenvalues, result.continuousEigenvalues);
        appendParticipation(scr, modelName(model),
                            criticalMode(result.continuousEigenvalues),
                            result.participationFactors, result.stateNames);
      }
      if (mRunTimeDomain &&
          (std::abs(scr - 1.325) < 1e-12 ||
           std::abs(scr - 1.35) < 1e-12)) {
        const String stabilityCase = scr < 1.34 ? "scr_unstable" : "scr_stable";
        constexpr Real responseTimeStep = 5e-6;
        constexpr Real finalTime = 1.5;
        for (const ModelKind model : models) {
          const TimeDomainResult result = runTimeDomainCase(
              model, powerFlow, op, mParameters.stableGainScale,
              responseTimeStep, finalTime, stabilityCase, 2e-5);
          timeDomainRecords.push_back(
              {scr,
               {modelName(model), stabilityCase,
                mParameters.stableGainScale, responseTimeStep, finalTime,
                mParameters.perturbationTime, mParameters.perturbationDuration,
                2e-5,
                result.prePerturbationActivePowerError,
                result.prePerturbationReactivePowerError, result.logPath}});
        }
      }
      std::cout << "Weak-grid extracted sweep: SCR = " << scr << '\n';
    } catch (const std::exception &error) {
      std::cerr << "Weak-grid SCR " << scr << " skipped: " << error.what()
                << '\n';
    }
  }
  mParameters = baseline;
  std::ofstream timeDomainStream(
      mOutputDirectory / "weak_grid_time_domain_manifest.csv");
  timeDomainStream
      << std::setprecision(std::numeric_limits<Real>::max_digits10)
      << "short_circuit_ratio,model,stability_case,kp_power_multiplier,"
         "time_step_s,time_step_us,final_time_s,perturbation_time_s,"
         "perturbation_duration_s,perturbation_relative,"
         "pre_perturbation_p_error,pre_perturbation_q_error,log_path\n";
  for (const auto &item : timeDomainRecords) {
    const auto &record = item.record;
    timeDomainStream
        << item.shortCircuitRatio << ',' << record.model << ','
        << record.stabilityCase << ',' << record.gainScale << ','
        << record.timeStep << ',' << 1e6 * record.timeStep << ','
        << record.finalTime << ',' << record.perturbationTime << ','
        << record.perturbationDuration << ',' << record.perturbationRelative
        << ',' << record.prePerturbationActivePowerError << ','
        << record.prePerturbationReactivePowerError << ',' << record.logPath
        << '\n';
  }
  std::cout << "Weak-grid modes and participation written to "
            << mOutputDirectory.string() << '\n';
}

void EMTDPPh3GFLStateSpaceValidation::runCandidateCheck() {
  std::filesystem::create_directories(mOutputDirectory);
  std::ofstream stream(mOutputDirectory / "benchmark_candidate_check.csv");
  stream << std::setprecision(std::numeric_limits<Real>::max_digits10)
         << "stability_case,model,kp_power_multiplier,max_real_lambda,"
            "critical_oscillatory_real,critical_frequency_hz\n";

  const OperatingPoint op = operatingPoint();
  const SystemTopology powerFlow =
      runPowerFlow("GFLValidation_CandidateCheck_PowerFlow");
  const std::array<ModelKind, 4> models = {
      ModelKind::EmtVariable, ModelKind::EmtSplit, ModelKind::DpVariable,
      ModelKind::DpSplit};
  constexpr Real timeStep = 1e-6;

  const auto append = [&stream](const String &stabilityCase,
                                const String &model, Real multiplier,
                                const VectorComp &modes) {
    Real maxReal = -std::numeric_limits<Real>::infinity();
    Real criticalReal = -std::numeric_limits<Real>::infinity();
    Real criticalFrequency = std::numeric_limits<Real>::quiet_NaN();
    for (Eigen::Index idx = 0; idx < modes.rows(); ++idx) {
      const Complex mode = modes(idx);
      if (!std::isfinite(mode.real()) || !std::isfinite(mode.imag()))
        continue;
      maxReal = std::max(maxReal, mode.real());
      if (mode.imag() > 2.0 * PI * 0.1 && mode.real() > criticalReal) {
        criticalReal = mode.real();
        criticalFrequency = mode.imag() / (2.0 * PI);
      }
    }
    stream << stabilityCase << ',' << model << ',' << multiplier << ','
           << maxReal << ',' << criticalReal << ',' << criticalFrequency
           << '\n';
    std::cout << std::setw(10) << stabilityCase << " | " << std::setw(22)
              << model << " | max Re(lambda) = " << std::setw(12) << maxReal
              << " 1/s | critical = " << criticalReal << " 1/s at "
              << criticalFrequency << " Hz\n";
  };

  for (const auto &[stabilityCase, multiplier] :
       std::array<std::pair<String, Real>, 2>{
           std::pair<String, Real>{"stable", mParameters.stableGainScale},
           std::pair<String, Real>{"unstable",
                                   mParameters.unstableGainScale}}) {
    const ReferenceResult reference = buildReferences(
        op, multiplier, timeStep, mEnableCurrentCrossCoupling);
    append(stabilityCase, "analytical no delay", multiplier,
           eigenvalues(reference.noDelayA));
    append(stabilityCase, "exact split companion", multiplier,
           logarithmicContinuousEigenvalues(eigenvalues(reference.splitAd),
                                            timeStep));
    for (const ModelKind model : models) {
      const ModalResult result = extractOneStep(
          model, powerFlow, op, multiplier, timeStep, "candidate_check",
          false, mEnableCurrentCrossCoupling);
      append(stabilityCase, modelName(model), multiplier,
             result.continuousEigenvalues);
    }
  }
  std::cout << "Candidate check written to "
            << (mOutputDirectory / "benchmark_candidate_check.csv").string()
            << "\n";
}

void EMTDPPh3GFLStateSpaceValidation::runEmtVariableDiagnostics() {
  std::filesystem::create_directories(mOutputDirectory);
  const OperatingPoint op = operatingPoint();
  const SystemTopology powerFlow =
      runPowerFlow("GFLValidation_EmtVariableDiagnostics_PowerFlow");

  std::ofstream modalStream(mOutputDirectory /
                            "emt_variable_modal_diagnostics.csv");
  modalStream << std::setprecision(std::numeric_limits<Real>::max_digits10)
              << "diagnostic,method,time_step_s,time_step_us,warmup_periods,"
                 "phase_deg,index,z_real,z_imag,lambda_real,lambda_imag\n";

  const auto appendModes = [&modalStream](
                               const String &diagnostic, const String &method,
                               Real timeStep, UInt warmupPeriods, Real phaseDeg,
                               const VectorComp &discrete,
                               const VectorComp &continuous) {
    for (Eigen::Index idx = 0; idx < continuous.rows(); ++idx) {
      const Complex z = idx < discrete.rows()
                            ? discrete(idx)
                            : Complex(std::numeric_limits<Real>::quiet_NaN(),
                                      std::numeric_limits<Real>::quiet_NaN());
      modalStream << diagnostic << ',' << method << ',' << timeStep << ','
                  << 1e6 * timeStep << ',' << warmupPeriods << ',' << phaseDeg
                  << ',' << idx << ',' << z.real() << ',' << z.imag() << ','
                  << continuous(idx).real() << ',' << continuous(idx).imag()
                  << '\n';
    }
  };

  const std::array<Real, 11> convergenceSteps = {
      0.1e-6, 0.125e-6, 0.25e-6, 0.5e-6, 1e-6,  2e-6,
      5e-6,   10e-6,   20e-6,   50e-6,  100e-6};
  for (const Real timeStep : convergenceSteps) {
    const ReferenceResult reference =
        buildReferences(op, mParameters.stableGainScale, timeStep);
    const VectorComp referenceZ = eigenvalues(reference.noDelayAd);
    appendModes("timestep_convergence", "analytical_no_delay", timeStep, 0,
                0.0, referenceZ,
                bilinearContinuousEigenvalues(referenceZ, timeStep));

    // Control experiment: express the same continuous no-delay reference in
    // native abc coordinates, freeze its time-periodic matrix at the beginning
    // of the step, apply trapezoidal integration, and transform the result back
    // to global dq0. This isolates the error caused solely by freezing a
    // rotating-frame Jacobian over one EMT step.
    const Matrix frozenAbcAd = frozenAbcTrapezoidalStep(
        reference.noDelayA, mParameters.omega, 0.0, timeStep);
    const VectorComp frozenAbcZ = eigenvalues(frozenAbcAd);
    appendModes("timestep_convergence", "analytical_frozen_abc_one_step",
                timeStep, 0, 0.0, frozenAbcZ,
                bilinearContinuousEigenvalues(frozenAbcZ, timeStep));
    if (timeStep <= 10e-6) {
      const UInt periodSteps = static_cast<UInt>(std::llround(
          1.0 / (mParameters.frequency * timeStep)));
      VectorComp frozenMultipliers(frozenAbcZ.rows());
      for (Eigen::Index idx = 0; idx < frozenAbcZ.rows(); ++idx)
        frozenMultipliers(idx) =
            std::pow(frozenAbcZ(idx), static_cast<Real>(periodSteps));
      appendModes("timestep_convergence", "analytical_frozen_abc_monodromy",
                  timeStep, 0, 0.0, frozenMultipliers,
                  logarithmicContinuousEigenvalues(
                      frozenMultipliers, 1.0 / mParameters.frequency));
    }

    for (const auto &[model, prefix] :
         std::array<std::pair<ModelKind, String>, 2>{
             std::pair<ModelKind, String>{ModelKind::EmtVariable, ""},
             std::pair<ModelKind, String>{ModelKind::EmtLegacyVariable,
                                          "legacy_"}}) {
      const ModalResult oneStep = extractOneStep(
          model, powerFlow, op, mParameters.stableGainScale, timeStep,
          "emt_variable_compact_diagnostics", true);
      appendModes("timestep_convergence",
                  prefix + "compact_history_one_step", timeStep, 0,
                  0.0, oneStep.discreteEigenvalues,
                  oneStep.continuousEigenvalues);

      if (timeStep <= 10e-6) {
        const VectorComp multipliers = calculateMonodromy(
            model, powerFlow, op, mParameters.stableGainScale, timeStep, 0,
            true);
        appendModes("timestep_convergence",
                    prefix + "compact_history_monodromy_initial",
                    timeStep, 0, 0.0, multipliers,
                    logarithmicContinuousEigenvalues(
                        multipliers, 1.0 / mParameters.frequency));
      }
      if (timeStep == 1e-6 || timeStep == 10e-6) {
        const VectorComp multipliers = calculateMonodromy(
            model, powerFlow, op, mParameters.stableGainScale, timeStep, 2,
            true);
        appendModes("timestep_convergence",
                    prefix + "compact_history_monodromy_warm",
                    timeStep, 2, 0.0, multipliers,
                    logarithmicContinuousEigenvalues(
                        multipliers, 1.0 / mParameters.frequency));
      }
    }

    const ModalResult augmented = extractOneStep(
        ModelKind::EmtVariable, powerFlow, op, mParameters.stableGainScale,
        timeStep, "emt_variable_augmented_diagnostics");
    appendModes("timestep_convergence", "augmented_physical_one_step",
                timeStep, 0, 0.0, augmented.discreteEigenvalues,
                augmented.continuousEigenvalues);
    if (timeStep <= 10e-6) {
      const VectorComp multipliers = calculateMonodromy(
          ModelKind::EmtVariable, powerFlow, op, mParameters.stableGainScale,
          timeStep);
      appendModes("timestep_convergence", "augmented_physical_monodromy",
                  timeStep, 0, 0.0, multipliers,
                  logarithmicContinuousEigenvalues(
                      multipliers, 1.0 / mParameters.frequency));
    }
  }

  // Differentiate the actual nonlinear EMT simulation step in the same
  // augmented coordinates used by the diagnostic contributor:
  //
  //   [grid-inductor history, inverter physical state, inverter u_previous].
  //
  // This is an end-to-end check: every perturbed column is obtained from a
  // fresh initialized simulation, including pre-step relinearization, MNA
  // solution and component post-step updates.
  constexpr Real mapValidationTimeStep = 1e-6;
  struct OneStepSample {
    Matrix initial;
    Matrix next;
    Matrix extractedAd;
  };
  const auto runAugmentedOneStep =
      [this, &powerFlow, &op](Int perturbIndex, Real perturbation,
                             UInt runIndex) {
        SystemHandles handles = buildSystem(
            ModelKind::EmtVariable, powerFlow, op,
            mParameters.stableGainScale, nullptr);
        handles.emtVariableBase->useAugmentedPhysicalStateExtraction();

        Simulation simulation(
            "GFLValidation_AugmentedMapFD_" + std::to_string(runIndex),
            Logger::Level::warn);
        simulation.setSystem(handles.system);
        simulation.setDomain(Domain::EMT);
        simulation.setSolverType(Solver::Type::MNA);
        simulation.doSystemMatrixRecomputation(true);
        simulation.doInitFromNodesAndTerminals(true);
        simulation.doStateSpaceExtraction(true);
        simulation.setTimeStep(mapValidationTimeStep);
        simulation.setFinalTime(mapValidationTimeStep);
        simulation.initialize();

        auto inverterState =
            handles.emtVariableBase->attributeTyped<Matrix>("x");
        auto inverterVoltage =
            handles.emtVariableBase->attributeTyped<Matrix>("v_intf");
        auto inductorVoltage =
            handles.emtGridInductor->attributeTyped<Matrix>("v_intf");
        auto inductorCurrent =
            handles.emtGridInductor->attributeTyped<Matrix>("i_intf");
        const Matrix &gridConductance =
            handles.emtGridInductor->getMNAConductance();

        const UInt inverterStateCount =
            handles.emtVariableBase->getStateCount();
        Matrix coordinate = Matrix::Zero(3 + inverterStateCount + 3, 1);
        coordinate.block(0, 0, 3, 1) =
            gridConductance * (**inductorVoltage) + (**inductorCurrent);
        coordinate.block(3, 0, inverterStateCount, 1) = **inverterState;
        coordinate.block(3 + inverterStateCount, 0, 3, 1) =
            **inverterVoltage;

        if (perturbIndex >= 0) {
          coordinate(perturbIndex, 0) += perturbation;
          if (perturbIndex < 3) {
            inductorCurrent->set(coordinate.block(0, 0, 3, 1) -
                                 gridConductance * (**inductorVoltage));
          } else if (perturbIndex <
                     static_cast<Int>(3 + inverterStateCount)) {
            inverterState->set(
                coordinate.block(3, 0, inverterStateCount, 1));
          } else {
            inverterVoltage->set(
                coordinate.block(3 + inverterStateCount, 0, 3, 1));
          }
        }

        simulation.start();
        simulation.step();

        Matrix next = Matrix::Zero(coordinate.rows(), 1);
        next.block(0, 0, 3, 1) =
            gridConductance * (**inductorVoltage) + (**inductorCurrent);
        next.block(3, 0, inverterStateCount, 1) = **inverterState;
        next.block(3 + inverterStateCount, 0, 3, 1) =
            **inverterVoltage;
        const Matrix extractedAd =
            simulation.getStateSpaceExtractor().getDiscreteStateMatrix();
        simulation.stop();
        return OneStepSample{coordinate, next, extractedAd};
      };

  const OneStepSample nominalMap = runAugmentedOneStep(-1, 0.0, 0);
  std::ofstream mapStream(
      mOutputDirectory / "emt_variable_one_step_map_validation.csv");
  mapStream << std::setprecision(std::numeric_limits<Real>::max_digits10)
            << "relative_step,max_absolute_error,max_relative_column_error,"
               "relative_frobenius_error\n";
  UInt mapRunIndex = 1;
  for (const Real relativeStep : {1e-4, 3e-5, 1e-5, 3e-6, 1e-6}) {
    Matrix numericalAd = Matrix::Zero(nominalMap.initial.rows(),
                                      nominalMap.initial.rows());
    for (Eigen::Index column = 0; column < nominalMap.initial.rows();
         ++column) {
      const Real step = relativeStep *
                        std::max(1.0, std::abs(nominalMap.initial(column, 0)));
      const OneStepSample plus =
          runAugmentedOneStep(static_cast<Int>(column), step, mapRunIndex++);
      const OneStepSample minus =
          runAugmentedOneStep(static_cast<Int>(column), -step, mapRunIndex++);
      numericalAd.col(column) = (plus.next - minus.next) / (2.0 * step);
    }

    const Matrix error = numericalAd - nominalMap.extractedAd;
    Real maxRelativeColumnError = 0.0;
    for (Eigen::Index column = 0; column < numericalAd.cols(); ++column) {
      maxRelativeColumnError =
          std::max(maxRelativeColumnError,
                   error.col(column).norm() /
                       std::max(1.0, numericalAd.col(column).norm()));
    }
    mapStream << relativeStep << ',' << error.cwiseAbs().maxCoeff() << ','
              << maxRelativeColumnError << ','
              << error.norm() / std::max(1.0, numericalAd.norm()) << '\n';
    if (relativeStep == 1e-5) {
      Matrix transformNow = Matrix::Identity(numericalAd.rows(),
                                             numericalAd.cols());
      Matrix transformNext = transformNow;
      const Real thetaNow = mParameters.omega * mapValidationTimeStep;
      const Real thetaNext =
          mParameters.omega * (2.0 * mapValidationTimeStep);
      for (const UInt offset : {0u, 11u, 14u, 17u}) {
        transformNow.block(offset, offset, 3, 3) =
            parkTransformDQ0(thetaNow);
        transformNext.block(offset, offset, 3, 3) =
            parkTransformDQ0(thetaNext);
      }
      const Matrix numericalGlobalDq0 =
          transformNext * numericalAd * transformNow.transpose();
      const VectorComp numericalZ = eigenvalues(numericalGlobalDq0);
      appendModes("one_step_map_validation",
                  "finite_difference_actual_one_step",
                  mapValidationTimeStep, 0, 0.0, numericalZ,
                  bilinearContinuousEigenvalues(numericalZ,
                                                mapValidationTimeStep));
    }
  }

  // Inspect the one-step tangent at eight points of the fundamental period,
  // both immediately and after two undisturbed warm-up periods.
  constexpr Real phaseTimeStep = 1e-6;
  const UInt stepsPerPeriod = static_cast<UInt>(
      std::llround(1.0 / (mParameters.frequency * phaseTimeStep)));
  for (const ModelKind model : {ModelKind::EmtVariable,
                                ModelKind::EmtLegacyVariable}) {
    for (const UInt warmupPeriods : {0u, 2u}) {
      SystemHandles handles = buildSystem(
          model, powerFlow, op, mParameters.stableGainScale, nullptr);
      handles.emtVariableBase->useAugmentedPhysicalStateExtraction(false);
      Simulation simulation("GFLValidation_EmtVariablePhase_" +
                                fileToken(model) + "_w" +
                                std::to_string(warmupPeriods),
                            Logger::Level::warn);
    simulation.setSystem(handles.system);
    simulation.setDomain(Domain::EMT);
    simulation.setSolverType(Solver::Type::MNA);
    simulation.doSystemMatrixRecomputation(true);
    simulation.doInitFromNodesAndTerminals(true);
    simulation.doStateSpaceExtraction(true);
    simulation.setTimeStep(phaseTimeStep);
    simulation.setFinalTime((warmupPeriods + 1) * stepsPerPeriod *
                            phaseTimeStep);
    simulation.initialize();
    simulation.start();
    for (UInt step = 1; step <= (warmupPeriods + 1) * stepsPerPeriod; ++step) {
      simulation.step();
      if (step <= warmupPeriods * stepsPerPeriod)
        continue;
      const UInt periodStep = step - warmupPeriods * stepsPerPeriod;
      if (periodStep % (stepsPerPeriod / 8) != 0)
        continue;
      StateSpaceModalAnalysis modal(simulation.getStateSpaceExtractor());
      modal.setAnalysisFrame(StateSpaceAnalysisFrame::GlobalDQ0);
      modal.setGlobalDq0Frame(mParameters.omega);
      modal.update();
      const Real phaseDeg =
          360.0 * static_cast<Real>(periodStep) / stepsPerPeriod;
      appendModes("phase_sweep", model == ModelKind::EmtVariable
                                     ? "compact_history_one_step"
                                     : "legacy_compact_history_one_step",
                  phaseTimeStep, warmupPeriods, phaseDeg,
                  modal.getDiscreteEigenvalues(),
                  modal.getContinuousEigenvalues());
    }
    simulation.stop();
    }
  }

  std::ofstream jacobianStream(mOutputDirectory /
                               "emt_variable_jacobian_validation.csv");
  jacobianStream
      << std::setprecision(std::numeric_limits<Real>::max_digits10)
      << "phase_deg,relative_step,matrix,max_absolute_error,"
         "relative_frobenius_error,max_column_relative_error\n";

  const auto writeJacobianMetric = [&jacobianStream](
                                       Real phaseDeg, Real relativeStep,
                                       const String &name,
                                       const Matrix &analytical,
                                       const Matrix &numerical) {
    const Matrix difference = analytical - numerical;
    const Real maxAbsolute = difference.cwiseAbs().maxCoeff();
    const Real relativeFrobenius =
        difference.norm() /
        std::max(numerical.norm(), std::numeric_limits<Real>::epsilon());
    Real maxColumnRelative = 0.0;
    for (Eigen::Index column = 0; column < numerical.cols(); ++column) {
      const Real denominator = std::max(
          numerical.col(column).norm(), std::numeric_limits<Real>::epsilon());
      maxColumnRelative = std::max(
          maxColumnRelative, difference.col(column).norm() / denominator);
    }
    jacobianStream << phaseDeg << ',' << relativeStep << ',' << name << ','
                   << maxAbsolute << ',' << relativeFrobenius << ','
                   << maxColumnRelative << '\n';
  };

  SystemHandles handles = buildSystem(
      ModelKind::EmtVariable, powerFlow, op, mParameters.stableGainScale,
      nullptr);
  Simulation simulation("GFLValidation_EmtVariableJacobian",
                        Logger::Level::warn);
  simulation.setSystem(handles.system);
  simulation.setDomain(Domain::EMT);
  simulation.setSolverType(Solver::Type::MNA);
  simulation.doSystemMatrixRecomputation(true);
  simulation.doInitFromNodesAndTerminals(true);
  simulation.setTimeStep(phaseTimeStep);
  simulation.setFinalTime(1.0 / mParameters.frequency);
  simulation.initialize();
  if (!handles.emtVariableInverter)
    throw std::logic_error("EMT-variable diagnostic handle is missing.");

  const auto validate = [&](Real phaseDeg, Real relativeStep) {
    const auto result =
        handles.emtVariableInverter->validateAnalyticalJacobians(relativeStep);
    writeJacobianMetric(phaseDeg, relativeStep, "A", result.analyticalA,
                        result.numericalA);
    writeJacobianMetric(phaseDeg, relativeStep, "B", result.analyticalB,
                        result.numericalB);
    writeJacobianMetric(phaseDeg, relativeStep, "C", result.analyticalC,
                        result.numericalC);
    writeJacobianMetric(phaseDeg, relativeStep, "D", result.analyticalD,
                        result.numericalD);
  };
  for (const Real relativeStep : {1e-4, 1e-5, 1e-6, 1e-7, 1e-8})
    validate(0.0, relativeStep);

  simulation.start();
  for (UInt step = 1; step <= stepsPerPeriod; ++step) {
    simulation.step();
    if (step % (stepsPerPeriod / 4) == 0) {
      const Real phaseDeg =
          360.0 * static_cast<Real>(step) / stepsPerPeriod;
      validate(phaseDeg, 1e-6);
    }
  }
  simulation.stop();

  std::cout << "EMT-variable modal diagnostics written to "
            << (mOutputDirectory / "emt_variable_modal_diagnostics.csv")
                   .string()
            << "\nEMT-variable Jacobian diagnostics written to "
            << (mOutputDirectory / "emt_variable_jacobian_validation.csv")
                   .string()
            << '\n';
}

} // namespace

#ifndef DP_SIM_GFL_CROSS_COUPLED_VALIDATION
int main(int argc, char **argv) {
  Bool runTimeDomain = false;
  Bool runBenchmarkScan = false;
  Bool runWeakGridStudy = false;
  Bool runCandidateCheck = false;
  Bool runEmtVariableDiagnostics = false;
  for (Int idx = 1; idx < argc; ++idx) {
    if (String(argv[idx]) == "--time-domain" || String(argv[idx]) == "--all")
      runTimeDomain = true;
    if (String(argv[idx]) == "--benchmark-scan")
      runBenchmarkScan = true;
    if (String(argv[idx]) == "--weak-grid-study")
      runWeakGridStudy = true;
    if (String(argv[idx]) == "--candidate-check")
      runCandidateCheck = true;
    if (String(argv[idx]) == "--emt-variable-diagnostics")
      runEmtVariableDiagnostics = true;
  }
  EMTDPPh3GFLStateSpaceValidation example(runTimeDomain);
  if (runBenchmarkScan)
    example.runBenchmarkSelectionScan();
  else if (runWeakGridStudy)
    example.runWeakGridStudy();
  else if (runCandidateCheck)
    example.runCandidateCheck();
  else if (runEmtVariableDiagnostics)
    example.runEmtVariableDiagnostics();
  else
    example.run();
  return 0;
}
#endif
