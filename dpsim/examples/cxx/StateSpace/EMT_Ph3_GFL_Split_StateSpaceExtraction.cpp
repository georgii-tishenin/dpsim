// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include <DPsim.h>

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
#include <string>
#include <vector>

#include <dpsim-models/EMT/EMT_Ph3_SSN_GFL.h>
#include <dpsim-models/EMT/EMT_Ph3_SSN_GFL_Split.h>

using namespace CPS;
using namespace DPsim;

namespace {

constexpr UInt ControllerStateCount = 8;
constexpr UInt ElectricalStateCount = 9;
constexpr UInt DelayStateCount = 3;
constexpr UInt NoDelayStateCount = ControllerStateCount + ElectricalStateCount;
constexpr UInt DelayedStateCount = NoDelayStateCount + DelayStateCount;

enum ControllerStateIndex : UInt {
  ThetaPLL = 0,
  PhiPLL = 1,
  PFiltered = 2,
  QFiltered = 3,
  PhiD = 4,
  PhiQ = 5,
  GammaD = 6,
  GammaQ = 7,
};

constexpr UInt VcOffset = 0;
constexpr UInt IfOffset = 3;
constexpr UInt IGridOffset = 6;

struct SteadyState {
  Complex sourceVoltage;
  Complex seriesVoltage;
  Complex pccVoltage;
  Complex filterVoltage;
  Complex filterCurrent;
  Complex gridCurrent;
  Complex converterVoltage;
  Matrix controllerState;
  Matrix electricalState;
  Matrix controllerMeasurement;
  Matrix delayedControllerOutput;
};

struct ExtractionResult {
  VectorComp discreteEigenvalues;
  Real extractionTime;
  UInt stateCount;
};

struct EigenvalueRecord {
  Real timeStep;
  String model;
  String representation;
  VectorComp values;
};

VectorComp eigenvalues(const Matrix &matrix) {
  Eigen::EigenSolver<Matrix> solver(matrix);
  if (solver.info() != Eigen::Success)
    throw std::runtime_error("Eigenvalue computation failed.");
  return solver.eigenvalues();
}

Real maxNearestDistance(const VectorComp &reference, const VectorComp &actual) {
  Real maximum = 0.0;
  for (Eigen::Index refIdx = 0; refIdx < reference.rows(); ++refIdx) {
    Real nearest = std::numeric_limits<Real>::max();
    for (Eigen::Index actualIdx = 0; actualIdx < actual.rows(); ++actualIdx)
      nearest =
          std::min(nearest, std::abs(reference(refIdx) - actual(actualIdx)));
    maximum = std::max(maximum, nearest);
  }
  return maximum;
}

Real symmetricEigenvalueDistance(const VectorComp &reference,
                                 const VectorComp &actual) {
  return std::max(maxNearestDistance(reference, actual),
                  maxNearestDistance(actual, reference));
}

void appendEigenvaluesToCsv(std::ofstream &stream,
                            const EigenvalueRecord &record) {
  for (Eigen::Index idx = 0; idx < record.values.rows(); ++idx) {
    stream << record.timeStep << ',' << 1e6 * record.timeStep << ','
           << record.model << ',' << record.representation << ',' << idx << ','
           << record.values(idx).real() << ',' << record.values(idx).imag()
           << '\n';
  }
}

void writeEigenvaluesCsv(const std::filesystem::path &path,
                         const std::vector<EigenvalueRecord> &records) {
  std::filesystem::create_directories(path.parent_path());
  std::ofstream stream(path);
  if (!stream.is_open())
    throw std::runtime_error("Could not open eigenvalue CSV output: " +
                             path.string());
  stream << std::setprecision(std::numeric_limits<Real>::max_digits10);
  stream << "time_step_s,time_step_us,model,representation,index,real,imag\n";
  for (const auto &record : records)
    appendEigenvaluesToCsv(stream, record);
  if (!stream)
    throw std::runtime_error("Could not write eigenvalue CSV output: " +
                             path.string());
}

Matrix dq0Vector(const Complex &value) {
  Matrix result = Matrix::Zero(3, 1);
  result(0, 0) = value.real();
  result(1, 0) = value.imag();
  return result;
}

Matrix localFromGlobalDq0(Real angle) {
  Matrix transform = Matrix::Identity(3, 3);
  const Real cosine = std::cos(angle);
  const Real sine = std::sin(angle);
  transform(0, 0) = cosine;
  transform(0, 1) = sine;
  transform(1, 0) = -sine;
  transform(1, 1) = cosine;
  return transform;
}

Matrix rotatingFrameDerivative(const Matrix &value, Real omega) {
  Matrix result = Matrix::Zero(3, 1);
  result(0, 0) = omega * value(1, 0);
  result(1, 0) = -omega * value(0, 0);
  return result;
}

Matrix numericalJacobian(const std::function<Matrix(const Matrix &)> &function,
                         const Matrix &operatingPoint) {
  const Matrix value = function(operatingPoint);
  Matrix jacobian(value.rows(), operatingPoint.rows());
  for (Eigen::Index column = 0; column < operatingPoint.rows(); ++column) {
    const Real step = 1e-6 * std::max(1.0, std::abs(operatingPoint(column, 0)));
    Matrix plus = operatingPoint;
    Matrix minus = operatingPoint;
    plus(column, 0) += step;
    minus(column, 0) -= step;
    jacobian.col(column) = (function(plus) - function(minus)) / (2.0 * step);
  }
  return jacobian;
}

Matrix trapezoidalStateMatrix(const Matrix &continuousA, Real timeStep) {
  const Matrix identity =
      Matrix::Identity(continuousA.rows(), continuousA.cols());
  const Matrix lhs = identity - 0.5 * timeStep * continuousA;
  const Matrix rhs = identity + 0.5 * timeStep * continuousA;
  return lhs.fullPivLu().solve(rhs);
}

Matrix trapezoidalInputMatrix(const Matrix &continuousA,
                              const Matrix &continuousB, Real timeStep) {
  const Matrix identity =
      Matrix::Identity(continuousA.rows(), continuousA.cols());
  const Matrix lhs = identity - 0.5 * timeStep * continuousA;
  return lhs.fullPivLu().solve(0.5 * timeStep * continuousB);
}

} // namespace

class EMTPh3SplitGFLStateSpaceExtractionExample {
public:
  EMTPh3SplitGFLStateSpaceExtractionExample(Real extractionTime = 0.1,
                                            Real finalTime = 0.1)
      : mFrequency(50.0), mOmega(2.0 * PI * mFrequency),
        mGridVoltageRmsLineToLine(400.0), mGridResistance(0.3),
        mGridInductance(0.1e-3), mLf(2e-3), mCf(10e-6), mRf(0.2), mRc(0.2),
        mKpPLL(0.25), mKiPLL(0.2), mOmegaCutoff(mOmega), mPRef(10000.0),
        mQRef(5000.0), mKpPowerCtrl(0.05), mKiPowerCtrl(0.2), mKpCurrCtrl(0.25),
        mKiCurrCtrl(1.0), mExtractionTime(extractionTime),
        mFinalTime(finalTime) {
    if (!(mExtractionTime > 0.0))
      throw std::invalid_argument("Extraction time must be positive.");
    if (mFinalTime < mExtractionTime)
      throw std::invalid_argument(
          "Final simulation time must not precede extraction time.");
  }

  void run() const {
    const String simName = "EMT_Ph3_GFL_Split_StateSpaceExtraction";
    const SteadyState steadyState = calculateSteadyState();
    const Matrix noDelayA =
        buildIndependentNoDelayContinuousMatrix(steadyState);
    const VectorComp noDelayLambda = eigenvalues(noDelayA);

    const std::array<Real, 4> timeSteps{1e-6, 10e-6, 100e-6, 1e-3};
    const std::array<String, 4> timeStepLabels{"1us", "10us", "100us", "1ms"};
    std::vector<EigenvalueRecord> records;

    std::cout
        << "\n============================================================\n"
        << "EMT Ph3 GFL state-space extraction comparison\n"
        << "============================================================\n"
        << "Topology: ideal source -> R/L grid -> SSN GFL\n"
        << "Analysis frame: global dq0\n"
        << "Reference: independent nonlinear dq0 equations, numerical "
           "linearization, trapezoidal discretization\n"
        << "Operating point: manual steady-state calculation\n"
        << "Time steps: 1 us, 10 us, 100 us, 1 ms\n"
        << "Extraction time: " << mExtractionTime << " s\n"
        << "Final simulation time: " << mFinalTime << " s\n"
        << "\nManual operating point (global dq phasors):\n"
        << "  source voltage    = " << steadyState.sourceVoltage << " V\n"
        << "  PCC voltage       = " << steadyState.pccVoltage << " V\n"
        << "  capacitor voltage = " << steadyState.filterVoltage << " V\n"
        << "  grid current      = " << steadyState.gridCurrent << " A\n"
        << "  filter current    = " << steadyState.filterCurrent << " A\n"
        << "  converter voltage = " << steadyState.converterVoltage << " V\n";

    for (UInt caseIdx = 0; caseIdx < timeSteps.size(); ++caseIdx) {
      const Real timeStep = timeSteps[caseIdx];
      const String suffix = "_dt_" + timeStepLabels[caseIdx];

      auto splitInverter =
          EMT::Ph3::SSN_GFL_Split::make("InverterSplit", Logger::Level::warn);
      setInverterParameters(splitInverter);
      const ExtractionResult extractedSplit =
          runExtraction(steadyState, splitInverter, simName + "_Split" + suffix,
                        timeStep, DelayedStateCount);

      auto noDelayInverter =
          EMT::Ph3::SSN_GFL::make("InverterNoDelay", Logger::Level::warn);
      setInverterParameters(noDelayInverter);
      const ExtractionResult extractedNoDelay = runExtraction(
          steadyState, noDelayInverter, simName + "_NoDelay" + suffix, timeStep,
          NoDelayStateCount);

      const Matrix independentDelayedAd =
          buildIndependentDelayedDiscreteMatrix(steadyState, timeStep);
      const Matrix independentNoDelayAd =
          trapezoidalStateMatrix(noDelayA, timeStep);
      const VectorComp independentDelayedZ = eigenvalues(independentDelayedAd);
      const VectorComp independentNoDelayZ = eigenvalues(independentNoDelayAd);

      records.push_back({timeStep, "extracted_ssn_gfl_split", "discrete_z",
                         extractedSplit.discreteEigenvalues});
      records.push_back({timeStep, "extracted_ssn_gfl", "discrete_z",
                         extractedNoDelay.discreteEigenvalues});
      records.push_back({timeStep, "independent_dq0_split_delay", "discrete_z",
                         independentDelayedZ});
      records.push_back({timeStep, "independent_dq0_no_delay", "discrete_z",
                         independentNoDelayZ});
      records.push_back({timeStep, "independent_dq0_no_delay",
                         "continuous_lambda", noDelayLambda});

      const Real delayedDifference = symmetricEigenvalueDistance(
          independentDelayedZ, extractedSplit.discreteEigenvalues);
      const Real noDelayDifference = symmetricEigenvalueDistance(
          independentNoDelayZ, extractedNoDelay.discreteEigenvalues);

      std::cout << "\nDelta t = " << 1e6 * timeStep << " us\n"
                << "  split states: " << extractedSplit.stateCount
                << " (8 controller + 6 plant + 3 grid + 3 delay)\n"
                << "  no-delay states: " << extractedNoDelay.stateCount
                << " (8 controller + 6 plant + 3 grid)\n"
                << "  split extraction vs independent delayed reference "
                   "max error: "
                << delayedDifference << "\n"
                << "  no-delay extraction vs independent no-delay reference "
                   "max error: "
                << noDelayDifference << "\n";
    }

    const std::filesystem::path eigenvalueCsvPath =
        std::filesystem::path("logs") / simName / "eigenvalues.csv";
    writeEigenvaluesCsv(eigenvalueCsvPath, records);
    std::cout << "\nEigenvalue CSV: " << eigenvalueCsvPath.string() << "\n";
  }

private:
  template <typename InverterType>
  void
  setInverterParameters(const std::shared_ptr<InverterType> &inverter) const {
    inverter->setParameters(mLf, mCf, mRf, mRc, mOmega, mKpPLL, mKiPLL,
                            mOmegaCutoff, mPRef, mQRef, mKpPowerCtrl,
                            mKiPowerCtrl, mKpCurrCtrl, mKiCurrCtrl);
  }

  SteadyState calculateSteadyState() const {
    const Complex imaginaryUnit(0.0, 1.0);
    const Complex sourceVoltage(mGridVoltageRmsLineToLine, 0.0);
    const Complex powerReference(mPRef, mQRef);
    const Complex totalImpedance(mGridResistance + mRc,
                                 mOmega * mGridInductance);

    Complex filterVoltage = sourceVoltage;
    Bool converged = false;
    for (Int iteration = 0; iteration < 100; ++iteration) {
      const Complex gridCurrent = std::conj(powerReference / filterVoltage);
      const Complex nextFilterVoltage =
          sourceVoltage + totalImpedance * gridCurrent;
      if (std::abs(nextFilterVoltage - filterVoltage) < 1e-12) {
        filterVoltage = nextFilterVoltage;
        converged = true;
        break;
      }
      filterVoltage = nextFilterVoltage;
    }
    if (!converged)
      throw std::runtime_error(
          "Manual steady-state calculation did not converge.");

    const Complex gridCurrent = std::conj(powerReference / filterVoltage);
    const Complex pccVoltage = filterVoltage - mRc * gridCurrent;
    const Complex seriesVoltage =
        sourceVoltage + imaginaryUnit * mOmega * mGridInductance * gridCurrent;
    const Complex filterCurrent =
        gridCurrent + imaginaryUnit * mOmega * mCf * filterVoltage;
    const Complex converterVoltage =
        filterVoltage + (mRf + imaginaryUnit * mOmega * mLf) * filterCurrent;

    const Real theta = std::arg(filterVoltage);
    const Complex localRotation = std::exp(-imaginaryUnit * theta);
    const Complex filterVoltageLocal = filterVoltage * localRotation;
    const Complex gridCurrentLocal = gridCurrent * localRotation;
    const Complex converterVoltageLocal = converterVoltage * localRotation;
    const Real pInitial = filterVoltageLocal.real() * gridCurrentLocal.real() +
                          filterVoltageLocal.imag() * gridCurrentLocal.imag();
    const Real qInitial = -filterVoltageLocal.real() * gridCurrentLocal.imag() +
                          filterVoltageLocal.imag() * gridCurrentLocal.real();

    Matrix controllerState = Matrix::Zero(ControllerStateCount, 1);
    controllerState(ThetaPLL, 0) = theta;
    controllerState(PhiPLL, 0) = 0.0;
    controllerState(PFiltered, 0) = pInitial;
    controllerState(QFiltered, 0) = qInitial;
    controllerState(PhiD, 0) =
        (gridCurrentLocal.real() + mKpPowerCtrl * (pInitial - mPRef)) /
        mKiPowerCtrl;
    controllerState(PhiQ, 0) =
        (gridCurrentLocal.imag() - mKpPowerCtrl * (qInitial - mQRef)) /
        mKiPowerCtrl;

    const Real currentReferenceD = mKpPowerCtrl * (mPRef - pInitial) +
                                   mKiPowerCtrl * controllerState(PhiD, 0);
    const Real currentReferenceQ = mKpPowerCtrl * (qInitial - mQRef) +
                                   mKiPowerCtrl * controllerState(PhiQ, 0);
    controllerState(GammaD, 0) =
        (converterVoltageLocal.real() +
         mKpCurrCtrl * (gridCurrentLocal.real() - currentReferenceD)) /
        mKiCurrCtrl;
    controllerState(GammaQ, 0) =
        (converterVoltageLocal.imag() +
         mKpCurrCtrl * (gridCurrentLocal.imag() - currentReferenceQ)) /
        mKiCurrCtrl;

    Matrix electricalState = Matrix::Zero(ElectricalStateCount, 1);
    electricalState.block(VcOffset, 0, 3, 1) = dq0Vector(filterVoltage);
    electricalState.block(IfOffset, 0, 3, 1) = dq0Vector(filterCurrent);
    electricalState.block(IGridOffset, 0, 3, 1) = dq0Vector(gridCurrent);

    Matrix measurement = Matrix::Zero(6, 1);
    measurement.block(0, 0, 3, 1) = electricalState.block(VcOffset, 0, 3, 1);
    measurement.block(3, 0, 3, 1) = electricalState.block(IGridOffset, 0, 3, 1);

    const SteadyState result{
        sourceVoltage,
        seriesVoltage,
        pccVoltage,
        filterVoltage,
        filterCurrent,
        gridCurrent,
        converterVoltage,
        controllerState,
        electricalState,
        measurement,
        dq0Vector(converterVoltage),
    };

    Matrix completeState = Matrix::Zero(NoDelayStateCount, 1);
    completeState.block(0, 0, ControllerStateCount, 1) = controllerState;
    completeState.block(ControllerStateCount, 0, ElectricalStateCount, 1) =
        electricalState;
    const Real residualNorm = noDelayDerivative(completeState).norm();
    if (residualNorm > 1e-6)
      throw std::runtime_error(
          "Manual operating point is not an equilibrium of the independent "
          "dq0 equations. Residual norm: " +
          std::to_string(residualNorm));
    return result;
  }

  Matrix controllerDerivative(const Matrix &controllerState,
                              const Matrix &measurement) const {
    Matrix derivative = Matrix::Zero(ControllerStateCount, 1);
    const Matrix localTransform =
        localFromGlobalDq0(controllerState(ThetaPLL, 0));
    const Matrix filterVoltageLocal =
        localTransform * measurement.block(0, 0, 3, 1);
    const Matrix gridCurrentLocal =
        localTransform * measurement.block(3, 0, 3, 1);

    const Real pInstantaneous =
        filterVoltageLocal(0, 0) * gridCurrentLocal(0, 0) +
        filterVoltageLocal(1, 0) * gridCurrentLocal(1, 0);
    const Real qInstantaneous =
        -filterVoltageLocal(0, 0) * gridCurrentLocal(1, 0) +
        filterVoltageLocal(1, 0) * gridCurrentLocal(0, 0);

    derivative(ThetaPLL, 0) =
        mKpPLL * filterVoltageLocal(1, 0) + mKiPLL * controllerState(PhiPLL, 0);
    derivative(PhiPLL, 0) = filterVoltageLocal(1, 0);
    derivative(PFiltered, 0) =
        mOmegaCutoff * (pInstantaneous - controllerState(PFiltered, 0));
    derivative(QFiltered, 0) =
        mOmegaCutoff * (qInstantaneous - controllerState(QFiltered, 0));
    derivative(PhiD, 0) = mPRef - controllerState(PFiltered, 0);
    derivative(PhiQ, 0) = controllerState(QFiltered, 0) - mQRef;

    const Real currentReferenceD =
        mKpPowerCtrl * (mPRef - controllerState(PFiltered, 0)) +
        mKiPowerCtrl * controllerState(PhiD, 0);
    const Real currentReferenceQ =
        mKpPowerCtrl * (controllerState(QFiltered, 0) - mQRef) +
        mKiPowerCtrl * controllerState(PhiQ, 0);
    derivative(GammaD, 0) = currentReferenceD - gridCurrentLocal(0, 0);
    derivative(GammaQ, 0) = currentReferenceQ - gridCurrentLocal(1, 0);
    return derivative;
  }

  Matrix controllerOutput(const Matrix &controllerState,
                          const Matrix &measurement) const {
    const Matrix localTransform =
        localFromGlobalDq0(controllerState(ThetaPLL, 0));
    const Matrix gridCurrentLocal =
        localTransform * measurement.block(3, 0, 3, 1);
    const Real currentReferenceD =
        mKpPowerCtrl * (mPRef - controllerState(PFiltered, 0)) +
        mKiPowerCtrl * controllerState(PhiD, 0);
    const Real currentReferenceQ =
        mKpPowerCtrl * (controllerState(QFiltered, 0) - mQRef) +
        mKiPowerCtrl * controllerState(PhiQ, 0);

    Matrix localOutput = Matrix::Zero(3, 1);
    localOutput(0, 0) =
        mKpCurrCtrl * (currentReferenceD - gridCurrentLocal(0, 0)) +
        mKiCurrCtrl * controllerState(GammaD, 0);
    localOutput(1, 0) =
        mKpCurrCtrl * (currentReferenceQ - gridCurrentLocal(1, 0)) +
        mKiCurrCtrl * controllerState(GammaQ, 0);
    return localTransform.transpose() * localOutput;
  }

  Matrix electricalDerivative(const Matrix &electricalState,
                              const Matrix &converterVoltage) const {
    const Matrix filterVoltage = electricalState.block(VcOffset, 0, 3, 1);
    const Matrix filterCurrent = electricalState.block(IfOffset, 0, 3, 1);
    const Matrix gridCurrent = electricalState.block(IGridOffset, 0, 3, 1);
    const Matrix sourceVoltage =
        dq0Vector(Complex(mGridVoltageRmsLineToLine, 0.0));

    Matrix derivative = Matrix::Zero(ElectricalStateCount, 1);
    derivative.block(VcOffset, 0, 3, 1) =
        (filterCurrent - gridCurrent) / mCf +
        rotatingFrameDerivative(filterVoltage, mOmega);
    derivative.block(IfOffset, 0, 3, 1) =
        (converterVoltage - filterVoltage - mRf * filterCurrent) / mLf +
        rotatingFrameDerivative(filterCurrent, mOmega);
    derivative.block(IGridOffset, 0, 3, 1) =
        (filterVoltage - sourceVoltage -
         (mGridResistance + mRc) * gridCurrent) /
            mGridInductance +
        rotatingFrameDerivative(gridCurrent, mOmega);
    return derivative;
  }

  Matrix noDelayDerivative(const Matrix &state) const {
    const Matrix controllerState = state.block(0, 0, ControllerStateCount, 1);
    const Matrix electricalState =
        state.block(ControllerStateCount, 0, ElectricalStateCount, 1);
    Matrix measurement = Matrix::Zero(6, 1);
    measurement.block(0, 0, 3, 1) = electricalState.block(VcOffset, 0, 3, 1);
    measurement.block(3, 0, 3, 1) = electricalState.block(IGridOffset, 0, 3, 1);

    Matrix derivative = Matrix::Zero(NoDelayStateCount, 1);
    derivative.block(0, 0, ControllerStateCount, 1) =
        controllerDerivative(controllerState, measurement);
    derivative.block(ControllerStateCount, 0, ElectricalStateCount, 1) =
        electricalDerivative(electricalState,
                             controllerOutput(controllerState, measurement));
    return derivative;
  }

  Matrix buildIndependentNoDelayContinuousMatrix(
      const SteadyState &steadyState) const {
    Matrix operatingPoint = Matrix::Zero(NoDelayStateCount, 1);
    operatingPoint.block(0, 0, ControllerStateCount, 1) =
        steadyState.controllerState;
    operatingPoint.block(ControllerStateCount, 0, ElectricalStateCount, 1) =
        steadyState.electricalState;
    return numericalJacobian(
        [this](const Matrix &state) { return noDelayDerivative(state); },
        operatingPoint);
  }

  Matrix buildIndependentDelayedDiscreteMatrix(const SteadyState &steadyState,
                                               Real timeStep) const {
    const Matrix controllerA = numericalJacobian(
        [this, &steadyState](const Matrix &controllerState) {
          return controllerDerivative(controllerState,
                                      steadyState.controllerMeasurement);
        },
        steadyState.controllerState);
    const Matrix controllerB = numericalJacobian(
        [this, &steadyState](const Matrix &measurement) {
          return controllerDerivative(steadyState.controllerState, measurement);
        },
        steadyState.controllerMeasurement);
    const Matrix controllerC = numericalJacobian(
        [this, &steadyState](const Matrix &controllerState) {
          return controllerOutput(controllerState,
                                  steadyState.controllerMeasurement);
        },
        steadyState.controllerState);
    const Matrix controllerD = numericalJacobian(
        [this, &steadyState](const Matrix &measurement) {
          return controllerOutput(steadyState.controllerState, measurement);
        },
        steadyState.controllerMeasurement);
    const Matrix electricalA = numericalJacobian(
        [this, &steadyState](const Matrix &electricalState) {
          return electricalDerivative(electricalState,
                                      steadyState.delayedControllerOutput);
        },
        steadyState.electricalState);
    const Matrix electricalB = numericalJacobian(
        [this, &steadyState](const Matrix &converterVoltage) {
          return electricalDerivative(steadyState.electricalState,
                                      converterVoltage);
        },
        steadyState.delayedControllerOutput);

    const Matrix controllerAd = trapezoidalStateMatrix(controllerA, timeStep);
    const Matrix controllerBd =
        trapezoidalInputMatrix(controllerA, controllerB, timeStep);
    const Matrix electricalAd = trapezoidalStateMatrix(electricalA, timeStep);
    const Matrix electricalBd =
        trapezoidalInputMatrix(electricalA, electricalB, timeStep);

    Matrix measurementFromElectrical = Matrix::Zero(6, ElectricalStateCount);
    measurementFromElectrical.block(0, VcOffset, 3, 3) = Matrix::Identity(3, 3);
    measurementFromElectrical.block(3, IGridOffset, 3, 3) =
        Matrix::Identity(3, 3);

    const Matrix controllerFromElectrical =
        controllerBd * measurementFromElectrical *
        (Matrix::Identity(ElectricalStateCount, ElectricalStateCount) +
         electricalAd);
    const Matrix controllerFromDelay =
        2.0 * controllerBd * measurementFromElectrical * electricalBd;
    const Matrix measurementNewFromElectrical =
        measurementFromElectrical * electricalAd;
    const Matrix measurementNewFromDelay =
        2.0 * measurementFromElectrical * electricalBd;

    Matrix result = Matrix::Zero(DelayedStateCount, DelayedStateCount);
    const UInt controllerOffset = 0;
    const UInt electricalOffset = ControllerStateCount;
    const UInt delayOffset = NoDelayStateCount;

    result.block(controllerOffset, controllerOffset, ControllerStateCount,
                 ControllerStateCount) = controllerAd;
    result.block(controllerOffset, electricalOffset, ControllerStateCount,
                 ElectricalStateCount) = controllerFromElectrical;
    result.block(controllerOffset, delayOffset, ControllerStateCount,
                 DelayStateCount) = controllerFromDelay;
    result.block(electricalOffset, electricalOffset, ElectricalStateCount,
                 ElectricalStateCount) = electricalAd;
    result.block(electricalOffset, delayOffset, ElectricalStateCount,
                 DelayStateCount) = 2.0 * electricalBd;
    result.block(delayOffset, controllerOffset, DelayStateCount,
                 ControllerStateCount) = controllerC * controllerAd;
    result.block(delayOffset, electricalOffset, DelayStateCount,
                 ElectricalStateCount) =
        controllerC * controllerFromElectrical +
        controllerD * measurementNewFromElectrical;
    result.block(delayOffset, delayOffset, DelayStateCount, DelayStateCount) =
        controllerC * controllerFromDelay +
        controllerD * measurementNewFromDelay;
    return result;
  }

  template <typename InverterType>
  ExtractionResult runExtraction(const SteadyState &steadyState,
                                 const std::shared_ptr<InverterType> &inverter,
                                 const String &simName, Real timeStep,
                                 UInt expectedStateCount) const {
    Logger::setLogDir("logs/" + simName);
    const UInt logDownsampling = std::max<UInt>(
        1, static_cast<UInt>(std::llround(mLogInterval / timeStep)));
    auto logger = DataLogger::make(simName, true, logDownsampling);
    auto system = createEmtSystem(steadyState, inverter, logger);
    const UInt extractionStep =
        stepCountForTime(mExtractionTime, timeStep, "extraction time");
    const UInt finalStep =
        stepCountForTime(mFinalTime, timeStep, "final simulation time");

    Simulation simulation(simName, Logger::Level::debug);
    simulation.setSystem(system);
    simulation.addLogger(logger);
    simulation.setDomain(Domain::EMT);
    simulation.setSolverType(Solver::Type::MNA);
    simulation.setTimeStep(timeStep);
    simulation.setFinalTime(mFinalTime);
    simulation.doStateSpaceExtraction(true);
    simulation.doInitFromNodesAndTerminals(true);
    simulation.start();

    ExtractionResult result;
    Bool extractionCaptured = false;
    for (UInt step = 1; step <= finalStep; ++step) {
      simulation.next();
      if (step != extractionStep)
        continue;
      const auto &extractor = simulation.getStateSpaceExtractor();
      if (extractor.getStateCount() != expectedStateCount)
        throw std::runtime_error("Unexpected extracted state count in " +
                                 simName + ".");
      StateSpaceModalAnalysis modalAnalysis(extractor);
      modalAnalysis.setAnalysisFrame(StateSpaceAnalysisFrame::GlobalDQ0);
      modalAnalysis.setGlobalDq0Frame(mOmega);
      modalAnalysis.update();
      result = {modalAnalysis.getDiscreteEigenvalues(),
                extractor.getLastExtractionTime(), extractor.getStateCount()};
      extractionCaptured = true;
    }
    simulation.stop();
    if (!extractionCaptured)
      throw std::runtime_error("State-space extraction was not captured in " +
                               simName + ".");
    return result;
  }

  UInt stepCountForTime(Real requestedTime, Real timeStep,
                        const String &description) const {
    const auto roundedSteps = std::llround(requestedTime / timeStep);
    if (roundedSteps < 1)
      throw std::invalid_argument(description +
                                  " must be at least one simulation step.");
    const Real alignedTime = static_cast<Real>(roundedSteps) * timeStep;
    const Real tolerance =
        1e-9 * std::max({1.0, std::abs(requestedTime), std::abs(timeStep)});
    if (std::abs(alignedTime - requestedTime) > tolerance)
      throw std::invalid_argument(description + " (" +
                                  std::to_string(requestedTime) +
                                  " s) must be an integer multiple of the "
                                  "time step (" +
                                  std::to_string(timeStep) + " s).");
    return static_cast<UInt>(roundedSteps);
  }

  template <typename InverterType>
  SystemTopology createEmtSystem(const SteadyState &steadyState,
                                 const std::shared_ptr<InverterType> &inverter,
                                 const DataLogger::Ptr &logger) const {
    auto nGrid = SimNode<Real>::make("nGrid", PhaseType::ABC);
    auto nSeries = SimNode<Real>::make("nSeries", PhaseType::ABC);
    auto nPcc = SimNode<Real>::make("nPcc", PhaseType::ABC);
    nGrid->setInitialVoltage(steadyState.sourceVoltage);
    nSeries->setInitialVoltage(steadyState.seriesVoltage);
    nPcc->setInitialVoltage(steadyState.pccVoltage);

    auto slack = EMT::Ph3::NetworkInjection::make("Slack");
    slack->setParameters(
        Math::singlePhaseVariableToThreePhase(steadyState.sourceVoltage),
        mFrequency);
    auto resistance = EMT::Ph3::Resistor::make("GridResistance");
    resistance->setParameters(
        Math::singlePhaseParameterToThreePhase(mGridResistance));
    auto inductance = EMT::Ph3::Inductor::make("GridInductance");
    inductance->setParameters(
        Math::singlePhaseParameterToThreePhase(mGridInductance));

    slack->connect({nGrid});
    inductance->connect({nGrid, nSeries});
    resistance->connect({nSeries, nPcc});
    inverter->connect({EMT::SimNode::GND, nPcc});

    logger->logAttribute("v_grid", nGrid->attribute("v"));
    logger->logAttribute("v_series", nSeries->attribute("v"));
    logger->logAttribute("v_pcc", nPcc->attribute("v"));
    logger->logAttribute("i_inv", inverter->attribute("i_intf"));
    logger->logAttribute("vc_d", inverter->attribute("vc_d"));
    logger->logAttribute("vc_q", inverter->attribute("vc_q"));
    logger->logAttribute("igrid_d", inverter->attribute("irc_d"));
    logger->logAttribute("igrid_q", inverter->attribute("irc_q"));
    logger->logAttribute("p_inst", inverter->attribute("p_inst"));
    logger->logAttribute("q_inst", inverter->attribute("q_inst"));
    logger->logAttribute("omega_pll", inverter->attribute("omega_pll"));
    addStateSignal(logger, inverter);

    return SystemTopology(
        mFrequency, SystemNodeList{nGrid, nSeries, nPcc},
        SystemComponentList{slack, inductance, resistance, inverter});
  }

  void
  addStateSignal(const DataLogger::Ptr &logger,
                 const std::shared_ptr<EMT::Ph3::SSN_GFL> &inverter) const {
    logger->logAttribute("state", inverter->attribute("x"));
  }

  void addStateSignal(
      const DataLogger::Ptr &logger,
      const std::shared_ptr<EMT::Ph3::SSN_GFL_Split> &inverter) const {
    logger->logAttribute("state", inverter->getSplitStateAttribute());
  }

  Real mFrequency;
  Real mOmega;
  Real mGridVoltageRmsLineToLine;
  Real mGridResistance;
  Real mGridInductance;
  Real mLf;
  Real mCf;
  Real mRf;
  Real mRc;
  Real mKpPLL;
  Real mKiPLL;
  Real mOmegaCutoff;
  Real mPRef;
  Real mQRef;
  Real mKpPowerCtrl;
  Real mKiPowerCtrl;
  Real mKpCurrCtrl;
  Real mKiCurrCtrl;
  Real mExtractionTime;
  Real mFinalTime;
  Real mLogInterval = 100e-6;
};

int main(int argc, char *argv[]) {
  try {
    if (argc > 3) {
      std::cerr << "Usage: " << argv[0]
                << " [extraction_time_s] [final_time_s]\n";
      return 1;
    }
    const Real extractionTime = argc >= 2 ? std::stod(argv[1]) : 0.1;
    const Real finalTime = argc >= 3 ? std::stod(argv[2]) : 0.1;
    EMTPh3SplitGFLStateSpaceExtractionExample example(extractionTime,
                                                      finalTime);
    example.run();
    return 0;
  } catch (const std::exception &exception) {
    std::cerr << "Error: " << exception.what() << '\n';
    return 1;
  }
}
