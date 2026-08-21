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
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/SVD>

using namespace CPS;
using namespace DPsim;

namespace {

constexpr UInt StateCount = 14;

enum StateIndex : UInt {
  Psi = 0,
  PhiPLL,
  PFiltered,
  QFiltered,
  PhiD,
  PhiQ,
  GammaD,
  GammaQ,
  VcD,
  VcQ,
  IfD,
  IfQ,
  IGridD,
  IGridQ,
};

struct SteadyState {
  Complex sourceVoltage;
  Complex midVoltage;
  Complex pccVoltage;
  Matrix state;
};

struct SimulationHandles {
  SystemTopology system;
  std::shared_ptr<DP::Ph1::AvVoltSourceInverterStateSpace> inverter;
  std::shared_ptr<DP::Ph1::Inductor> gridInductance;
};

struct ExtractionResult {
  VectorComp discreteEigenvalues;
  Matrix discreteStateMatrix;
  std::vector<String> stateNames;
  Matrix simulationState;
  Matrix previousInverterState;
  Complex terminalVoltage;
  Complex inverterInterfaceCurrent;
  Complex gridCurrent;
  Complex nonlinearGridCurrent;
  Int iterationCount;
  Real iterationStateResidual;
  Real iterationInputResidual;
  Real extractionTime;
  UInt stateCount;
  Matrix fullNortonAdmittance;
  Matrix scalarNortonAdmittance;
  Matrix nortonStructureDefect;
  UInt discreteModelRevision;
  UInt mnaStampRevision;
  UInt extractionStampRevision;
};

struct EigenvalueRecord {
  Real timeStep;
  String model;
  VectorComp values;
};

struct MatrixRecord {
  Real timeStep;
  String solutionMode;
  String matrixName;
  String coordinateSystem;
  Matrix values;
  std::vector<String> stateNames;
};

struct DiagnosticRecord {
  Real timeStep;
  String solutionMode;
  Real extractionTime;
  Int iterationCount;
  Real iterationStateResidual;
  Real iterationInputResidual;
  Real stateStepDistance;
  Real operatingPointDistance;
  Real continuousResidual;
  Real controllerResidual;
  Real filterResidual;
  Real gridResidual;
  Complex nonlinearGridCurrentMismatch;
  Complex mnaKclMismatch;
  Real nortonDefectNorm;
  Real nortonDefectMax;
  UInt discreteModelRevision;
  UInt mnaStampRevision;
  UInt extractionStampRevision;
  Real extractedCondition;
  Real extractedCayleyCondition;
  Real referenceCondition;
  Real referenceCayleyCondition;
  Real extractedMaxEntry;
  Real referenceMaxEntry;
  Real extractedEigenBackwardError;
  Real referenceEigenBackwardError;
};

VectorComp eigenvalues(const Matrix &matrix) {
  Eigen::EigenSolver<Matrix> solver(matrix);
  if (solver.info() != Eigen::Success)
    throw std::runtime_error("Eigenvalue computation failed.");
  return solver.eigenvalues();
}

Real directedEigenvalueDistance(const VectorComp &reference,
                                const VectorComp &actual) {
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

Real eigenvalueDistance(const VectorComp &reference,
                        const VectorComp &actual) {
  return std::max(directedEigenvalueDistance(reference, actual),
                  directedEigenvalueDistance(actual, reference));
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

Matrix trapezoidalStateMatrix(const Matrix &continuousA, Real timeStep) {
  const Matrix identity =
      Matrix::Identity(continuousA.rows(), continuousA.cols());
  return (identity - 0.5 * timeStep * continuousA)
      .fullPivLu()
      .solve(identity + 0.5 * timeStep * continuousA);
}

Real matrixConditionNumber(const Matrix &matrix) {
  if (matrix.rows() == 0 || matrix.cols() == 0)
    return 0.0;

  Eigen::JacobiSVD<Matrix> decomposition(matrix);
  const auto &singularValues = decomposition.singularValues();
  const Real smallest = singularValues(singularValues.rows() - 1);
  if (smallest == 0.0)
    return std::numeric_limits<Real>::infinity();
  return singularValues(0) / smallest;
}

Real matrixMaxAbsEntry(const Matrix &matrix) {
  return matrix.cwiseAbs().maxCoeff();
}

Real eigenvalueBackwardError(const Matrix &matrix) {
  Eigen::EigenSolver<Matrix> solver(matrix, true);
  if (solver.info() != Eigen::Success)
    throw std::runtime_error("Eigenvector computation failed.");

  const VectorComp values = solver.eigenvalues();
  const MatrixComp vectors = solver.eigenvectors();
  const Real matrixNorm = matrix.norm();
  Real maximum = 0.0;
  for (Eigen::Index index = 0; index < values.rows(); ++index) {
    const VectorComp vector = vectors.col(index);
    const Real denominator =
        (matrixNorm + std::abs(values(index))) * vector.norm();
    if (denominator > 0.0) {
      const Real residual =
          (matrix.cast<Complex>() * vector - values(index) * vector).norm();
      maximum = std::max(maximum, residual / denominator);
    }
  }
  return maximum;
}

Real scaledOperatingPointDistance(const Matrix &reference,
                                  const Matrix &actual) {
  Matrix difference(reference.rows(), 1);
  for (Eigen::Index index = 0; index < reference.rows(); ++index)
    difference(index, 0) =
        (actual(index, 0) - reference(index, 0)) /
        std::max(1.0, std::abs(reference(index, 0)));
  return difference.norm() / std::sqrt(static_cast<Real>(reference.rows()));
}

void writeEigenvaluesCsv(const std::filesystem::path &path,
                         const std::vector<EigenvalueRecord> &records) {
  std::filesystem::create_directories(path.parent_path());
  std::ofstream stream(path);
  if (!stream.is_open())
    throw std::runtime_error("Could not open eigenvalue CSV output: " +
                             path.string());
  stream << std::setprecision(std::numeric_limits<Real>::max_digits10)
         << "time_step_s,time_step_us,model,index,real,imag\n";
  for (const auto &record : records) {
    for (Eigen::Index index = 0; index < record.values.rows(); ++index)
      stream << record.timeStep << ',' << 1e6 * record.timeStep << ','
             << record.model << ',' << index << ','
             << record.values(index).real() << ','
             << record.values(index).imag() << '\n';
  }
  if (!stream)
    throw std::runtime_error("Could not write eigenvalue CSV output: " +
                             path.string());
}

void writeMatricesCsv(const std::filesystem::path &path,
                      const std::vector<MatrixRecord> &records) {
  std::filesystem::create_directories(path.parent_path());
  std::ofstream stream(path);
  if (!stream.is_open())
    throw std::runtime_error("Could not open matrix CSV output: " +
                             path.string());
  stream << std::setprecision(std::numeric_limits<Real>::max_digits10)
         << "time_step_s,time_step_us,solution_mode,matrix_name,"
            "coordinate_system,row,column,row_state,column_state,value\n";
  for (const auto &record : records) {
    for (Eigen::Index row = 0; row < record.values.rows(); ++row) {
      for (Eigen::Index column = 0; column < record.values.cols(); ++column) {
        stream << record.timeStep << ',' << 1e6 * record.timeStep << ','
               << record.solutionMode << ',' << record.matrixName << ','
               << record.coordinateSystem << ',' << row << ',' << column
               << ',' << record.stateNames[static_cast<std::size_t>(row)]
               << ','
               << record.stateNames[static_cast<std::size_t>(column)] << ','
               << record.values(row, column) << '\n';
      }
    }
  }
  if (!stream)
    throw std::runtime_error("Could not write matrix CSV output: " +
                             path.string());
}

void writeDiagnosticsCsv(const std::filesystem::path &path,
                         const std::vector<DiagnosticRecord> &records) {
  std::filesystem::create_directories(path.parent_path());
  std::ofstream stream(path);
  if (!stream.is_open())
    throw std::runtime_error("Could not open diagnostics CSV output: " +
                             path.string());
  stream << std::setprecision(std::numeric_limits<Real>::max_digits10)
         << "time_step_s,time_step_us,solution_mode,extraction_time_s,"
            "iteration_count,iteration_state_residual,"
            "iteration_input_residual,state_step_distance,"
            "operating_point_distance,continuous_residual,"
            "controller_residual,filter_residual,grid_residual,"
            "nonlinear_grid_current_mismatch_re,"
            "nonlinear_grid_current_mismatch_im,"
            "nonlinear_grid_current_mismatch_abs,mna_kcl_mismatch_re,"
            "mna_kcl_mismatch_im,mna_kcl_mismatch_abs,"
            "norton_structure_defect_norm,norton_structure_defect_max,"
            "discrete_model_revision,mna_stamp_revision,"
            "extraction_stamp_revision,extracted_condition,"
            "extracted_cayley_condition,reference_condition,"
            "reference_cayley_condition,extracted_max_entry,"
            "reference_max_entry,extracted_eigen_backward_error,"
            "reference_eigen_backward_error\n";
  for (const auto &record : records) {
    stream << record.timeStep << ',' << 1e6 * record.timeStep << ','
           << record.solutionMode << ',' << record.extractionTime << ','
           << record.iterationCount << ',' << record.iterationStateResidual
           << ',' << record.iterationInputResidual << ','
           << record.stateStepDistance << ',' << record.operatingPointDistance
           << ',' << record.continuousResidual << ','
           << record.controllerResidual << ',' << record.filterResidual << ','
           << record.gridResidual << ','
           << record.nonlinearGridCurrentMismatch.real() << ','
           << record.nonlinearGridCurrentMismatch.imag() << ','
           << std::abs(record.nonlinearGridCurrentMismatch) << ','
           << record.mnaKclMismatch.real() << ','
           << record.mnaKclMismatch.imag() << ','
           << std::abs(record.mnaKclMismatch) << ','
           << record.nortonDefectNorm << ',' << record.nortonDefectMax << ','
           << record.discreteModelRevision << ',' << record.mnaStampRevision
           << ',' << record.extractionStampRevision << ','
           << record.extractedCondition << ','
           << record.extractedCayleyCondition << ','
           << record.referenceCondition << ','
           << record.referenceCayleyCondition << ','
           << record.extractedMaxEntry << ',' << record.referenceMaxEntry
           << ',' << record.extractedEigenBackwardError << ','
           << record.referenceEigenBackwardError << '\n';
  }
  if (!stream)
    throw std::runtime_error("Could not write diagnostics CSV output: " +
                             path.string());
}

std::vector<String> physicalStateNames() {
  return {"psi",       "phi_pll", "p_filtered", "q_filtered", "phi_d",
          "phi_q",     "gamma_d", "gamma_q",    "vc_d",       "vc_q",
          "if_d",      "if_q",    "i_grid_d",   "i_grid_q"};
}

} // namespace

class DPPh1InverterStateSpaceExtractionExample {
public:
  DPPh1InverterStateSpaceExtractionExample(Real extractionTime = 0.1,
                                           Real finalTime = 0.1)
      : mFrequency(50.0), mOmega(2.0 * PI * mFrequency),
        mSourceVoltage(RMS3PH_TO_PEAK1PH * 400.0, 0.0), mGridResistance(0.3),
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
    const String baseName = "DP_Ph1_Inverter_StateSpaceExtraction";
    const SteadyState manualPoint = calculateSteadyState();
    const Matrix manualA = buildIndependentDqJacobian(manualPoint.state);
    const Real manualResidual = independentDqDerivative(manualPoint.state).norm();
    const std::array<Real, 4> timeSteps{1e-6, 10e-6, 100e-6, 1e-3};
    const std::array<String, 4> labels{"1us", "10us", "100us", "1ms"};
    const std::array<Bool, 2> iterationModes{false, true};
    const std::array<Bool, 2> fullNortonMatrixModes{false, true};
    std::vector<EigenvalueRecord> eigenvalueRecords;
    std::vector<MatrixRecord> matrixRecords;
    std::vector<DiagnosticRecord> diagnosticRecords;

    std::cout
        << "\n============================================================\n"
        << "DP Ph1 inverter state-space extraction comparison\n"
        << "============================================================\n"
        << "Topology: ideal DP source -> R/L grid -> averaged inverter\n"
        << "Analysis frame: native DP/global synchronous dq\n"
        << "Reference: independent nonlinear dq equations, numerical "
           "linearization, trapezoidal discretization\n"
        << "Norton stamp modes: complex scalar and full real 2x2 matrix\n"
        << "Solution modes: non-iterative and iterative\n"
        << "Time steps: 1 us, 10 us, 100 us, 1 ms\n"
        << "Extraction time: " << mExtractionTime << " s\n"
        << "Final simulation time: " << mFinalTime << " s\n"
        << "Manual operating-point residual norm: " << manualResidual << "\n";

    for (UInt index = 0; index < timeSteps.size(); ++index) {
      const Real timeStep = timeSteps[index];
      const Matrix manualAd = trapezoidalStateMatrix(manualA, timeStep);
      const VectorComp manualZ =
          eigenvalues(manualAd);

      std::cout << "\nDelta t = " << 1e6 * timeStep << " us\n";
      for (const Bool fullNortonMatrix : fullNortonMatrixModes) {
        const String stampMode =
            fullNortonMatrix ? "full_matrix" : "complex_scalar";
        for (const Bool iterativeSolution : iterationModes) {
          const String iterationMode =
              iterativeSolution ? "iterative" : "non_iterative";
          const String mode = stampMode + "_" + iterationMode;
          const String simName =
              baseName + "_" + mode + "_dt_" + labels[index];
          const ExtractionResult extracted =
              runExtraction(manualPoint, simName, timeStep,
                            iterativeSolution, fullNortonMatrix);
        const Matrix simulationA =
            buildIndependentDqJacobian(extracted.simulationState);
        const Matrix simulationAd =
            trapezoidalStateMatrix(simulationA, timeStep);
        const VectorComp simulationZ = eigenvalues(simulationAd);
        const Matrix residual =
            independentDqDerivative(extracted.simulationState);
        const Complex nonlinearCurrentMismatch =
            extracted.nonlinearGridCurrent - extracted.gridCurrent;
        const Complex mnaKclMismatch =
            extracted.inverterInterfaceCurrent + extracted.gridCurrent;
        const Real stateStepDistance = scaledOperatingPointDistance(
            extracted.previousInverterState,
            extracted.simulationState.block(0, 0, 12, 1));
        const Real operatingPointDistance = scaledOperatingPointDistance(
            manualPoint.state, extracted.simulationState);
        const Matrix extractedCayley =
            extracted.discreteStateMatrix +
            Matrix::Identity(extracted.discreteStateMatrix.rows(),
                             extracted.discreteStateMatrix.cols());
        const Matrix referenceCayley =
            manualAd + Matrix::Identity(manualAd.rows(), manualAd.cols());
        const Real extractedCondition =
            matrixConditionNumber(extracted.discreteStateMatrix);
        const Real extractedCayleyCondition =
            matrixConditionNumber(extractedCayley);
        const Real referenceCondition = matrixConditionNumber(manualAd);
        const Real referenceCayleyCondition =
            matrixConditionNumber(referenceCayley);

        eigenvalueRecords.push_back(
            {timeStep, "extracted_dpsim_" + mode,
             extracted.discreteEigenvalues});
        eigenvalueRecords.push_back(
            {timeStep, "independent_dq_manual_pf_" + mode, manualZ});
        eigenvalueRecords.push_back(
            {timeStep, "independent_dq_simulation_operating_point_" + mode,
             simulationZ});

        matrixRecords.push_back(
            {timeStep, mode, "extracted_ad", "companion_history",
             extracted.discreteStateMatrix, extracted.stateNames});
        matrixRecords.push_back(
            {timeStep, mode, "analytical_manual_pf_ad", "physical_dq",
             manualAd, physicalStateNames()});
        matrixRecords.push_back(
            {timeStep, mode, "analytical_simulation_op_ad", "physical_dq",
             simulationAd, physicalStateNames()});
        matrixRecords.push_back(
            {timeStep, mode, "norton_full_w", "packed_real_dq",
             extracted.fullNortonAdmittance, {"u_re", "u_im"}});
        matrixRecords.push_back(
            {timeStep, mode, "norton_complex_scalar_w", "packed_real_dq",
             extracted.scalarNortonAdmittance, {"u_re", "u_im"}});
        matrixRecords.push_back(
            {timeStep, mode, "norton_structure_defect", "packed_real_dq",
             extracted.nortonStructureDefect, {"u_re", "u_im"}});

        diagnosticRecords.push_back(
            {timeStep,
             mode,
             extracted.extractionTime,
             extracted.iterationCount,
             extracted.iterationStateResidual,
             extracted.iterationInputResidual,
             stateStepDistance,
             operatingPointDistance,
             residual.norm(),
             residual.block(0, 0, 8, 1).norm(),
             residual.block(8, 0, 4, 1).norm(),
             residual.block(12, 0, 2, 1).norm(),
             nonlinearCurrentMismatch,
             mnaKclMismatch,
             extracted.nortonStructureDefect.norm(),
             extracted.nortonStructureDefect.cwiseAbs().maxCoeff(),
             extracted.discreteModelRevision,
             extracted.mnaStampRevision,
             extracted.extractionStampRevision,
             extractedCondition,
             extractedCayleyCondition,
             referenceCondition,
             referenceCayleyCondition,
             matrixMaxAbsEntry(extracted.discreteStateMatrix),
             matrixMaxAbsEntry(manualAd),
             eigenvalueBackwardError(extracted.discreteStateMatrix),
             eigenvalueBackwardError(manualAd)});

        std::cout
            << "  " << mode << ":\n"
            << "    extraction time used: " << extracted.extractionTime
            << " s\n"
            << "    extraction states: " << extracted.stateCount << "\n"
            << "    additional MNA iterations: " << extracted.iterationCount
            << "\n"
            << "    final iteration state residual: "
            << extracted.iterationStateResidual << "\n"
            << "    final iteration input residual: "
            << extracted.iterationInputResidual << "\n"
            << "    Norton structure defect norm / max: "
            << extracted.nortonStructureDefect.norm() << " / "
            << extracted.nortonStructureDefect.cwiseAbs().maxCoeff() << "\n"
            << "    model/MNA/extraction revisions: "
            << extracted.discreteModelRevision << " / "
            << extracted.mnaStampRevision << " / "
            << extracted.extractionStampRevision
            << (extracted.discreteModelRevision == extracted.mnaStampRevision &&
                        extracted.discreteModelRevision ==
                            extracted.extractionStampRevision
                    ? " (match)\n"
                    : " (MISMATCH)\n")
            << "    extracted vs manual-PF dq reference max error: "
            << eigenvalueDistance(manualZ, extracted.discreteEigenvalues)
            << "\n"
            << "    extracted vs simulation-OP dq reference max error: "
            << eigenvalueDistance(simulationZ,
                                  extracted.discreteEigenvalues)
            << "\n"
            << "    manual-PF vs simulation-OP reference max error: "
            << eigenvalueDistance(manualZ, simulationZ) << "\n"
            << "    scaled state change in extraction step: "
            << stateStepDistance << "\n"
            << "    scaled simulation/manual operating-point distance: "
            << operatingPointDistance << "\n"
            << "    continuous residual (controller/filter/grid): "
            << residual.block(0, 0, 8, 1).norm() << " / "
            << residual.block(8, 0, 4, 1).norm() << " / "
            << residual.block(12, 0, 2, 1).norm() << "\n"
            << "    |nonlinear Irc - grid current|: "
            << std::abs(nonlinearCurrentMismatch) << " A\n"
            << "    |inverter MNA current + grid current|: "
            << std::abs(mnaKclMismatch) << " A\n"
            << "    cond(Ad), cond(Ad + I): " << extractedCondition << " / "
            << extractedCayleyCondition << "\n"
            << "    reference cond(Ad), cond(Ad + I): "
            << referenceCondition << " / " << referenceCayleyCondition
            << "\n"
            << "    eigenvalue backward error (extracted/reference): "
            << eigenvalueBackwardError(extracted.discreteStateMatrix) << " / "
            << eigenvalueBackwardError(manualAd) << "\n"
            << "    simulation log: logs/" << simName << '/' << simName
            << ".csv\n";
        }
      }
    }

    const std::filesystem::path outputDirectory =
        std::filesystem::path("logs") / baseName;
    const auto eigenvaluePath = outputDirectory / "eigenvalues.csv";
    const auto matrixPath = outputDirectory / "matrices.csv";
    const auto diagnosticPath = outputDirectory / "diagnostics.csv";
    writeEigenvaluesCsv(eigenvaluePath, eigenvalueRecords);
    writeMatricesCsv(matrixPath, matrixRecords);
    writeDiagnosticsCsv(diagnosticPath, diagnosticRecords);
    std::cout << "\nEigenvalue CSV: " << eigenvaluePath.string()
              << "\nMatrix CSV: " << matrixPath.string()
              << "\nDiagnostics CSV: " << diagnosticPath.string() << "\n";
  }

private:
  SimulationHandles createSystem(const SteadyState &steadyState,
                                 const DataLogger::Ptr &logger,
                                 Bool iterativeSolution,
                                 Bool fullNortonMatrix) const {
    auto nGrid = SimNode<Complex>::make("nGrid");
    auto nMid = SimNode<Complex>::make("nMid");
    auto nPcc = SimNode<Complex>::make("nPcc");
    nGrid->setInitialVoltage(steadyState.sourceVoltage);
    nMid->setInitialVoltage(steadyState.midVoltage);
    nPcc->setInitialVoltage(steadyState.pccVoltage);

    auto slack = DP::Ph1::VoltageSource::make("Slack");
    slack->setParameters(mSourceVoltage, 0.0);
    auto gridResistance = DP::Ph1::Resistor::make("GridResistance");
    gridResistance->setParameters(mGridResistance);
    auto gridInductance = DP::Ph1::Inductor::make("GridInductance");
    gridInductance->setParameters(mGridInductance);
    auto inverter = DP::Ph1::AvVoltSourceInverterStateSpace::make(
        "Inverter", Logger::Level::warn);
    inverter->setParameters(mLf, mCf, mRf, mRc, mOmega, mKpPLL, mKiPLL,
                            mOmegaCutoff, mPRef, mQRef, mKpPowerCtrl,
                            mKiPowerCtrl, mKpCurrCtrl, mKiCurrCtrl);
    inverter->setIterativeSolution(iterativeSolution);
    inverter->setNortonAdmittanceMode(
        fullNortonMatrix
            ? DP::Ph1::MixedVTypeVariableSSNComp::NortonAdmittanceMode::
                  FullRealMatrix
            : DP::Ph1::MixedVTypeVariableSSNComp::NortonAdmittanceMode::
                  ComplexScalar);

    slack->connect(SimNode<Complex>::List{SimNode<Complex>::GND, nGrid});
    gridResistance->connect(SimNode<Complex>::List{nGrid, nMid});
    gridInductance->connect(SimNode<Complex>::List{nMid, nPcc});
    inverter->connect(SimNode<Complex>::List{SimNode<Complex>::GND, nPcc});

    logger->logAttribute("v_grid", nGrid->attribute("v"));
    logger->logAttribute("v_mid", nMid->attribute("v"));
    logger->logAttribute("v_pcc", nPcc->attribute("v"));
    logger->logAttribute("i_grid_resistance",
                         gridResistance->attribute("i_intf"));
    logger->logAttribute("i_grid_inductance",
                         gridInductance->attribute("i_intf"));
    logger->logAttribute("i_inv", inverter->attribute("i_intf"));
    logger->logAttribute("inverter_state", inverter->attribute("x"));
    logger->logAttribute("vc_d", inverter->attribute("vc_d"));
    logger->logAttribute("vc_q", inverter->attribute("vc_q"));
    logger->logAttribute("irc_d", inverter->attribute("irc_d"));
    logger->logAttribute("irc_q", inverter->attribute("irc_q"));
    logger->logAttribute("p_inst", inverter->attribute("p_inst"));
    logger->logAttribute("q_inst", inverter->attribute("q_inst"));
    logger->logAttribute("omega_pll", inverter->attribute("omega_pll"));
    logger->logAttribute("iteration_count",
                         inverter->attribute("iteration_count"));
    logger->logAttribute("iteration_state_residual",
                         inverter->attribute("iteration_state_residual"));
    logger->logAttribute("iteration_input_residual",
                         inverter->attribute("iteration_input_residual"));

    SystemTopology system(
        mFrequency, SystemNodeList{nGrid, nMid, nPcc},
        SystemComponentList{slack, gridResistance, gridInductance, inverter});
    return {system, inverter, gridInductance};
  }

  ExtractionResult runExtraction(const SteadyState &steadyState,
                                 const String &simName,
                                 Real timeStep,
                                 Bool iterativeSolution,
                                 Bool fullNortonMatrix) const {
    Logger::setLogDir("logs/" + simName);
    const UInt logDownsampling = std::max<UInt>(
        1, static_cast<UInt>(std::llround(mLogInterval / timeStep)));
    auto logger = DataLogger::make(simName, true, logDownsampling);
    auto handles = createSystem(steadyState, logger, iterativeSolution,
                                fullNortonMatrix);
    const UInt extractionStep =
        stepCountForTime(mExtractionTime, timeStep, "extraction time");
    const UInt finalStep =
        stepCountForTime(mFinalTime, timeStep, "final simulation time");

    Simulation simulation(simName, Logger::Level::warn);
    simulation.setSystem(handles.system);
    simulation.addLogger(logger);
    simulation.setDomain(Domain::DP);
    simulation.setSolverType(Solver::Type::MNA);
    simulation.setTimeStep(timeStep);
    simulation.setFinalTime(mFinalTime);
    simulation.doStateSpaceExtraction(true);
    simulation.doSystemMatrixRecomputation(true);
    simulation.doInitFromNodesAndTerminals(true);
    simulation.start();

    ExtractionResult result;
    Bool captured = false;
    Matrix previousInverterState;
    for (UInt step = 1; step <= finalStep; ++step) {
      if (step == extractionStep)
        previousInverterState =
            handles.inverter->attributeTyped<Matrix>("x")->get();
      simulation.next();
      if (step != extractionStep)
        continue;

      const auto &extractor = simulation.getStateSpaceExtractor();
      if (extractor.getStateCount() != StateCount)
        throw std::runtime_error("Unexpected extracted state count in " +
                                 simName + ".");
      StateSpaceModalAnalysis modalAnalysis(extractor);
      modalAnalysis.update();

      const Matrix inverterState =
          handles.inverter->attributeTyped<Matrix>("x")->get();
      const MatrixComp inductorCurrent =
          handles.gridInductance->attributeTyped<MatrixComp>("i_intf")->get();
      const MatrixComp terminalVoltageMatrix =
          handles.inverter->attributeTyped<MatrixComp>("v_intf")->get();
      const MatrixComp inverterCurrentMatrix =
          handles.inverter->attributeTyped<MatrixComp>("i_intf")->get();
      const Complex terminalVoltage = terminalVoltageMatrix(0, 0);
      const Complex inverterInterfaceCurrent = inverterCurrentMatrix(0, 0);
      const Complex gridCurrent = inductorCurrent(0, 0);
      const Complex capacitorVoltage(inverterState(VcD, 0),
                                     inverterState(VcQ, 0));
      const Complex nonlinearGridCurrent =
          (capacitorVoltage - terminalVoltage) / mRc;
      Matrix simulationState = Matrix::Zero(StateCount, 1);
      simulationState.block(0, 0, 12, 1) = inverterState;
      simulationState(IGridD, 0) = gridCurrent.real();
      simulationState(IGridQ, 0) = gridCurrent.imag();

      result = {
          modalAnalysis.getDiscreteEigenvalues(),
          extractor.getDiscreteStateMatrix(),
          extractor.getMetadata().stateNames,
          simulationState,
          previousInverterState,
          terminalVoltage,
          inverterInterfaceCurrent,
          gridCurrent,
          nonlinearGridCurrent,
          handles.inverter->attributeTyped<Int>("iteration_count")->get(),
          handles.inverter
              ->attributeTyped<Real>("iteration_state_residual")
              ->get(),
          handles.inverter
              ->attributeTyped<Real>("iteration_input_residual")
              ->get(),
          extractor.getLastExtractionTime(),
          extractor.getStateCount(),
          handles.inverter->getFullNortonAdmittance(),
          handles.inverter->getComplexScalarNortonAdmittanceMatrix(),
          handles.inverter->getNortonAdmittanceStructureDefect(),
          handles.inverter->getDiscreteModelRevision(),
          handles.inverter->getMnaStampRevision(),
          handles.inverter->getExtractionStampRevision()};
      captured = true;
    }
    simulation.stop();
    if (!captured)
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
                                  "time step (" + std::to_string(timeStep) +
                                  " s).");
    return static_cast<UInt>(roundedSteps);
  }

  SteadyState calculateSteadyState() const {
    const Complex j(0.0, 1.0);
    const Complex powerReference(mPRef, mQRef);
    const Complex totalImpedance(mGridResistance + mRc,
                                 mOmega * mGridInductance);
    Complex filterVoltage = mSourceVoltage;
    Bool converged = false;
    for (Int iteration = 0; iteration < 100; ++iteration) {
      const Complex gridCurrent =
          std::conj(powerReference / filterVoltage);
      const Complex nextFilterVoltage =
          mSourceVoltage + totalImpedance * gridCurrent;
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

    const Complex gridCurrent =
        std::conj(powerReference / filterVoltage);
    const Complex pccVoltage = filterVoltage - mRc * gridCurrent;
    const Complex midVoltage =
        mSourceVoltage + mGridResistance * gridCurrent;
    const Complex filterCurrent =
        gridCurrent + j * mOmega * mCf * filterVoltage;
    const Complex converterVoltage =
        filterVoltage + (mRf + j * mOmega * mLf) * filterCurrent;
    const Real psi = std::arg(filterVoltage);
    const Complex rotation = std::exp(-j * psi);
    const Complex filterVoltageLocal = filterVoltage * rotation;
    const Complex gridCurrentLocal = gridCurrent * rotation;
    const Complex converterVoltageLocal = converterVoltage * rotation;
    const Real pInitial =
        filterVoltageLocal.real() * gridCurrentLocal.real() +
        filterVoltageLocal.imag() * gridCurrentLocal.imag();
    const Real qInitial =
        -filterVoltageLocal.real() * gridCurrentLocal.imag() +
        filterVoltageLocal.imag() * gridCurrentLocal.real();

    Matrix state = Matrix::Zero(StateCount, 1);
    state(Psi, 0) = psi;
    state(PhiPLL, 0) = 0.0;
    state(PFiltered, 0) = pInitial;
    state(QFiltered, 0) = qInitial;
    state(PhiD, 0) =
        (gridCurrentLocal.real() +
         mKpPowerCtrl * (pInitial - mPRef)) /
        mKiPowerCtrl;
    state(PhiQ, 0) =
        (gridCurrentLocal.imag() -
         mKpPowerCtrl * (qInitial - mQRef)) /
        mKiPowerCtrl;
    const Real currentReferenceD =
        mKpPowerCtrl * (mPRef - pInitial) + mKiPowerCtrl * state(PhiD, 0);
    const Real currentReferenceQ =
        mKpPowerCtrl * (qInitial - mQRef) + mKiPowerCtrl * state(PhiQ, 0);
    state(GammaD, 0) =
        (converterVoltageLocal.real() +
         mKpCurrCtrl * (gridCurrentLocal.real() - currentReferenceD)) /
        mKiCurrCtrl;
    state(GammaQ, 0) =
        (converterVoltageLocal.imag() +
         mKpCurrCtrl * (gridCurrentLocal.imag() - currentReferenceQ)) /
        mKiCurrCtrl;
    state(VcD, 0) = filterVoltage.real();
    state(VcQ, 0) = filterVoltage.imag();
    state(IfD, 0) = filterCurrent.real();
    state(IfQ, 0) = filterCurrent.imag();
    state(IGridD, 0) = gridCurrent.real();
    state(IGridQ, 0) = gridCurrent.imag();

    const SteadyState result{mSourceVoltage, midVoltage, pccVoltage, state};
    const Real residual = independentDqDerivative(result.state).norm();
    if (residual > 1e-6)
      throw std::runtime_error(
          "Manual operating point is not an equilibrium of the independent "
          "dq equations. Residual norm: " +
          std::to_string(residual));
    return result;
  }

  // The electrical states are expressed directly in the global synchronous
  // dq frame used by DP. Controller measurements are rotated into the PLL
  // frame. This model is independent of the simulator's assembled matrices.
  Matrix independentDqDerivative(const Matrix &state) const {
    const Real cosine = std::cos(state(Psi, 0));
    const Real sine = std::sin(state(Psi, 0));
    const Real vcDLocal =
        cosine * state(VcD, 0) + sine * state(VcQ, 0);
    const Real vcQLocal =
        -sine * state(VcD, 0) + cosine * state(VcQ, 0);
    const Real gridCurrentDLocal =
        cosine * state(IGridD, 0) + sine * state(IGridQ, 0);
    const Real gridCurrentQLocal =
        -sine * state(IGridD, 0) + cosine * state(IGridQ, 0);
    const Real pInstantaneous =
        vcDLocal * gridCurrentDLocal + vcQLocal * gridCurrentQLocal;
    const Real qInstantaneous =
        -vcDLocal * gridCurrentQLocal + vcQLocal * gridCurrentDLocal;
    const Real currentReferenceD =
        mKpPowerCtrl * (mPRef - state(PFiltered, 0)) +
        mKiPowerCtrl * state(PhiD, 0);
    const Real currentReferenceQ =
        mKpPowerCtrl * (state(QFiltered, 0) - mQRef) +
        mKiPowerCtrl * state(PhiQ, 0);
    const Real converterVoltageDLocal =
        mKpCurrCtrl * (currentReferenceD - gridCurrentDLocal) +
        mKiCurrCtrl * state(GammaD, 0);
    const Real converterVoltageQLocal =
        mKpCurrCtrl * (currentReferenceQ - gridCurrentQLocal) +
        mKiCurrCtrl * state(GammaQ, 0);
    const Real converterVoltageD =
        cosine * converterVoltageDLocal - sine * converterVoltageQLocal;
    const Real converterVoltageQ =
        sine * converterVoltageDLocal + cosine * converterVoltageQLocal;

    Matrix derivative = Matrix::Zero(StateCount, 1);
    derivative(Psi, 0) =
        mKpPLL * vcQLocal + mKiPLL * state(PhiPLL, 0);
    derivative(PhiPLL, 0) = vcQLocal;
    derivative(PFiltered, 0) =
        mOmegaCutoff * (pInstantaneous - state(PFiltered, 0));
    derivative(QFiltered, 0) =
        mOmegaCutoff * (qInstantaneous - state(QFiltered, 0));
    derivative(PhiD, 0) = mPRef - state(PFiltered, 0);
    derivative(PhiQ, 0) = state(QFiltered, 0) - mQRef;
    derivative(GammaD, 0) = currentReferenceD - gridCurrentDLocal;
    derivative(GammaQ, 0) = currentReferenceQ - gridCurrentQLocal;

    derivative(VcD, 0) =
        (state(IfD, 0) - state(IGridD, 0)) / mCf +
        mOmega * state(VcQ, 0);
    derivative(VcQ, 0) =
        (state(IfQ, 0) - state(IGridQ, 0)) / mCf -
        mOmega * state(VcD, 0);
    derivative(IfD, 0) =
        (converterVoltageD - state(VcD, 0) - mRf * state(IfD, 0)) /
            mLf +
        mOmega * state(IfQ, 0);
    derivative(IfQ, 0) =
        (converterVoltageQ - state(VcQ, 0) - mRf * state(IfQ, 0)) /
            mLf -
        mOmega * state(IfD, 0);
    derivative(IGridD, 0) =
        (state(VcD, 0) - mSourceVoltage.real() -
         (mGridResistance + mRc) * state(IGridD, 0)) /
            mGridInductance +
        mOmega * state(IGridQ, 0);
    derivative(IGridQ, 0) =
        (state(VcQ, 0) - mSourceVoltage.imag() -
         (mGridResistance + mRc) * state(IGridQ, 0)) /
            mGridInductance -
        mOmega * state(IGridD, 0);
    return derivative;
  }

  Matrix buildIndependentDqJacobian(const Matrix &operatingPoint) const {
    return numericalJacobian(
        [this](const Matrix &state) { return independentDqDerivative(state); },
        operatingPoint);
  }

  Real mFrequency;
  Real mOmega;
  Complex mSourceVoltage;
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
    DPPh1InverterStateSpaceExtractionExample example(extractionTime,
                                                      finalTime);
    example.run();
    return 0;
  } catch (const std::exception &exception) {
    std::cerr << "Error: " << exception.what() << '\n';
    return 1;
  }
}
