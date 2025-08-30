#include <dpsim-models/EMT/EMT_Ph1_CurrentControlledTorqueSource.h>

using namespace CPS;

EMT::Ph1::CurrentControlledTorqueSource::CurrentControlledTorqueSource(
    String uid, String name, Logger::Level logLevel)
    : MNASimPowerComp<Real>(uid, name, true, true, logLevel),
      mCoefficientCurrent(mAttributes->create<Real>("coeffI")) {

  setVirtualNodeNumber(1);
  setTerminalNumber(4);

  SPDLOG_LOGGER_INFO(mSLog, "Create {} {}", this->type(), name);
  **mIntfVoltage = Matrix::Zero(1, 1);
  **mIntfCurrent = Matrix::Zero(1, 1);
}

void EMT::Ph1::CurrentControlledTorqueSource::mnaCompInitialize(
    Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
  updateMatrixNodeIndices();
}

void EMT::Ph1::CurrentControlledTorqueSource::mnaCompApplySystemMatrixStamp(
    SparseMatrixRow &systemMatrix) {
  int m = matrixNodeIndex(0);
  int n = matrixNodeIndex(1);
  int p = matrixNodeIndex(2);
  int q = matrixNodeIndex(3);

  int imn = mVirtualNodes[0]->matrixNodeIndex();

  double alpha = calculateAlpha();

  // current controlled torque source equations
  //m
  if (terminalNotGrounded(0)) {
    Math::setMatrixElement(systemMatrix, m, imn, 1);
    Math::setMatrixElement(systemMatrix, imn, m, 1);
  }
  //n
  if (terminalNotGrounded(1)) {
    Math::setMatrixElement(systemMatrix, n, imn, -1);
    Math::setMatrixElement(systemMatrix, imn, n, -1);
  }
  //p
  if (terminalNotGrounded(2)) {
    Math::setMatrixElement(systemMatrix, p, imn, alpha);
  }
  //q
  if (terminalNotGrounded(3)) {
    Math::setMatrixElement(systemMatrix, q, imn, -alpha);
  }
}

double EMT::Ph1::CurrentControlledTorqueSource::calculateAlpha() {
  updateCoefficientCurrent();
  double alpha = mCoefficient * **mCoefficientCurrent;
  if (mIsNegative) {
    alpha = -alpha;
  }
  return alpha;
}

void EMT::Ph1::CurrentControlledTorqueSource::updateCoefficientCurrent() {
  if (mInductor == nullptr) {
    SPDLOG_LOGGER_ERROR(mSLog, "CurrentControlledTorqueSource has no inductor "
                               "set for current reference");
    throw std::runtime_error("No inductor set for current reference");
  }
  auto inductorCurrent = mInductor->intfCurrent();
  **mCoefficientCurrent = inductorCurrent(0, 0);
  **mCoefficientCurrent *=
      -1; // interface current in inductor is from terminal 1 to terminal 0
}

void EMT::Ph1::CurrentControlledTorqueSource::stampBranchNodeIncidenceMatrix(
    UInt branchIdx, Matrix &branchNodeIncidenceMatrix) {}