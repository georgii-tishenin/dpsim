#include <dpsim-models/EMT/EMT_Ph1_CurrentControlledTorqueSource.h>

using namespace CPS;

EMT::Ph1::CurrentControlledTorqueSource::CurrentControlledTorqueSource(
    String uid, String name, Logger::Level logLevel)
    : MNASimPowerComp<Real>(uid, name, true, true, logLevel),
      mCoefficientCurrent(mAttributes->create<Real>("coeffI")){

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

  double alpha = mCoefficient * **mCoefficientCurrent;
  if (mIsNegative) {
    alpha = -alpha;
  }

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

void EMT::Ph1::CurrentControlledTorqueSource::mnaCompAddPostStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes,
    Attribute<Matrix>::Ptr &leftVector) {

  // TODO
}

void EMT::Ph1::CurrentControlledTorqueSource::mnaCompPostStep(
    Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector) {
  // TODO
}

void EMT::Ph1::CurrentControlledTorqueSource::stampBranchNodeIncidenceMatrix(
    UInt branchIdx, Matrix &branchNodeIncidenceMatrix) {}