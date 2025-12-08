#include <dpsim-models/EMT/EMT_Ph1_TimeLaggingTorqueSource.h>

using namespace CPS;

EMT::Ph1::TimeLaggingTorqueSource::TimeLaggingTorqueSource(
    String uid, String name, Logger::Level logLevel)
    : MNASimPowerComp<Real>(uid, name, true, true, logLevel),
      mCurrent(mAttributes->create<Real>("i")) {
  setVirtualNodeNumber(1);
  setTerminalNumber(4);
  **mIntfVoltage = Matrix::Zero(1, 1);
  **mIntfCurrent = Matrix::Zero(1, 1);
}

void EMT::Ph1::TimeLaggingTorqueSource::mnaCompInitialize(
    Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
  updateMatrixNodeIndices();
  mM = matrixNodeIndex(0);
  mN = matrixNodeIndex(1);
  mP = matrixNodeIndex(2);
  mQ = matrixNodeIndex(3);
  mImn = mVirtualNodes[0]->matrixNodeIndex();
}

void EMT::Ph1::TimeLaggingTorqueSource::mnaCompApplySystemMatrixStamp(
    SparseMatrixRow &systemMatrix) {
  if (terminalNotGrounded(0)) {
    Math::addToMatrixElement(systemMatrix, mM, mImn, 1);
    Math::addToMatrixElement(systemMatrix, mImn, mM, 1);
  }
  if (terminalNotGrounded(1)) {
    Math::addToMatrixElement(systemMatrix, mN, mImn, -1);
    Math::addToMatrixElement(systemMatrix, mImn, mN, -1);
  }
}

void EMT::Ph1::TimeLaggingTorqueSource::mnaCompAddPreStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes) {
  prevStepDependencies.push_back(mCurrent);
  modifiedAttributes.push_back(mRightVector);
}

void EMT::Ph1::TimeLaggingTorqueSource::mnaCompAddPostStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes,
    Attribute<Matrix>::Ptr &leftVector) {
  attributeDependencies.push_back(leftVector);
  modifiedAttributes.push_back(mCurrent);
}

void EMT::Ph1::TimeLaggingTorqueSource::mnaCompPreStep(Real time,
                                                       Int timeStepCount) {
  if (terminalNotGrounded(2)) {
    Math::setVectorElement(**mRightVector, mP, -**mCurrent);
  }
  if (terminalNotGrounded(3)) {
    Math::setVectorElement(**mRightVector, mQ, **mCurrent);
  }
}

void EMT::Ph1::TimeLaggingTorqueSource::mnaCompPostStep(
    Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector) {
  (**mCurrent) = Math::realFromVectorElement(**leftVector, mImn);
}

void EMT::Ph1::TimeLaggingTorqueSource::stampBranchNodeIncidenceMatrix(
    UInt branchIdx, Matrix &branchNodeIncidenceMatrix) {
  if (terminalNotGrounded(0)) {
    branchNodeIncidenceMatrix(branchIdx, mM) = 1.0;
  }
  if (terminalNotGrounded(1)) {
    branchNodeIncidenceMatrix(branchIdx, mN) = -1.0;
  }
}
