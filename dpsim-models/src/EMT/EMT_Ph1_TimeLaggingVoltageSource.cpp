#include <dpsim-models/EMT/EMT_Ph1_TimeLaggingVoltageSource.h>

using namespace CPS;

EMT::Ph1::TimeLaggingVoltageSource::TimeLaggingVoltageSource(
    String uid, String name, Logger::Level logLevel)
    : MNASimPowerComp<Real>(uid, name, true, true, logLevel),
      mVoltage(mAttributes->create<Real>("v")) {
  setVirtualNodeNumber(1);
  setTerminalNumber(2);
  **mIntfVoltage = Matrix::Zero(1, 1);
  **mIntfCurrent = Matrix::Zero(1, 1);
}

void EMT::Ph1::TimeLaggingVoltageSource::mnaCompInitialize(
    Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
  updateMatrixNodeIndices();
  mP = matrixNodeIndex(0);
  mQ = matrixNodeIndex(1);
  mIpq = mVirtualNodes[0]->matrixNodeIndex();
}

void EMT::Ph1::TimeLaggingVoltageSource::mnaCompApplySystemMatrixStamp(
    SparseMatrixRow &systemMatrix) {
  if (terminalNotGrounded(0)) {
    Math::addToMatrixElement(systemMatrix, mP, mIpq, 1);
    Math::addToMatrixElement(systemMatrix, mIpq, mP, 1);
  }
  if (terminalNotGrounded(1)) {
    Math::addToMatrixElement(systemMatrix, mQ, mIpq, -1);
    Math::addToMatrixElement(systemMatrix, mIpq, mQ, -1);
  }
}

void EMT::Ph1::TimeLaggingVoltageSource::mnaCompAddPreStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes) {
  prevStepDependencies.push_back(mRightVector);
  modifiedAttributes.push_back(mRightVector);
  modifiedAttributes.push_back(mVoltage);
}

void EMT::Ph1::TimeLaggingVoltageSource::mnaCompPreStep(Real time,
                                                        Int timeStepCount) {
  if (mNodeM == nullptr || mNodeN == nullptr) {
    throw std::runtime_error(
        "mNodeM or mNodeN is null in TimeLaggingVoltageSource.");
  }
  Real vM = mNodeM->voltage()(0, 0);
  Real vN = mNodeN->voltage()(0, 0);
  (**mVoltage) = vM - vN;
  Math::setVectorElement(**mRightVector, mIpq, **mVoltage);
}

void EMT::Ph1::TimeLaggingVoltageSource::stampBranchNodeIncidenceMatrix(
    UInt branchIdx, Matrix &branchNodeIncidenceMatrix) {
  if (terminalNotGrounded(0)) {
    branchNodeIncidenceMatrix(branchIdx, mP) = 1.0;
  }
  if (terminalNotGrounded(1)) {
    branchNodeIncidenceMatrix(branchIdx, mQ) = -1.0;
  }
}
