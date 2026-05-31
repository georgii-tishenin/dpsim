#include <dpsim-models/EMT/EMT_Ph1_CurrentControlledVoltageSource.h>

using namespace CPS;

EMT::Ph1::CurrentControlledVoltageSource::CurrentControlledVoltageSource(
    String uid, String name, Logger::Level logLevel)
    : MNASimPowerComp<Real>(uid, name, true, true, logLevel) {

  setVirtualNodeNumber(2);
  setTerminalNumber(4);

  SPDLOG_LOGGER_INFO(mSLog, "Create {} {}", this->type(), name);
  **mIntfVoltage = Matrix::Zero(1, 1);
  **mIntfCurrent = Matrix::Zero(1, 1);
}



void EMT::Ph1::CurrentControlledVoltageSource::mnaCompInitialize(
    Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
  updateMatrixNodeIndices();
}

void EMT::Ph1::CurrentControlledVoltageSource::mnaCompApplySystemMatrixStamp(
    SparseMatrixRow &systemMatrix) { 

  int m = matrixNodeIndex(0);
  int n = matrixNodeIndex(1);
  int p = matrixNodeIndex(2);
  int q = matrixNodeIndex(3);

  int imn = mVirtualNodes[0]->matrixNodeIndex();
  int ipq = mVirtualNodes[1]->matrixNodeIndex();

  // current controlled voltage source equations
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
    Math::setMatrixElement(systemMatrix, p, ipq, 1);
    Math::setMatrixElement(systemMatrix, ipq, p, 1);
  }
  //q
  if (terminalNotGrounded(3)) {
    Math::setMatrixElement(systemMatrix, q, ipq, -1);
    Math::setMatrixElement(systemMatrix, ipq, q, -1);
  }

  double alpha = mCoefficient;

  Math::setMatrixElement(systemMatrix, ipq, imn, -alpha);
}

void EMT::Ph1::CurrentControlledVoltageSource::mnaCompAddPostStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes,
    Attribute<Matrix>::Ptr &leftVector) {
  attributeDependencies.push_back(leftVector);
}

void EMT::Ph1::CurrentControlledVoltageSource::mnaCompPostStep(
    Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector) {
  if (mVoltageReferenceNode == nullptr) {
    throw std::runtime_error(
        "Voltage Reference node is null in CurrentControlledVoltageSource.");
  }

}

void EMT::Ph1::CurrentControlledVoltageSource::stampBranchNodeIncidenceMatrix(
    UInt branchIdx, Matrix &branchNodeIncidenceMatrix) {}