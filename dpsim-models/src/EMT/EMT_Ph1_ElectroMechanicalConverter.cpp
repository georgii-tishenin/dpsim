#include <dpsim-models/EMT/EMT_Ph1_ElectroMechanicalConverter.h>
#include <dpsim-models/EMT/EMT_Ph1_Inductor.h>

using namespace CPS;

EMT::Ph1::ElectroMechanicalConverter::ElectroMechanicalConverter(
    String uid, String name, Logger::Level logLevel)
    : MNASimPowerComp<Real>(uid, name, true, true, logLevel),
    mFlux(mAttributes->create<Real>("flux")) {

  setVirtualNodeNumber(1);

  setTerminalNumber(4);

  SPDLOG_LOGGER_INFO(mSLog, "Create {} {}", this->type(), name);
  **mIntfVoltage = Matrix::Zero(1, 1);
  **mIntfCurrent = Matrix::Zero(1, 1);
}

void EMT::Ph1::ElectroMechanicalConverter::setInitialFlux(Real flux) {
  **mFlux = flux;
}

void EMT::Ph1::ElectroMechanicalConverter::initializeFromNodesAndTerminals(
    Real frequency) {
}

void EMT::Ph1::ElectroMechanicalConverter::mnaCompInitialize(
    Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
  updateMatrixNodeIndices();
  mTimeStep = timeStep;
}

void EMT::Ph1::ElectroMechanicalConverter::mnaCompApplySystemMatrixStamp(
    SparseMatrixRow &systemMatrix) {

  // Ideal transformer equations
  if (terminalNotGrounded(0)) {
    Math::setMatrixElement(systemMatrix, mVirtualNodes[0]->matrixNodeIndex(),
                           matrixNodeIndex(0), 1 / (**mFlux));
    Math::setMatrixElement(systemMatrix, matrixNodeIndex(0),
                           mVirtualNodes[0]->matrixNodeIndex(), 1 / (**mFlux));
  }
  if (terminalNotGrounded(1)) {
    Math::setMatrixElement(systemMatrix, mVirtualNodes[0]->matrixNodeIndex(),
                           matrixNodeIndex(1), -1 / (**mFlux));
    Math::setMatrixElement(systemMatrix, matrixNodeIndex(1),
                           mVirtualNodes[0]->matrixNodeIndex(), -1 / (**mFlux));
  }
  if (terminalNotGrounded(2)) {
    Math::setMatrixElement(systemMatrix, mVirtualNodes[0]->matrixNodeIndex(),
                           matrixNodeIndex(2), -1);
    Math::setMatrixElement(systemMatrix, matrixNodeIndex(2),
                           mVirtualNodes[0]->matrixNodeIndex(), -1);
  }
  if (terminalNotGrounded(3)) {
    Math::setMatrixElement(systemMatrix, mVirtualNodes[0]->matrixNodeIndex(),
                           matrixNodeIndex(3), 1);
    Math::setMatrixElement(systemMatrix, matrixNodeIndex(3),
                           mVirtualNodes[0]->matrixNodeIndex(), 1);
  }

  if (terminalNotGrounded(0)) {
    SPDLOG_LOGGER_INFO(mSLog, "Add {:s} to system at ({:d},{:d})",
                       Logger::complexToString(Complex(1.0 / (**mFlux), 0)),
                       mVirtualNodes[0]->matrixNodeIndex(), matrixNodeIndex(0));
    SPDLOG_LOGGER_INFO(mSLog, "Add {:s} to system at ({:d},{:d})",
                       Logger::complexToString(Complex(1.0 / (**mFlux), 0)),
                       matrixNodeIndex(0), mVirtualNodes[0]->matrixNodeIndex());
  }
  if (terminalNotGrounded(1)) {
    SPDLOG_LOGGER_INFO(mSLog, "Add {:s} to system at ({:d},{:d})",
                       Logger::complexToString(Complex(-1.0 / (**mFlux), 0)),
                       mVirtualNodes[0]->matrixNodeIndex(), matrixNodeIndex(1));
    SPDLOG_LOGGER_INFO(mSLog, "Add {:s} to system at ({:d},{:d})",
                       Logger::complexToString(Complex(-1.0 / (**mFlux), 0)),
                       matrixNodeIndex(1), mVirtualNodes[0]->matrixNodeIndex());
  }
  if (terminalNotGrounded(2)) {
    SPDLOG_LOGGER_INFO(mSLog, "Add {:s} to system at ({:d},{:d})",
                       Logger::complexToString(-1.0),
                       mVirtualNodes[0]->matrixNodeIndex(), matrixNodeIndex(2));
    SPDLOG_LOGGER_INFO(mSLog, "Add {:s} to system at ({:d},{:d})",
                       Logger::complexToString(-1.0), matrixNodeIndex(2),
                       mVirtualNodes[0]->matrixNodeIndex());
  }
  if (terminalNotGrounded(3)) {
    SPDLOG_LOGGER_INFO(mSLog, "Add {:s} to system at ({:d},{:d})",
                       Logger::complexToString(-1.0),
                       mVirtualNodes[0]->matrixNodeIndex(), matrixNodeIndex(3));
    SPDLOG_LOGGER_INFO(mSLog, "Add {:s} to system at ({:d},{:d})",
                       Logger::complexToString(-1.0), matrixNodeIndex(3),
                       mVirtualNodes[0]->matrixNodeIndex());
  }
}

void EMT::Ph1::ElectroMechanicalConverter::mnaCompAddPostStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes,
    Attribute<Matrix>::Ptr &leftVector) {
  attributeDependencies.push_back(leftVector);
  modifiedAttributes.push_back(mIntfVoltage);
  modifiedAttributes.push_back(mIntfCurrent);
}

void EMT::Ph1::ElectroMechanicalConverter::mnaCompPostStep(
    Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector) {
  mnaCompUpdateVoltage(**leftVector);
  mnaCompUpdateCurrent(**leftVector);

  if (mInductor != nullptr) {

    mOldVoltage = mVoltage;
    // Update the actual voltage of the terminal
    auto idx = mInductor->matrixNodeIndex(0);
    mVoltage = Math::realFromVectorElement(**leftVector, idx);
  }
}

void EMT::Ph1::ElectroMechanicalConverter::mnaCompUpdateVoltage(
    const Matrix &leftVector) {
}

void EMT::Ph1::ElectroMechanicalConverter::mnaCompUpdateCurrent(
    const Matrix &leftVector) {}

void EMT::Ph1::ElectroMechanicalConverter::stampBranchNodeIncidenceMatrix(
    UInt branchIdx, Matrix &branchNodeIncidenceMatrix) {
}

void EMT::Ph1::ElectroMechanicalConverter::mnaCompPreStep(
    Real time, Int timeStepCount) {
  **mFlux = **mFlux + (mTimeStep / 2) * (mOldVoltage + mVoltage);
  mnaCompApplyRightSideVectorStamp(**mRightVector);
}

void EMT::Ph1::ElectroMechanicalConverter::mnaCompAddPreStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes) {
  // add pre-step dependencies of component itself
  prevStepDependencies.push_back(mIntfCurrent);
  prevStepDependencies.push_back(mIntfVoltage);
  modifiedAttributes.push_back(mRightVector);
}