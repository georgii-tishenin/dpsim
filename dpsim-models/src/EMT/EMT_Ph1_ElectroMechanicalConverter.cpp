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

void EMT::Ph1::ElectroMechanicalConverter::mnaCompInitialize(
    Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
  updateMatrixNodeIndices();
  mTimeStep = timeStep;
}

void EMT::Ph1::ElectroMechanicalConverter::mnaCompApplySystemMatrixStamp(
    SparseMatrixRow &systemMatrix) {
  if (**mFlux == 0) {
    **mFlux = 1e-12; // Avoid division by zero
  }

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
}

void EMT::Ph1::ElectroMechanicalConverter::mnaCompAddPostStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes,
    Attribute<Matrix>::Ptr &leftVector) {
  prevStepDependencies.push_back(mFlux);
  attributeDependencies.push_back(leftVector);
  modifiedAttributes.push_back(mFlux);
}

void EMT::Ph1::ElectroMechanicalConverter::mnaCompPostStep(
    Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector) {
  if (mVoltageReferenceNode == nullptr) {
    throw std::runtime_error(
        "mFluxNode is null in ElectroMechanicalConverter.");
  }
  mOldVoltage = mVoltage;
  mVoltage = mVoltageReferenceNode->voltage()(0, 0);
  **mFlux = **mFlux + (mTimeStep / 2) * (mOldVoltage + mVoltage);
}

void EMT::Ph1::ElectroMechanicalConverter::stampBranchNodeIncidenceMatrix(
    UInt branchIdx, Matrix &branchNodeIncidenceMatrix) {}