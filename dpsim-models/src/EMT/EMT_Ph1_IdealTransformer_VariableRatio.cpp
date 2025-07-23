#include <dpsim-models/EMT/EMT_Ph1_IdealTransformer_VariableRatio.h>
#include <dpsim-models/EMT/EMT_Ph1_Inductor.h>

using namespace CPS;

EMT::Ph1::IdealTransformerVariableRatio::IdealTransformerVariableRatio(
    String uid, String name, Logger::Level logLevel)
    : MNASimPowerComp<Real>(uid, name, true, true, logLevel),
      mN(std::make_shared<Real>()), mRatio(mAttributes->create<Real>("Ratio")) {

  setVirtualNodeNumber(1);

  setTerminalNumber(4);

  SPDLOG_LOGGER_INFO(mSLog, "Create {} {}", this->type(), name);
  **mIntfVoltage = Matrix::Zero(1, 1);
  **mIntfCurrent = Matrix::Zero(1, 1);
}

void EMT::Ph1::IdealTransformerVariableRatio::setParameters(Real N) {

  *mN = N;

  SPDLOG_LOGGER_INFO(mSLog, "Turns Ratio={} [ ] ", std::abs(N));

  mParametersSet = true;
}

void EMT::Ph1::IdealTransformerVariableRatio::initializeFromNodesAndTerminals(
    Real frequency) {
}

void EMT::Ph1::IdealTransformerVariableRatio::mnaCompInitialize(
    Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
  updateMatrixNodeIndices();
  mTimeStep = timeStep;
}

void EMT::Ph1::IdealTransformerVariableRatio::mnaCompApplySystemMatrixStamp(
    SparseMatrixRow &systemMatrix) {

  // Ideal transformer equations
  if (terminalNotGrounded(0)) {
    Math::setMatrixElement(systemMatrix, mVirtualNodes[0]->matrixNodeIndex(),
                           matrixNodeIndex(0), 1 / (*mN));
    Math::setMatrixElement(systemMatrix, matrixNodeIndex(0),
                           mVirtualNodes[0]->matrixNodeIndex(), 1 / (*mN));
  }
  if (terminalNotGrounded(1)) {
    Math::setMatrixElement(systemMatrix, mVirtualNodes[0]->matrixNodeIndex(),
                           matrixNodeIndex(1), -1 / (*mN));
    Math::setMatrixElement(systemMatrix, matrixNodeIndex(1),
                           mVirtualNodes[0]->matrixNodeIndex(), -1 / (*mN));
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
                       Logger::complexToString(Complex(1.0 / (*mN), 0)),
                       mVirtualNodes[0]->matrixNodeIndex(), matrixNodeIndex(0));
    SPDLOG_LOGGER_INFO(mSLog, "Add {:s} to system at ({:d},{:d})",
                       Logger::complexToString(Complex(1.0 / (*mN), 0)),
                       matrixNodeIndex(0), mVirtualNodes[0]->matrixNodeIndex());
  }
  if (terminalNotGrounded(1)) {
    SPDLOG_LOGGER_INFO(mSLog, "Add {:s} to system at ({:d},{:d})",
                       Logger::complexToString(Complex(-1.0 / (*mN), 0)),
                       mVirtualNodes[0]->matrixNodeIndex(), matrixNodeIndex(1));
    SPDLOG_LOGGER_INFO(mSLog, "Add {:s} to system at ({:d},{:d})",
                       Logger::complexToString(Complex(-1.0 / (*mN), 0)),
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

void EMT::Ph1::IdealTransformerVariableRatio::mnaCompAddPostStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes,
    Attribute<Matrix>::Ptr &leftVector) {
  attributeDependencies.push_back(leftVector);
  modifiedAttributes.push_back(mIntfVoltage);
  modifiedAttributes.push_back(mIntfCurrent);
}

void EMT::Ph1::IdealTransformerVariableRatio::mnaCompPostStep(
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

void EMT::Ph1::IdealTransformerVariableRatio::mnaCompUpdateVoltage(
    const Matrix &leftVector) {
}

void EMT::Ph1::IdealTransformerVariableRatio::mnaCompUpdateCurrent(
    const Matrix &leftVector) {}

void EMT::Ph1::IdealTransformerVariableRatio::stampBranchNodeIncidenceMatrix(
    UInt branchIdx, Matrix &branchNodeIncidenceMatrix) {
}

void EMT::Ph1::IdealTransformerVariableRatio::mnaCompPreStep(
    Real time, Int timeStepCount) {
  *mN = *mN + (mTimeStep / 2) * (mOldVoltage + mVoltage);
  mRatio->set(*mN);
  mnaCompApplyRightSideVectorStamp(**mRightVector);
}

void EMT::Ph1::IdealTransformerVariableRatio::mnaCompAddPreStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes) {
  // add pre-step dependencies of component itself
  prevStepDependencies.push_back(mIntfCurrent);
  prevStepDependencies.push_back(mIntfVoltage);
  modifiedAttributes.push_back(mRightVector);
}