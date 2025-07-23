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

  // Component parameters are referred to higher voltage side.
  // Switch terminals to have terminal 0 at higher voltage side
  // if transformer is connected the other way around.
  if (Math::abs(*mN) < 1.) {
    if (*mN == 0) {
      *mN = 0.00001;
    } else {
      *mN = 1. / *mN;
      std::shared_ptr<SimTerminal<Real>> tmp = mTerminals[0];
      mTerminals[0] = mTerminals[1];
      mTerminals[1] = tmp;
      SPDLOG_LOGGER_INFO(mSLog, "Switching terminals to have first terminal at "
                                "higher voltage side. Updated parameters: ");
      SPDLOG_LOGGER_INFO(mSLog, "Turns Ratio = {} [ ]", std::abs(*mN));
    }
  }

  // Set initial voltage of virtual node in between
  mVirtualNodes[0]->setInitialVoltage(initialSingleVoltage(1) * *mN);

  // Log initialization results
  SPDLOG_LOGGER_INFO(
      mSLog,
      "--- Initialization ---\n"
      "Terminal 0 voltage: {:s}\n"
      "Terminal 1 voltage: {:s}\n"
      "Virtual Node voltage: {:s}\n"
      "--- Initialization finished ---",
      Logger::phasorToString(initialSingleVoltage(0)),
      Logger::phasorToString(initialSingleVoltage(1)),
      Logger::phasorToString(mVirtualNodes[0]->initialSingleVoltage()));
}

void EMT::Ph1::IdealTransformerVariableRatio::mnaCompInitialize(
    Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
  updateMatrixNodeIndices();
  //**mRightVector = Matrix::Zero(0, 0); //If you uncomment this line, mnaCompPreStep will not be called every time step.
}

void EMT::Ph1::IdealTransformerVariableRatio::mnaCompApplySystemMatrixStamp(
    SparseMatrixRow &systemMatrix) {

  //    *mN = *mN + ( mTimeStep / 2 ) * (mOldVoltage + mVoltage);
  //    mRatio->set(*mN);

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

  //updateVoltages(**leftVector);

  if (mInductor != nullptr) {

    mOldVoltage = mVoltage;
    // Update the actual voltage of the terminal
    //auto idx = mInductor->terminal(0)->matrixNodeIndex();
    auto idx = mInductor->matrixNodeIndex(0);
    mVoltage = Math::realFromVectorElement(**leftVector, idx);
  }
}

void EMT::Ph1::IdealTransformerVariableRatio::mnaCompUpdateVoltage(
    const Matrix &leftVector) {
  // v1 - v0
  /*
  (**mIntfVoltage)(0, 0) = 0;
  (**mIntfVoltage)(0, 0) =
      Math::realFromVectorElement(leftVector, matrixNodeIndex(1));
  (**mIntfVoltage)(0, 0) = (**mIntfVoltage)(0, 0) -
                           Math::realFromVectorElement(
                               leftVector, mVirtualNodes[0]->matrixNodeIndex());
  SPDLOG_LOGGER_DEBUG(mSLog, "Voltage {:s}",
                      Logger::phasorToString((**mIntfVoltage)(0, 0)));
                      */
}

void EMT::Ph1::IdealTransformerVariableRatio::mnaCompUpdateCurrent(
    const Matrix &leftVector) {}

void EMT::Ph1::IdealTransformerVariableRatio::stampBranchNodeIncidenceMatrix(
    UInt branchIdx, Matrix &branchNodeIncidenceMatrix) {

  /*
  UInt primaryBranchIdx = branchIdx - 1;
  UInt secondaryBranchIdx = branchIdx;

  if (terminalNotGrounded(0)) {
    branchNodeIncidenceMatrix(primaryBranchIdx, matrixNodeIndex(0)) = 1.0;
  }
    
    branchNodeIncidenceMatrix(primaryBranchIdx, mVirtualNodes[0]->matrixNodeIndex()) = -1.0;
  
  if (terminalNotGrounded(1)) {
  branchNodeIncidenceMatrix(secondaryBranchIdx, matrixNodeIndex(1)) = -1;
 }

  branchNodeIncidenceMatrix(secondaryBranchIdx, mVirtualNodes[0]->matrixNodeIndex()) = 1;
  */
}

void EMT::Ph1::IdealTransformerVariableRatio::setTimeStep(Real timeStep) {
  mTimeStep = timeStep;
}

/*
 // A function which serves to update the voltages which will be used in the mnaCompPreStep to calculate the new Ratio as a discrete integration
 // of the voltage source declared as a pointer in the header file
 void EMT::Ph1::IdealTransformerVariableRatio::updateVoltages(const Matrix &leftVector) {
   if (mInductor != nullptr) {

    mOldVoltage = mVoltage;
    // Update the actual voltage of the terminal
    auto idx = mInductor->terminal(0)->matrixNodeIndex();
    mVoltage = Math::realFromVectorElement(**leftVector, idx);
    //mVoltage = mVoltageSource->mVoltageRef->get().real();
   }
 }
*/

void EMT::Ph1::IdealTransformerVariableRatio::mnaCompPreStep(
    Real time, Int timeStepCount) {

  if (mNegativeSpeedTermVoltageFlag == true) {
    *mN = *mN - (mTimeStep / 2) * (mOldVoltage + mVoltage);
  } else {
    // If the voltage is not negative, we can use the voltage as it is
    *mN = *mN + (mTimeStep / 2) * (mOldVoltage + mVoltage);
  }

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