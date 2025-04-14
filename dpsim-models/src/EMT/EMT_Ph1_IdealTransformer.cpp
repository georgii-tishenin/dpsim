#include <dpsim-models/EMT/EMT_Ph1_IdealTransformer.h>

using namespace CPS;



EMT::Ph1::IdealTransformer::IdealTransformer(String uid, String name,
                                  Logger::Level logLevel)
      : MNASimPowerComp<Real>(uid, name, true, true, logLevel), 
       mN(std::make_shared<Real>()) {

  setVirtualNodeNumber(1);

  setTerminalNumber(2);

  SPDLOG_LOGGER_INFO(mSLog, "Create {} {}", this->type(), name);
  **mIntfVoltage = Matrix::Zero(1, 1);
  **mIntfCurrent = Matrix::Zero(1, 1);
}



void EMT::Ph1::IdealTransformer::setParameters(Real N) {

 *mN = N;

  SPDLOG_LOGGER_INFO(mSLog, "Turns Ratio={} [ ] ",
                     std::abs(N));

  mParametersSet = true;
  }


  void EMT::Ph1::IdealTransformer::initializeFromNodesAndTerminals(Real frequency) {

  // Component parameters are referred to higher voltage side.
  // Switch terminals to have terminal 0 at higher voltage side
  // if transformer is connected the other way around.
  if (Math::abs(*mN) < 1.) {
    *mN = 1. / *mN;
    std::shared_ptr<SimTerminal<Real>> tmp = mTerminals[0];
    mTerminals[0] = mTerminals[1];
    mTerminals[1] = tmp;
    SPDLOG_LOGGER_INFO(mSLog, "Switching terminals to have first terminal at "
                              "higher voltage side. Updated parameters: ");
    SPDLOG_LOGGER_INFO(mSLog, "Turns Ratio = {} [ ]",
                       std::abs(*mN));
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


void EMT::Ph1::IdealTransformer::mnaCompInitialize(Real omega, Real timeStep,
                                           Attribute<Matrix>::Ptr leftVector) {
  updateMatrixNodeIndices();
  **mRightVector = Matrix::Zero(0, 0);
}

void EMT::Ph1::IdealTransformer::mnaCompApplySystemMatrixStamp(
    SparseMatrixRow &systemMatrix) {
  // Ideal transformer equations
  if (terminalNotGrounded(0)) {
    Math::setMatrixElement(systemMatrix, mVirtualNodes[0]->matrixNodeIndex(),
                           matrixNodeIndex(0), 1/(*mN));
    Math::setMatrixElement(systemMatrix, matrixNodeIndex(0),
                           mVirtualNodes[0]->matrixNodeIndex(),
                           1/(*mN));
  }
  if (terminalNotGrounded(1)) {
    Math::setMatrixElement(systemMatrix, mVirtualNodes[0]->matrixNodeIndex(),
                           matrixNodeIndex(1), -1);
    Math::setMatrixElement(systemMatrix, matrixNodeIndex(1),
                           mVirtualNodes[0]->matrixNodeIndex(), -1);
  }


  if (terminalNotGrounded(0)) {
    SPDLOG_LOGGER_INFO(mSLog, "Add {:s} to system at ({:d},{:d})",
                       Logger::complexToString(Complex(-1.0, 0)),
                       mVirtualNodes[0]->matrixNodeIndex(),
                       matrixNodeIndex(0));
    SPDLOG_LOGGER_INFO(mSLog, "Add {:s} to system at ({:d},{:d})",
                       Logger::complexToString(Complex(-1.0, 0)),
                       matrixNodeIndex(0),
                       mVirtualNodes[0]->matrixNodeIndex());
  }
  if (terminalNotGrounded(1)) {
    SPDLOG_LOGGER_INFO(mSLog, "Add {:s} to system at ({:d},{:d})",
                       Logger::complexToString(*mN), mVirtualNodes[0]->matrixNodeIndex(),
                       matrixNodeIndex(1));
    SPDLOG_LOGGER_INFO(mSLog, "Add {:s} to system at ({:d},{:d})",
                       Logger::complexToString(*mN),
                       matrixNodeIndex(1), mVirtualNodes[0]->matrixNodeIndex());
  }
}


void EMT::Ph1::IdealTransformer::mnaCompAddPostStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes,
    Attribute<Matrix>::Ptr &leftVector) {
  attributeDependencies.push_back(leftVector);
  modifiedAttributes.push_back(mIntfVoltage);
  modifiedAttributes.push_back(mIntfCurrent);
}


void EMT::Ph1::IdealTransformer::mnaCompPostStep(Real time, Int timeStepCount,
                                         Attribute<Matrix>::Ptr &leftVector) {
  mnaCompUpdateVoltage(**leftVector);
  mnaCompUpdateCurrent(**leftVector);
}


void EMT::Ph1::IdealTransformer::mnaCompUpdateVoltage(const Matrix &leftVector) {
  // v1 - v0
  (**mIntfVoltage)(0, 0) = 0;
  (**mIntfVoltage)(0, 0) =
      Math::realFromVectorElement(leftVector, matrixNodeIndex(1));
  (**mIntfVoltage)(0, 0) = (**mIntfVoltage)(0, 0) -
                           Math::realFromVectorElement(
                               leftVector, mVirtualNodes[0]->matrixNodeIndex());
  SPDLOG_LOGGER_DEBUG(mSLog, "Voltage {:s}",
                      Logger::phasorToString((**mIntfVoltage)(0, 0)));
}


void EMT::Ph1::IdealTransformer::mnaCompUpdateCurrent(const Matrix &leftVector) {
}


void EMT::Ph1::IdealTransformer::stampBranchNodeIncidenceMatrix(UInt branchIdx,
                                      Matrix &branchNodeIncidenceMatrix) {

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

 double EMT::Ph1::IdealTransformer::getNumberOfBranches() {
   return 2;
 }