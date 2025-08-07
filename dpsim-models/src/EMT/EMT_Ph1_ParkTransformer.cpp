
#include <dpsim-models/EMT/EMT_Ph1_InertiaMoment.h>
#include <dpsim-models/EMT/EMT_Ph1_ParkTransformer.h>

using namespace CPS;

EMT::Ph1::ParkTransformer::ParkTransformer(String uid, String name,
                                           Logger::Level logLevel)
    : MNASimPowerComp<Real>(uid, name, true, true, logLevel),
      mOmega(mAttributes->create<Real>("omega")),
      mTheta(mAttributes->create<Real>("theta")) {

  setVirtualNodeNumber(3);
  setTerminalNumber(6);

  SPDLOG_LOGGER_INFO(mSLog, "Create {} {}", this->type(), name);
  **mIntfVoltage = Matrix::Zero(3, 1);
  **mIntfCurrent = Matrix::Zero(3, 1);
}

void EMT::Ph1::ParkTransformer::setIsOmegaConstant(bool isOmegaConstant) {
  mIsOmegaConstant = isOmegaConstant;
}

void EMT::Ph1::ParkTransformer::setInitialValues(Real omega, Real theta) {
  mOldOmega = omega;
  **mOmega = omega;
  mInitialTheta = theta;
  **mTheta = theta;
}

void EMT::Ph1::ParkTransformer::mnaCompInitialize(
    Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
  updateMatrixNodeIndices();
  mTimeStep = timeStep;
}

void EMT::Ph1::ParkTransformer::mnaCompAddPreStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes) {
  // add pre-step dependencies of component itself
  prevStepDependencies.push_back(mTheta);
  attributeDependencies.push_back(mOmega);
  modifiedAttributes.push_back(mTheta);
}

void EMT::Ph1::ParkTransformer::mnaCompPreStep(Real time, Int timeStepCount) {
  if (mIsOmegaConstant) {
    **mTheta = std::fmod(mInitialTheta + **mOmega * time, 2 * M_PI);
  } else {
    **mTheta = **mTheta + (mTimeStep / 2) * (**mOmega + mOldOmega);
    **mTheta = std::fmod(**mTheta, 2 * M_PI);
  }
}

void EMT::Ph1::ParkTransformer::mnaCompApplySystemMatrixStamp(
    SparseMatrixRow &systemMatrix) {

  // Park transformation equations

  // Entries related to virtual node 1
  if (terminalNotGrounded(0)) {

    Math::setMatrixElement(systemMatrix, mVirtualNodes[0]->matrixNodeIndex(),
                           matrixNodeIndex(0), 1); // virtual node 1 / node a
    Math::setMatrixElement(systemMatrix, matrixNodeIndex(0),
                           mVirtualNodes[0]->matrixNodeIndex(),
                           1); // node a / virtual node 1

    Math::setMatrixElement(
        systemMatrix, mVirtualNodes[0]->matrixNodeIndex(), matrixNodeIndex(3),
        -sqrt(2.0 / 3.0) * cos(**mTheta)); //  virtual node 1 / node d
    Math::setMatrixElement(
        systemMatrix, matrixNodeIndex(3), mVirtualNodes[0]->matrixNodeIndex(),
        -sqrt(2.0 / 3.0) * cos(**mTheta)); // node d / virtual node 1

    Math::setMatrixElement(
        systemMatrix, mVirtualNodes[0]->matrixNodeIndex(), matrixNodeIndex(4),
        +sqrt(2.0 / 3.0) * sin(**mTheta)); // virtual node 1 / node q
    Math::setMatrixElement(
        systemMatrix, matrixNodeIndex(4), mVirtualNodes[0]->matrixNodeIndex(),
        +sqrt(2.0 / 3.0) * sin(**mTheta)); // node q / virtual node 1

    Math::setMatrixElement(systemMatrix, mVirtualNodes[0]->matrixNodeIndex(),
                           matrixNodeIndex(5),
                           -sqrt(1.0 / 3.0)); // virtual node 1 / node 0
    Math::setMatrixElement(systemMatrix, matrixNodeIndex(5),
                           mVirtualNodes[0]->matrixNodeIndex(),
                           -sqrt(1.0 / 3.0)); // node 0 / virtual node 1
  }

  // Entries related to virtual node 2
  if (terminalNotGrounded(1)) {
    Math::setMatrixElement(systemMatrix, mVirtualNodes[1]->matrixNodeIndex(),
                           matrixNodeIndex(1), 1); // node b / virtual node 2
    Math::setMatrixElement(systemMatrix, matrixNodeIndex(1),
                           mVirtualNodes[1]->matrixNodeIndex(),
                           1); // virtual node 2 / node b

    Math::setMatrixElement(
        systemMatrix, mVirtualNodes[1]->matrixNodeIndex(), matrixNodeIndex(3),
        -sqrt(2.0 / 3.0) *
            cos(**mTheta - 2 * M_PI / 3)); // virtual node 2 / node d
    Math::setMatrixElement(
        systemMatrix, matrixNodeIndex(3), mVirtualNodes[1]->matrixNodeIndex(),
        -sqrt(2.0 / 3.0) *
            cos(**mTheta - 2 * M_PI / 3)); // node d / virtual node 2

    Math::setMatrixElement(
        systemMatrix, mVirtualNodes[1]->matrixNodeIndex(), matrixNodeIndex(4),
        +sqrt(2.0 / 3.0) *
            sin(**mTheta - 2 * M_PI / 3)); // virtual node 2 / node q
    Math::setMatrixElement(
        systemMatrix, matrixNodeIndex(4), mVirtualNodes[1]->matrixNodeIndex(),
        +sqrt(2.0 / 3.0) *
            sin(**mTheta - 2 * M_PI / 3)); // node q / virtual node 2

    Math::setMatrixElement(systemMatrix, mVirtualNodes[1]->matrixNodeIndex(),
                           matrixNodeIndex(5),
                           -sqrt(1.0 / 3.0)); // virtual node 2 / node 0
    Math::setMatrixElement(systemMatrix, matrixNodeIndex(5),
                           mVirtualNodes[1]->matrixNodeIndex(),
                           -sqrt(1.0 / 3.0)); // node 0 / virtual node 2
  }

  // Entries related to virtual node 3
  if (terminalNotGrounded(2)) {
    Math::setMatrixElement(systemMatrix, mVirtualNodes[2]->matrixNodeIndex(),
                           matrixNodeIndex(2), 1); // virtual node 3 / node c
    Math::setMatrixElement(systemMatrix, matrixNodeIndex(2),
                           mVirtualNodes[2]->matrixNodeIndex(),
                           1); // node c / virtual node 3

    Math::setMatrixElement(
        systemMatrix, mVirtualNodes[2]->matrixNodeIndex(), matrixNodeIndex(3),
        -sqrt(2.0 / 3.0) *
            cos(**mTheta + 2 * M_PI / 3)); // virtual node 3 / node d
    Math::setMatrixElement(
        systemMatrix, matrixNodeIndex(3), mVirtualNodes[2]->matrixNodeIndex(),
        -sqrt(2.0 / 3.0) *
            cos(**mTheta + 2 * M_PI / 3)); // node d / virtual node 3

    Math::setMatrixElement(
        systemMatrix, mVirtualNodes[2]->matrixNodeIndex(), matrixNodeIndex(4),
        +sqrt(2.0 / 3.0) *
            sin(**mTheta + 2 * M_PI / 3)); // virtual node 3 / node q
    Math::setMatrixElement(
        systemMatrix, matrixNodeIndex(4), mVirtualNodes[2]->matrixNodeIndex(),
        +sqrt(2.0 / 3.0) *
            sin(**mTheta + 2 * M_PI / 3)); // node q / virtual node 3

    Math::setMatrixElement(systemMatrix, mVirtualNodes[2]->matrixNodeIndex(),
                           matrixNodeIndex(5),
                           -sqrt(1.0 / 3.0)); // virtual node 3 / node 0
    Math::setMatrixElement(systemMatrix, matrixNodeIndex(5),
                           mVirtualNodes[2]->matrixNodeIndex(),
                           -sqrt(1.0 / 3.0)); // node 0 / virtual node 3
  }
}

void EMT::Ph1::ParkTransformer::mnaCompAddPostStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes,
    Attribute<Matrix>::Ptr &leftVector) {
  attributeDependencies.push_back(leftVector);
  modifiedAttributes.push_back(mOmega);
}

void EMT::Ph1::ParkTransformer::mnaCompPostStep(
    Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector) {
  if (mIsOmegaConstant == false) {
    mOldOmega = **mOmega;
    **mOmega = mOmegaReferenceNode->voltage()(0, 0);

  }
}

void EMT::Ph1::ParkTransformer::stampBranchNodeIncidenceMatrix(
    UInt branchIdx, Matrix &branchNodeIncidenceMatrix) {}