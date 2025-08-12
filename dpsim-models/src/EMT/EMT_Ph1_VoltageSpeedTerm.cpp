#include <dpsim-models/EMT/EMT_Ph1_VoltageSpeedTerm.h>

using namespace CPS;

EMT::Ph1::VoltageSpeedTerm::VoltageSpeedTerm(
    String uid, String name, Logger::Level logLevel)
    : MNASimPowerComp<Real>(uid, name, true, true, logLevel),
      mOmega(mAttributes->create<Real>("omega")){

  setVirtualNodeNumber(2);
  setTerminalNumber(4);

  SPDLOG_LOGGER_INFO(mSLog, "Create {} {}", this->type(), name);
  **mIntfVoltage = Matrix::Zero(1, 1);
  **mIntfCurrent = Matrix::Zero(1, 1);
}

void EMT::Ph1::VoltageSpeedTerm::setInitialOmega(Real omega) {
  **mOmega = omega;
}

void EMT::Ph1::VoltageSpeedTerm::setInductance(Real inductance) {
  mInductance = inductance;
}

void EMT::Ph1::VoltageSpeedTerm::mnaCompInitialize(
    Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
  updateMatrixNodeIndices();
}

void EMT::Ph1::VoltageSpeedTerm::mnaCompApplySystemMatrixStamp(
    SparseMatrixRow &systemMatrix) {

  int m = matrixNodeIndex(0);
  int n = matrixNodeIndex(1);
  int p = matrixNodeIndex(2);
  int q = matrixNodeIndex(3);

  int imn = mVirtualNodes[0]->matrixNodeIndex();
  int ipq = mVirtualNodes[1]->matrixNodeIndex();

  // Ideal transformer equations
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

  double alpha = **mOmega * mInductance;
  if (mIsNegative) {
    alpha = -alpha;
  }

  Math::setMatrixElement(systemMatrix, ipq, imn, -alpha);
}

void EMT::Ph1::VoltageSpeedTerm::stampBranchNodeIncidenceMatrix(
    UInt branchIdx, Matrix &branchNodeIncidenceMatrix) {}