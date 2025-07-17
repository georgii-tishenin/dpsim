/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim-models/EMT/EMT_Ph1_BackEMF.h>
#include <dpsim-models/EMT/EMT_Ph1_InertiaMoment.h>

using namespace CPS;

EMT::Ph1::BackEMF::BackEMF(String uid, String name,
                                       Logger::Level logLevel)
    : MNASimPowerComp<Real>(uid, name, true, true, logLevel),
      mFieldFluxLinkageRef(mAttributes->create<Real>("L_ref")) {
  setVirtualNodeNumber(1);
  setTerminalNumber(2);
  **mIntfVoltage = Matrix::Zero(1, 1);
  **mIntfCurrent = Matrix::Zero(1, 1);
}

void EMT::Ph1::BackEMF::setParameters(Real FieldFluxLinkage) {
  **mFieldFluxLinkageRef = FieldFluxLinkage;;

  mParametersSet = true;
}

SimPowerComp<Real>::Ptr EMT::Ph1::BackEMF::clone(String name) {
  auto copy = BackEMF::make(name, mLogLevel);
  copy->setParameters(**mFieldFluxLinkageRef);
  return copy;
}

void EMT::Ph1::BackEMF::mnaCompInitialize(
    Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
  updateMatrixNodeIndices();
  (**mIntfVoltage)(0, 0) =
      **mFieldFluxLinkageRef * (**(mInertiaMoment->mIntfVoltage))(0, 0);

  mTimeStep = timeStep;
}

void EMT::Ph1::BackEMF::mnaCompApplySystemMatrixStamp(
    SparseMatrixRow &systemMatrix) {
  if (terminalNotGrounded(0)) {
    Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0),
                             mVirtualNodes[0]->matrixNodeIndex(), -1);
    Math::addToMatrixElement(systemMatrix, mVirtualNodes[0]->matrixNodeIndex(),
                             matrixNodeIndex(0), -1);
  }
  if (terminalNotGrounded(1)) {
    Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1),
                             mVirtualNodes[0]->matrixNodeIndex(), 1);
    Math::addToMatrixElement(systemMatrix, mVirtualNodes[0]->matrixNodeIndex(),
                             matrixNodeIndex(1), 1);
  }

  if (terminalNotGrounded(0)) {
    SPDLOG_LOGGER_INFO(mSLog, "Add {:f} to system at ({:d},{:d})", -1.,
                       matrixNodeIndex(0), mVirtualNodes[0]->matrixNodeIndex());
    SPDLOG_LOGGER_INFO(mSLog, "Add {:f} to system at ({:d},{:d})", -1.,
                       mVirtualNodes[0]->matrixNodeIndex(), matrixNodeIndex(0));
  }
  if (terminalNotGrounded(1)) {
    SPDLOG_LOGGER_INFO(mSLog, "Add {:f} to system at ({:d},{:d})", 1.,
                       matrixNodeIndex(1), mVirtualNodes[0]->matrixNodeIndex());
    SPDLOG_LOGGER_INFO(mSLog, "Add {:f} to system at ({:d},{:d})", 1.,
                       mVirtualNodes[0]->matrixNodeIndex(), matrixNodeIndex(1));
  }
}

void EMT::Ph1::BackEMF::mnaCompApplyRightSideVectorStamp(
    Matrix &rightVector) {
  Math::setVectorElement(rightVector, mVirtualNodes[0]->matrixNodeIndex(),
                         (**mIntfVoltage)(0, 0));
}

void EMT::Ph1::BackEMF::updateVoltage(Real time) {
  Real FieldFlux = mFieldFluxLinkageRef->get();  //Field flux linkage
  Real omega = (**(mInertiaMoment->mIntfVoltage))(0, 0);  // Angular speed of the machine is derived from the inertia moment
    (**mIntfVoltage)(0, 0) =
        omega * FieldFlux;

}

void EMT::Ph1::BackEMF::mnaCompAddPreStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes) {
  attributeDependencies.push_back(mFieldFluxLinkageRef);
  modifiedAttributes.push_back(mRightVector);
  modifiedAttributes.push_back(mIntfVoltage);
}

void EMT::Ph1::BackEMF::mnaCompPreStep(Real time, Int timeStepCount) {
  updateVoltage(time);
  mnaCompApplyRightSideVectorStamp(**mRightVector);
}

void EMT::Ph1::BackEMF::mnaCompAddPostStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes,
    Attribute<Matrix>::Ptr &leftVector) {
  attributeDependencies.push_back(leftVector);
  modifiedAttributes.push_back(mIntfCurrent);
}

void EMT::Ph1::BackEMF::mnaCompPostStep(
    Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector) {
  mnaCompUpdateCurrent(**leftVector);
}

void EMT::Ph1::BackEMF::mnaCompUpdateCurrent(const Matrix &leftVector) {
  (**mIntfCurrent)(0, 0) = Math::realFromVectorElement(
      leftVector, mVirtualNodes[0]->matrixNodeIndex());
}

void EMT::Ph1::BackEMF::stampBranchNodeIncidenceMatrix(
    UInt branchIdx, Matrix &branchNodeIncidenceMatrix) {
  if (terminalNotGrounded(0)) {
    branchNodeIncidenceMatrix(branchIdx, matrixNodeIndex(0)) = 1.0;
  }
  if (terminalNotGrounded(1)) {
    branchNodeIncidenceMatrix(branchIdx, matrixNodeIndex(1)) = -1.0;
  }
}
