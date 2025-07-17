/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim-models/EMT/EMT_Ph1_CrossCoupling.h>
#include <dpsim-models/EMT/EMT_Ph1_Inductor.h>
#include <dpsim-models/EMT/EMT_Ph1_InertiaMoment.h>

using namespace CPS;

EMT::Ph1::CrossCoupling::CrossCoupling(String uid, String name,
                                       Logger::Level logLevel)
    : MNASimPowerComp<Real>(uid, name, true, true, logLevel),
      mAngularSpeedRef(mAttributes->create<Real>("w_src")) {
  setVirtualNodeNumber(1);
  setTerminalNumber(2);
  **mIntfVoltage = Matrix::Zero(1, 1);
  **mIntfCurrent = Matrix::Zero(1, 1);
}

void EMT::Ph1::CrossCoupling::setParameters(Real omega) {
  **mAngularSpeedRef = omega;
  mConstantSpeed = 1;  // Store the constant speed for later use

  mParametersSet = true;
}

SimPowerComp<Real>::Ptr EMT::Ph1::CrossCoupling::clone(String name) {
  auto copy = CrossCoupling::make(name, mLogLevel);
  copy->setParameters(**mAngularSpeedRef);
  return copy;
}

void EMT::Ph1::CrossCoupling::mnaCompInitialize(
    Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
  updateMatrixNodeIndices();
  if (mConstantSpeed == 0) {
  (**mIntfVoltage)(0, 0) =
     (**(mInertiaMoment->mIntfVoltage))(0, 0) * mFluxLinkage;

     }

  else {
    (**mIntfVoltage)(0, 0) = mConstantSpeed * mFluxLinkage;
  }
  mTimeStep = timeStep;
}

void EMT::Ph1::CrossCoupling::mnaCompApplySystemMatrixStamp(
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

void EMT::Ph1::CrossCoupling::mnaCompApplyRightSideVectorStamp(
    Matrix &rightVector) {
  Math::setVectorElement(rightVector, mVirtualNodes[0]->matrixNodeIndex(),
                         (**mIntfVoltage)(0, 0));
}

void EMT::Ph1::CrossCoupling::updateVoltage(Real time) {

}

void EMT::Ph1::CrossCoupling::mnaCompAddPreStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes) {
  modifiedAttributes.push_back(mRightVector);
  modifiedAttributes.push_back(mIntfVoltage);
}

void EMT::Ph1::CrossCoupling::mnaCompPreStep(Real time, Int timeStepCount) {
  //updateVoltage(time);
  mnaCompApplyRightSideVectorStamp(**mRightVector);


  

    if (mNegativeSpeedTermVoltageFlag == true) {
            mFluxLinkage = mFluxLinkage - ( mTimeStep / 2 ) * (mOldVoltage + mVoltage);
    } else {
      // If the voltage is not negative, we can use the voltage as it is
            mFluxLinkage = mFluxLinkage + ( mTimeStep / 2 ) * (mOldVoltage + mVoltage);

    }
}

void EMT::Ph1::CrossCoupling::mnaCompAddPostStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes,
    Attribute<Matrix>::Ptr &leftVector) {
  attributeDependencies.push_back(leftVector);
  modifiedAttributes.push_back(mIntfCurrent);
}

void EMT::Ph1::CrossCoupling::mnaCompPostStep(
    Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector) {
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

void EMT::Ph1::CrossCoupling::mnaCompUpdateCurrent(const Matrix &leftVector) {
  (**mIntfCurrent)(0, 0) = Math::realFromVectorElement(
      leftVector, mVirtualNodes[0]->matrixNodeIndex());
}

void EMT::Ph1::CrossCoupling::stampBranchNodeIncidenceMatrix(
    UInt branchIdx, Matrix &branchNodeIncidenceMatrix) {
  if (terminalNotGrounded(0)) {
    branchNodeIncidenceMatrix(branchIdx, matrixNodeIndex(0)) = 1.0;
  }
  if (terminalNotGrounded(1)) {
    branchNodeIncidenceMatrix(branchIdx, matrixNodeIndex(1)) = -1.0;
  }
}
