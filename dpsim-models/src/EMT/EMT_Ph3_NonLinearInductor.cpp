/* Copyright 2017-2020 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim-models/EMT/EMT_Ph3_NonLinearInductor.h>

using namespace CPS;

EMT::Ph3::NonLinearInductor::NonLinearInductor(String uid, String name,
                                       Logger::Level logLevel)
    : MNASimPowerComp<Real>(uid, name, true, true, logLevel) {
  mPhaseType = PhaseType::ABC;
  setVirtualNodeNumber(0);
  setTerminalNumber(2);
  **mIntfVoltage = Matrix::Zero(3, 1);
  **mIntfCurrent = Matrix::Zero(3, 1);
}

void EMT::Ph3::NonLinearInductor::initializeFromNodesAndTerminals(Real frequency) {

}


void EMT::Ph3::NonLinearInductor::mnaCompInitialize(
    Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
  updateMatrixNodeIndices();
}

void EMT::Ph3::NonLinearInductor::mnaCompApplyRightSideVectorStamp(
    Matrix &rightVector) {
  if (terminalNotGrounded(0)) {
    Math::setVectorElement(rightVector, matrixNodeIndex(1, 0),
                           +(**mIntfCurrent)(0, 0));
    Math::setVectorElement(rightVector, matrixNodeIndex(1, 1),
                           +(**mIntfCurrent)(1, 0));
    Math::setVectorElement(rightVector, matrixNodeIndex(1, 2),
                           +(**mIntfCurrent)(2, 0));
  }
  if (terminalNotGrounded(1)) {
    Math::setVectorElement(rightVector, matrixNodeIndex(0, 0),
                           -(**mIntfCurrent)(0, 0));
    Math::setVectorElement(rightVector, matrixNodeIndex(0, 1),
                           -(**mIntfCurrent)(1, 0));
    Math::setVectorElement(rightVector, matrixNodeIndex(0, 2),
                           -(**mIntfCurrent)(2, 0));
  }
}

void EMT::Ph3::NonLinearInductor::updateCurrent(Real time) {
 
  for(int i=0; i<3; i++) {
    mFlux(i,0) = mFlux(i,0) + 0.5 * (mOldVoltages(i,0) + mNewVoltages(i,0)) * mTimeStep;
    (**mIntfCurrent)(i, 0) = mPieceWiseCharacteristic->getCurrent(mFlux(i,0));
  }
  SPDLOG_LOGGER_DEBUG(mSLog, "\nUpdate current: {:s}",
                      Logger::matrixToString(**mIntfCurrent));
}

void EMT::Ph3::NonLinearInductor::mnaCompAddPreStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes) {
  modifiedAttributes.push_back(mRightVector);
  modifiedAttributes.push_back(mIntfVoltage);
}

void EMT::Ph3::NonLinearInductor::mnaCompPreStep(Real time, Int timeStepCount) {
  updateCurrent(time);
  mnaCompApplyRightSideVectorStamp(**mRightVector);

}

void EMT::Ph3::NonLinearInductor::mnaCompAddPostStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes,
    Attribute<Matrix>::Ptr &leftVector) {
  attributeDependencies.push_back(leftVector);
  modifiedAttributes.push_back(mIntfVoltage);
};

void EMT::Ph3::NonLinearInductor::mnaCompPostStep(
    Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector) {
  mnaCompUpdateVoltage(**leftVector);


  for(int i=0; i<3; i++) {
    mOldVoltages(i,0) = mNewVoltages(i,0);
    mNewVoltages(i,0) = (**mIntfVoltage)(i,0);
  }

}

void EMT::Ph3::NonLinearInductor::mnaCompUpdateVoltage(const Matrix &leftVector) {
  // v1 - v0
  **mIntfVoltage = Matrix::Zero(3, 1);
  if (terminalNotGrounded(1)) {
    (**mIntfVoltage)(0, 0) =
        Math::realFromVectorElement(leftVector, matrixNodeIndex(1, 0));
    (**mIntfVoltage)(1, 0) =
        Math::realFromVectorElement(leftVector, matrixNodeIndex(1, 1));
    (**mIntfVoltage)(2, 0) =
        Math::realFromVectorElement(leftVector, matrixNodeIndex(1, 2));
  }
  if (terminalNotGrounded(0)) {
    (**mIntfVoltage)(0, 0) =
        (**mIntfVoltage)(0, 0) -
        Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 0));
    (**mIntfVoltage)(1, 0) =
        (**mIntfVoltage)(1, 0) -
        Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 1));
    (**mIntfVoltage)(2, 0) =
        (**mIntfVoltage)(2, 0) -
        Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 2));
  }
}


void EMT::Ph3::NonLinearInductor::setTimeStep(Real timeStep) {
  mTimeStep = timeStep;
}