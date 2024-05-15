/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim-models/DP/DP_Ph1_VSIVoltageControlDQ.h>

using namespace CPS;

DP::Ph1::VSIVoltageControlDQ::VSIVoltageControlDQ(String uid, String name,
                                                  Logger::Level logLevel,
                                                  Bool modelAsCurrentSource,
                                                  Bool withInterfaceResistor)
    : CompositePowerComp<Complex>(uid, name, true, true, logLevel),
      VSIVoltageSourceInverterDQ<Complex>(this->mSLog, mAttributes,
                                          modelAsCurrentSource,
                                          withInterfaceResistor) {

  setTerminalNumber(1);
  setVirtualNodeNumber(this->determineNumberOfVirtualNodes());

  **mIntfVoltage = MatrixComp::Zero(1, 1);
  **mIntfCurrent = MatrixComp::Zero(1, 1);
}

void DP::Ph1::VSIVoltageControlDQ::createSubComponents() {
  // voltage source
  if (!mModelAsCurrentSource) {
    mSubCtrledVoltageSource =
        DP::Ph1::VoltageSource::make(**mName + "_src", mLogLevel);
    addMNASubComponent(mSubCtrledVoltageSource, MNA_SUBCOMP_TASK_ORDER::NO_TASK,
                       MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);
  }

  // RL Element as part of the LC filter
  mSubFilterRL = DP::Ph1::ResIndSeries::make(**mName + "_FilterRL", mLogLevel);
  addMNASubComponent(mSubFilterRL, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT,
                     MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);
  mSubFilterRL->setParameters(mRf, mLf);

  // Capacitor as part of the LC filter
  mSubCapacitorF = DP::Ph1::Capacitor::make(**mName + "_CapF", mLogLevel);
  addMNASubComponent(mSubCapacitorF, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT,
                     MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);
  mSubCapacitorF->setParameters(mCf);

  // optinal: interface resistor
  if (mWithInterfaceResistor) {
    mSubResistorC = DP::Ph1::Resistor::make(**mName + "_ResC", mLogLevel);
    mSubResistorC->setParameters(mRc);
    addMNASubComponent(mSubResistorC,
                       MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT,
                       MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, false);
  }
}

void DP::Ph1::VSIVoltageControlDQ::connectSubComponents() {
  // TODO: COULD WE MOVE THIS FUNCTION TO THE BASE CLASS?
  if (!mModelAsCurrentSource)
    mSubCtrledVoltageSource->connect({SimNode::GND, mVirtualNodes[0]});
  if (mWithInterfaceResistor) {
    // with interface resistor
    mSubFilterRL->connect({mVirtualNodes[1], mVirtualNodes[0]});
    mSubCapacitorF->connect({SimNode::GND, mVirtualNodes[1]});
    mSubResistorC->connect({mTerminals[0]->node(), mVirtualNodes[1]});
  } else {
    // without interface resistor
    mSubFilterRL->connect({mTerminals[0]->node(), mVirtualNodes[0]});
    mSubCapacitorF->connect({SimNode::GND, mTerminals[0]->node()});
  }
}

void DP::Ph1::VSIVoltageControlDQ::initializeFromNodesAndTerminals(
    Real frequency) {
  // terminal powers in consumer system -> convert to generator system
  **mPower = -terminal(0)->singlePower();

  // set initial interface quantities --> Current flowing into the inverter is positive
  (**mIntfVoltage)(0, 0) = initialSingleVoltage(0);
  (**mIntfCurrent)(0, 0) = std::conj(**mPower / (**mIntfVoltage)(0, 0));

  // initialize filter variables and set initial voltage of virtual nodes
  initializeFilterVariables((**mIntfVoltage)(0, 0), (**mIntfCurrent)(0, 0),
                            mVirtualNodes);

  // calculate initial source value
  (**mSourceValue)(0, 0) =
      Math::rotatingFrame2to1(**mSourceValue_dq, **mThetaSys, **mThetaInv);

  // Create & Initialize electrical subcomponents
  this->connectSubComponents();
  for (auto subcomp : mSubComponents) {
    subcomp->initialize(mFrequencies);
    subcomp->initializeFromNodesAndTerminals(frequency);
  }

  // TODO: droop
  **mOmega = mOmegaNom;

  SPDLOG_LOGGER_INFO(mSLog,
                     "\n--- Initialization from powerflow ---"
                     "\nTerminal 0 connected to {} = sim node {}"
                     "\nInverter terminal voltage: {}[V]"
                     "\nInverter output current: {}[A]",
                     mTerminals[0]->node()->name(),
                     mTerminals[0]->node()->matrixNodeIndex(),
                     Logger::phasorToString((**mIntfVoltage)(0, 0)),
                     Logger::phasorToString((**mIntfCurrent)(0, 0)));
  mSLog->flush();
}

void DP::Ph1::VSIVoltageControlDQ::mnaParentInitialize(
    Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
  this->updateMatrixNodeIndices();
  mTimeStep = timeStep;
  if (mWithControl)
    mVSIController->initialize(**mSourceValue_dq, **mVcap_dq, **mIfilter_dq,
                               mTimeStep, mModelAsCurrentSource);
}

void DP::Ph1::VSIVoltageControlDQ::mnaParentAddPreStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes) {
  modifiedAttributes.push_back(mRightVector);
  prevStepDependencies.push_back(mIntfVoltage);
}

void DP::Ph1::VSIVoltageControlDQ::mnaParentPreStep(Real time,
                                                    Int timeStepCount) {
  // get measurements
  **mVcap_dq = Math::rotatingFrame2to1((**mSubCapacitorF->mIntfVoltage)(0, 0),
                                       **mThetaInv, **mThetaSys);
  **mIfilter_dq = Math::rotatingFrame2to1((**mSubFilterRL->mIntfCurrent)(0, 0),
                                          **mThetaInv, **mThetaSys);

  // TODO: droop
  //if (mWithDroop)
  //	mDroop->signalStep(time, timeStepCount);

  //  VCO Step
  **mThetaInv = **mThetaInv + mTimeStep * **mOmega;

  // Update nominal system angle
  **mThetaSys = **mThetaSys + mTimeStep * mOmegaNom;

  //
  if (mWithControl)
    **mSourceValue_dq = mVSIController->step(**mVcap_dq, **mIfilter_dq);

  // Transformation interface backward
  (**mSourceValue)(0, 0) =
      Math::rotatingFrame2to1(**mSourceValue_dq, **mThetaSys, **mThetaInv);

  // set reference voltage of voltage source
  if (!mModelAsCurrentSource) {
    // pre-step of voltage source
    **mSubCtrledVoltageSource->mVoltageRef = (**mSourceValue)(0, 0);
    std::dynamic_pointer_cast<MNAInterface>(mSubCtrledVoltageSource)
        ->mnaPreStep(time, timeStepCount);
  }

  // stamp right side vector
  mnaApplyRightSideVectorStamp(**mRightVector);
}

void DP::Ph1::VSIVoltageControlDQ::mnaParentApplyRightSideVectorStamp(
    Matrix &rightVector) {
  if (mModelAsCurrentSource)
    Math::addToVectorElement(**mRightVector,
                             mVirtualNodes[0]->matrixNodeIndex(),
                             (**mSourceValue)(0, 0));
}

void DP::Ph1::VSIVoltageControlDQ::mnaParentAddPostStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes,
    Attribute<Matrix>::Ptr &leftVector) {
  attributeDependencies.push_back(leftVector);
  modifiedAttributes.push_back(mIntfVoltage);
  modifiedAttributes.push_back(mIntfCurrent);
}

void DP::Ph1::VSIVoltageControlDQ::mnaParentPostStep(
    Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector) {
  mnaCompUpdateCurrent(**leftVector);
  mnaCompUpdateVoltage(**leftVector);
  updatePower();
}

void DP::Ph1::VSIVoltageControlDQ::mnaCompUpdateCurrent(
    const Matrix &leftvector) {
  **mIntfCurrent =
      mSubCapacitorF->mIntfCurrent->get() + mSubFilterRL->mIntfCurrent->get();
}

void DP::Ph1::VSIVoltageControlDQ::mnaCompUpdateVoltage(
    const Matrix &leftVector) {
  (**mIntfVoltage)(0, 0) =
      Math::complexFromVectorElement(leftVector, matrixNodeIndex(0));
}

void DP::Ph1::VSIVoltageControlDQ::updatePower() {
  **mPower = (**mIntfVoltage)(0, 0) * std::conj((**mIntfCurrent)(0, 0));
}
