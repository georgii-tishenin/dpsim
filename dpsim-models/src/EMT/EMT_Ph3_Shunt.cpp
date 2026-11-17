/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim-models/EMT/EMT_Ph3_SSN_Primitives.h>
#include <dpsim-models/EMT/EMT_Ph3_Shunt.h>
#include <dpsim-models/MNAStampUtils.h>

#include <stdexcept>

using namespace CPS;

EMT::Ph3::Shunt::Shunt(String uid, String name, Logger::Level logLevel)
    : CompositePowerComp<Real>(uid, name, true, true, logLevel) {
  mOutputUpdates = mAttributes->create<UInt>("output_updates");
  **mOutputUpdates = 0;
  mTerminalCurrent0 = mAttributes->create<Matrix>("i_terminal_0");
  **mTerminalCurrent0 = Matrix::Zero(3, 1);
  mPhaseType = PhaseType::ABC;
  mConductance = Matrix::Zero(3, 3);
  mSusceptance = Matrix::Zero(3, 3);
  SPDLOG_LOGGER_INFO(mSLog, "Create {} of type {}", this->type(), name);
  setTerminalNumber(1);
}

EMT::Ph3::Shunt::Shunt(String uid, String name, Bool useSubcomponents,
                       Logger::Level logLevel)
    : Shunt(uid, name, logLevel) {
  mUseSubcomponents = useSubcomponents;
}

void EMT::Ph3::Shunt::setParameters(Real conductance, Real susceptance) {
  this->setParameters(Math::singlePhaseParameterToThreePhase(conductance),
                      Math::singlePhaseParameterToThreePhase(susceptance));
}

void EMT::Ph3::Shunt::setParameters(Matrix conductance, Matrix susceptance) {
  mConductance = conductance;
  mSusceptance = susceptance;
  SPDLOG_LOGGER_INFO(mSLog, "Conductance={} [S] Susceptance={} [Ohm] ",
                     Logger::matrixToString(mConductance),
                     Logger::matrixToString(mSusceptance));
  mParametersSet = true;
}

/// MNA Section
void EMT::Ph3::Shunt::initializeParentFromNodesAndTerminals(Real frequency) {

  // If no dynamic-domain admittance was configured, derive it from the
  // power-flow operating point copied to this terminal by
  // SystemTopology::initWithPowerflow().  Terminal power uses the network
  // injection sign here: a positive conductance consumes active power.
  if (!mParametersSet) {
    const Complex voltage = mTerminals[0]->initialSingleVoltage();
    if (std::norm(voltage) <= 0.0)
      throw std::invalid_argument(
          "Shunt requires a non-zero terminal voltage when initialized from "
          "power-flow data.");

    const Complex terminalPower = mTerminals[0]->singlePower();
    const Real conductance = -terminalPower.real() / std::norm(voltage);
    const Real susceptance = terminalPower.imag() / std::norm(voltage);
    setParameters(conductance, susceptance);
  }

  // Static calculation
  Real omega = 2. * PI * frequency;
  MatrixComp admittance = Matrix::Zero(3, 3);
  admittance << Complex(mConductance(0, 0), mSusceptance(0, 0)), 0, 0, 0,
      Complex(mConductance(1, 1), mSusceptance(1, 1)), 0, 0, 0,
      Complex(mConductance(2, 2), mSusceptance(2, 2));
  MatrixComp vInitABC = MatrixComp::Zero(3, 1);
  vInitABC(0, 0) = RMS3PH_TO_PEAK1PH * mTerminals[0]->initialSingleVoltage();
  vInitABC(1, 0) = vInitABC(0, 0) * SHIFT_TO_PHASE_B;
  vInitABC(2, 0) = vInitABC(0, 0) * SHIFT_TO_PHASE_C;
  **mIntfVoltage = vInitABC.real();
  // The shunt interface current follows the component's terminal-current
  // convention (positive into the component), opposite to the current from
  // the network node into the shunt admittance.
  **mIntfCurrent = -(admittance * vInitABC).real();
  **mTerminalCurrent0 = -**mIntfCurrent;

  if (!mUseSubcomponents) {
    MatrixComp susceptiveAdmittance =
        MatrixComp(Complex(0, 1) * mSusceptance) * vInitABC;
    if (mSusceptance(0, 0) > 0)
      mCapBranchCurrent = susceptiveAdmittance.real();
    if (mSusceptance(0, 0) < 0)
      mIndBranchCurrent = susceptiveAdmittance.real();
    return;
  }

  // Create series rl sub component
  if (mConductance(0, 0) > 0) {
    mSubResistor = makePassivePrimitive<EMT::Ph3::Resistor>(
        mUseElementarySSN, **mName + "_Res", mLogLevel);
    mSubResistor->connect({mTerminals[0]->node(), SimNode::GND});
    mSubResistor->setParameters(mConductance.inverse());
    mSubResistor->initialize(mFrequencies);
    mSubResistor->initializeFromNodesAndTerminals(frequency);
    addMNASubComponent(mSubResistor, MNA_SUBCOMP_TASK_ORDER::NO_TASK,
                       MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, false);
  }

  if (mSusceptance(0, 0) > 0) {
    mSubCapacitor = makePassivePrimitive<EMT::Ph3::Capacitor>(
        mUseElementarySSN, **mName + "_cap", mLogLevel);
    mSubCapacitor->setParameters(mSusceptance * (1 / omega));
    mSubCapacitor->connect({mTerminals[0]->node(), SimNode::GND});
    mSubCapacitor->initialize(mFrequencies);
    mSubCapacitor->initializeFromNodesAndTerminals(frequency);
    addMNASubComponent(mSubCapacitor,
                       MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT,
                       MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);
  }

  if (mSusceptance(0, 0) < 0) {
    mSubInductor = makePassivePrimitive<EMT::Ph3::Inductor>(
        mUseElementarySSN, **mName + "_ind", mLogLevel);
    mSubInductor->setParameters(-mSusceptance.inverse() * (1 / omega));
    mSubInductor->connect({mTerminals[0]->node(), SimNode::GND});
    mSubInductor->initialize(mFrequencies);
    mSubInductor->initializeFromNodesAndTerminals(frequency);
    addMNASubComponent(mSubInductor, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT,
                       MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);
  }

  SPDLOG_LOGGER_INFO(mSLog,
                     "\n--- Initialization from powerflow ---"
                     "\nVoltage across: {:s}"
                     "\nCurrent: {:s}"
                     "\nTerminal voltage: {:s}"
                     "\n--- Initialization from powerflow finished ---",
                     Logger::phasorToString((**mIntfVoltage)(0, 0)),
                     Logger::phasorToString((**mIntfCurrent)(0, 0)),
                     Logger::phasorToString(initialSingleVoltage(0)));
}

void EMT::Ph3::Shunt::mnaParentInitialize(Real omega, Real timeStep,
                                          Attribute<Matrix>::Ptr leftVector) {
  if (mUseSubcomponents)
    return;

  if (mSusceptance(0, 0) > 0)
    mCapEquivalentConductance = (2.0 / (omega * timeStep)) * mSusceptance;
  if (mSusceptance(0, 0) < 0)
    mIndEquivalentConductance = (-0.5 * omega * timeStep) * mSusceptance;
  mCapEquivalentCurrent =
      -mCapBranchCurrent - mCapEquivalentConductance * **mIntfVoltage;
  mIndEquivalentCurrent =
      mIndEquivalentConductance * **mIntfVoltage + mIndBranchCurrent;
}

void EMT::Ph3::Shunt::mnaParentApplySystemMatrixStamp(
    SparseMatrixRow &systemMatrix) {
  if (mUseSubcomponents)
    return;
  const Matrix conductance =
      mConductance + mCapEquivalentConductance + mIndEquivalentConductance;
  MNAStampUtils::stampConductanceMatrix(conductance, systemMatrix,
                                        matrixNodeIndex(0), 0,
                                        terminalNotGrounded(0), false, mSLog);
}

void EMT::Ph3::Shunt::mnaParentApplyRightSideVectorStamp(Matrix &rightVector) {
  if (mUseSubcomponents || !terminalNotGrounded(0))
    return;
  for (Int phase = 0; phase < 3; ++phase) {
    Math::setVectorElement(rightVector, matrixNodeIndex(0, phase),
                           mCapEquivalentCurrent(phase, 0) +
                               mIndEquivalentCurrent(phase, 0));
  }
}

void EMT::Ph3::Shunt::mnaParentAddPreStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes) {
  modifiedAttributes.push_back(mRightVector);
}

void EMT::Ph3::Shunt::mnaParentPreStep(Real time, Int timeStepCount) {
  PassiveOperationProfile::Scope operation(mOperationProfile,
                                           PassiveOperationProfile::RhsLocal);
  if (!mUseSubcomponents) {
    mCapEquivalentCurrent =
        -mCapBranchCurrent - mCapEquivalentConductance * **mIntfVoltage;
    mIndEquivalentCurrent =
        mIndEquivalentConductance * **mIntfVoltage + mIndBranchCurrent;
  }
  mnaCompApplyRightSideVectorStamp(**mRightVector);
}

void EMT::Ph3::Shunt::mnaParentAddPostStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes,
    Attribute<Matrix>::Ptr &leftVector) {
  attributeDependencies.push_back(leftVector);
  modifiedAttributes.push_back(mIntfVoltage);
  modifiedAttributes.push_back(mIntfCurrent);
  modifiedAttributes.push_back(mOutputUpdates);
  modifiedAttributes.push_back(mTerminalCurrent0);
}

void EMT::Ph3::Shunt::mnaParentPostStep(Real time, Int timeStepCount,
                                        Attribute<Matrix>::Ptr &leftVector) {
  PassiveOperationProfile::Scope operation(
      mOperationProfile, PassiveOperationProfile::VoltageRecovery);
  ++**mOutputUpdates;
  this->mnaUpdateVoltage(**leftVector);
  operation.next(PassiveOperationProfile::StateAdvance);
  this->mnaUpdateCurrent(**leftVector);
  operation.next(PassiveOperationProfile::Output);
  **mTerminalCurrent0 = -**mIntfCurrent;
}

void EMT::Ph3::Shunt::mnaCompUpdateVoltage(const Matrix &leftVector) {
  (**mIntfVoltage)(0, 0) =
      Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 0));
  (**mIntfVoltage)(1, 0) =
      Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 1));
  (**mIntfVoltage)(2, 0) =
      Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 2));
}

void EMT::Ph3::Shunt::mnaCompUpdateCurrent(const Matrix &leftVector) {
  if (!mUseSubcomponents) {
    mCapBranchCurrent =
        mCapEquivalentConductance * **mIntfVoltage + mCapEquivalentCurrent;
    mIndBranchCurrent =
        mIndEquivalentConductance * **mIntfVoltage + mIndEquivalentCurrent;
    **mIntfCurrent = -(
        (mConductance + mCapEquivalentConductance + mIndEquivalentConductance) *
            **mIntfVoltage +
        mCapEquivalentCurrent + mIndEquivalentCurrent);
    return;
  }
  (**mIntfCurrent).setZero();

  if (mConductance(0, 0) > 0)
    **mIntfCurrent += mSubResistor->intfCurrent();
  if (mSusceptance(0, 0) > 0)
    **mIntfCurrent += mSubCapacitor->intfCurrent();
  if (mSusceptance(0, 0) < 0)
    **mIntfCurrent += mSubInductor->intfCurrent();
}
