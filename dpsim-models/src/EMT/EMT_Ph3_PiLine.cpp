/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim-models/EMT/EMT_Ph3_PiLine.h>
#include <dpsim-models/EMT/EMT_Ph3_SSN_Primitives.h>
#include <dpsim-models/MNAStampUtils.h>

using namespace CPS;

EMT::Ph3::PiLine::PiLine(String uid, String name, Logger::Level logLevel)
    : PiLine(uid, name, true, logLevel) {}

EMT::Ph3::PiLine::PiLine(String uid, String name, Bool useSubcomponents,
                         Logger::Level logLevel)
    : Base::Ph3::PiLine(mAttributes), CompositePowerComp<Real>(uid, name, true,
                                                               true, logLevel) {
  mUseSubcomponents = useSubcomponents;
  mDirectState = mAttributes->create<Matrix>("direct_state");
  **mDirectState = Matrix::Zero(18, 1);
  mOutputUpdates = mAttributes->create<UInt>("output_updates");
  **mOutputUpdates = 0;
  mTerminalCurrent0 = mAttributes->create<Matrix>("i_terminal_0");
  **mTerminalCurrent0 = Matrix::Zero(3, 1);
  mTerminalCurrent1 = mAttributes->create<Matrix>("i_terminal_1");
  **mTerminalCurrent1 = Matrix::Zero(3, 1);
  mPhaseType = PhaseType::ABC;
  // The direct path still keeps the series midpoint as an MNA unknown. This
  // preserves the exact conventional R-L companion equations while removing
  // the separate child components and their scheduler tasks.
  setVirtualNodeNumber(1);
  setTerminalNumber(2);

  SPDLOG_LOGGER_INFO(mSLog, "Create {} {}", this->type(), name);
  **mIntfVoltage = Matrix::Zero(3, 1);
  **mIntfCurrent = Matrix::Zero(3, 1);

  mSLog->flush();
}

/// DEPRECATED: Delete method
SimPowerComp<Real>::Ptr EMT::Ph3::PiLine::clone(String name) {
  auto copy = PiLine::make(name, mUseSubcomponents, mLogLevel);
  copy->setUseDefaultShuntConductance(mUseDefaultShuntConductance);
  copy->setParameters(**mSeriesRes, **mSeriesInd, **mParallelCap,
                      **mParallelCond);
  copy->setUseElementarySSN(mUseElementarySSN);
  return copy;
}

void EMT::Ph3::PiLine::createSubComponents() {
  if (mSubCompCreated)
    return;
  mSubCompCreated = true;

  mHasSeriesResistance = !(**mSeriesRes).isZero(0.);
  if (!mHasSeriesResistance)
    setVirtualNodeNumber(0);

  if (!mUseSubcomponents) {
    if (mUseDefaultShuntConductance && (**mParallelCond)(0, 0) <= 0.) {
      **mParallelCond = Matrix::Zero(3, 3);
      (**mParallelCond).diagonal().setConstant(1e-6);
    }
    return;
  }

  // By default there is always a small conductance to ground to
  // avoid problems with floating nodes.
  Matrix defaultParallelCond = Matrix::Zero(3, 3);
  defaultParallelCond << 1e-6, 0, 0, 0, 1e-6, 0, 0, 0, 1e-6;
  if (mUseDefaultShuntConductance && (**mParallelCond)(0, 0) <= 0.)
    **mParallelCond = defaultParallelCond;

  // Create series sub components
  // Do not approximate a zero resistance by a small resistance: connect the
  // ordinary inductor directly to terminal 0 and omit the unused MNA node.
  if (mHasSeriesResistance) {
    mSubSeriesResistor = makePassivePrimitive<EMT::Ph3::Resistor>(
        mUseElementarySSN, **mName + "_res", mLogLevel);
    mSubSeriesResistor->setParameters(**mSeriesRes);
    mSubSeriesResistor->connect({mTerminals[0]->node(), mVirtualNodes[0]});
    addMNASubComponent(mSubSeriesResistor,
                       MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT,
                       MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, false);
  }

  mSubSeriesInductor = makePassivePrimitive<EMT::Ph3::Inductor>(
      mUseElementarySSN, **mName + "_ind", mLogLevel);
  mSubSeriesInductor->setParameters(**mSeriesInd);
  mSubSeriesInductor->connect(
      {mHasSeriesResistance ? mVirtualNodes[0] : mTerminals[0]->node(),
       mTerminals[1]->node()});
  addMNASubComponent(mSubSeriesInductor,
                     MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT,
                     MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);

  // Create parallel sub components
  if (!(**mParallelCond).isZero(0.)) {
    mSubParallelResistor0 = makePassivePrimitive<EMT::Ph3::Resistor>(
        mUseElementarySSN, **mName + "_con0", mLogLevel);
    mSubParallelResistor0->setParameters(2. * (**mParallelCond).inverse());
    mSubParallelResistor0->connect(
        SimNode::List{SimNode::GND, mTerminals[0]->node()});
    addMNASubComponent(mSubParallelResistor0,
                       MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT,
                       MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, false);

    mSubParallelResistor1 = makePassivePrimitive<EMT::Ph3::Resistor>(
        mUseElementarySSN, **mName + "_con1", mLogLevel);
    mSubParallelResistor1->setParameters(2. * (**mParallelCond).inverse());
    mSubParallelResistor1->connect(
        SimNode::List{SimNode::GND, mTerminals[1]->node()});
    addMNASubComponent(mSubParallelResistor1,
                       MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT,
                       MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, false);
  }

  if ((**mParallelCap)(0, 0) > 0) {
    mSubParallelCapacitor0 = makePassivePrimitive<EMT::Ph3::Capacitor>(
        mUseElementarySSN, **mName + "_cap0", mLogLevel);
    mSubParallelCapacitor0->setParameters(**mParallelCap / 2.);
    mSubParallelCapacitor0->connect(
        SimNode::List{SimNode::GND, mTerminals[0]->node()});
    addMNASubComponent(mSubParallelCapacitor0,
                       MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT,
                       MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);

    mSubParallelCapacitor1 = makePassivePrimitive<EMT::Ph3::Capacitor>(
        mUseElementarySSN, **mName + "_cap1", mLogLevel);
    mSubParallelCapacitor1->setParameters(**mParallelCap / 2.);
    mSubParallelCapacitor1->connect(
        SimNode::List{SimNode::GND, mTerminals[1]->node()});
    addMNASubComponent(mSubParallelCapacitor1,
                       MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT,
                       MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);
  }
}

void EMT::Ph3::PiLine::initializeParentFromNodesAndTerminals(Real frequency) {
  // Physical port currents use the complete Pi network, including both end
  // shunts; i_intf remains the historical series-current coordinate.
  MatrixComp terminalPhasor[2];
  for (UInt terminal = 0; terminal < 2; ++terminal) {
    terminalPhasor[terminal] = MatrixComp::Zero(3, 1);
    terminalPhasor[terminal](0) =
        RMS3PH_TO_PEAK1PH * initialSingleVoltage(terminal);
    terminalPhasor[terminal](1) =
        terminalPhasor[terminal](0) * SHIFT_TO_PHASE_B;
    terminalPhasor[terminal](2) =
        terminalPhasor[terminal](0) * SHIFT_TO_PHASE_C;
  }
  const MatrixComp seriesAdmittance =
      ((**mSeriesRes).cast<Complex>() +
       Complex(0., 2. * PI * frequency) * (**mSeriesInd).cast<Complex>())
          .inverse();
  const MatrixComp shuntAdmittance =
      ((**mParallelCond).cast<Complex>() +
       Complex(0., 2. * PI * frequency) * (**mParallelCap).cast<Complex>()) /
      2.;
  const MatrixComp seriesCurrent =
      seriesAdmittance * (terminalPhasor[1] - terminalPhasor[0]);
  **mTerminalCurrent0 =
      (-seriesCurrent + shuntAdmittance * terminalPhasor[0]).real();
  **mTerminalCurrent1 =
      (seriesCurrent + shuntAdmittance * terminalPhasor[1]).real();
  if (!mUseSubcomponents) {
    const Real omega = 2. * PI * frequency;
    MatrixComp v0 = MatrixComp::Zero(3, 1), v1 = MatrixComp::Zero(3, 1);
    v0(0) = RMS3PH_TO_PEAK1PH * initialSingleVoltage(0);
    v1(0) = RMS3PH_TO_PEAK1PH * initialSingleVoltage(1);
    v0(1) = v0(0) * SHIFT_TO_PHASE_B;
    v0(2) = v0(0) * SHIFT_TO_PHASE_C;
    v1(1) = v1(0) * SHIFT_TO_PHASE_B;
    v1(2) = v1(0) * SHIFT_TO_PHASE_C;
    const MatrixComp voltage = v1 - v0;
    const MatrixComp impedance =
        (**mSeriesRes).cast<Complex>() +
        Complex(0., omega) * (**mSeriesInd).cast<Complex>();
    const MatrixComp current = impedance.inverse() * voltage;
    const MatrixComp halfCap = ((**mParallelCap) / 2.).cast<Complex>();
    **mIntfVoltage = voltage.real();
    **mIntfCurrent = current.real();
    if (mHasSeriesResistance)
      mVirtualNodes[0]->setInitialVoltage(
          PEAK1PH_TO_RMS3PH *
          (v0 + (**mSeriesRes).cast<Complex>() * current).real());
    (**mDirectState).setZero();
    (**mDirectState).block<3, 1>(0, 0) = current.real();
    (**mDirectState).block<3, 1>(3, 0) =
        (voltage - (**mSeriesRes).cast<Complex>() * current).real();
    (**mDirectState).block<3, 1>(6, 0) = v0.real();
    (**mDirectState).block<3, 1>(9, 0) =
        (Complex(0., omega) * halfCap * v0).real();
    (**mDirectState).block<3, 1>(12, 0) = v1.real();
    (**mDirectState).block<3, 1>(15, 0) =
        (Complex(0., omega) * halfCap * v1).real();
    return;
  }
  // Static calculation
  Real omega = 2. * PI * frequency;
  MatrixComp impedance = MatrixComp::Zero(3, 3);
  impedance << Complex((**mSeriesRes)(0, 0), omega * (**mSeriesInd)(0, 0)),
      Complex((**mSeriesRes)(0, 1), omega * (**mSeriesInd)(0, 1)),
      Complex((**mSeriesRes)(0, 2), omega * (**mSeriesInd)(0, 2)),
      Complex((**mSeriesRes)(1, 0), omega * (**mSeriesInd)(1, 0)),
      Complex((**mSeriesRes)(1, 1), omega * (**mSeriesInd)(1, 1)),
      Complex((**mSeriesRes)(1, 2), omega * (**mSeriesInd)(1, 2)),
      Complex((**mSeriesRes)(2, 0), omega * (**mSeriesInd)(2, 0)),
      Complex((**mSeriesRes)(2, 1), omega * (**mSeriesInd)(2, 1)),
      Complex((**mSeriesRes)(2, 2), omega * (**mSeriesInd)(2, 2));

  MatrixComp vInitABC = MatrixComp::Zero(3, 1);
  vInitABC(0, 0) = RMS3PH_TO_PEAK1PH * initialSingleVoltage(1) -
                   RMS3PH_TO_PEAK1PH * initialSingleVoltage(0);
  vInitABC(1, 0) = vInitABC(0, 0) * SHIFT_TO_PHASE_B;
  vInitABC(2, 0) = vInitABC(0, 0) * SHIFT_TO_PHASE_C;
  MatrixComp iInit = impedance.inverse() * vInitABC;
  **mIntfCurrent = iInit.real();
  **mIntfVoltage = vInitABC.real();

  // Initialization of virtual node
  // Initial voltage of phase B,C is set after A
  MatrixComp vInitTerm0 = MatrixComp::Zero(3, 1);
  vInitTerm0(0, 0) = RMS3PH_TO_PEAK1PH * initialSingleVoltage(0);
  vInitTerm0(1, 0) = vInitTerm0(0, 0) * SHIFT_TO_PHASE_B;
  vInitTerm0(2, 0) = vInitTerm0(0, 0) * SHIFT_TO_PHASE_C;

  if (mHasSeriesResistance)
    mVirtualNodes[0]->setInitialVoltage(PEAK1PH_TO_RMS3PH *
                                        (vInitTerm0 + **mSeriesRes * iInit));

  SPDLOG_LOGGER_DEBUG(mSLog,
                      "\n--debug--"
                      "\n seriesRes: {:s}"
                      "\n seriesInd: {:s}"
                      "\n Impedance: {:s}"
                      "\n vInit: {:s}"
                      "\n iInit: {:s}",
                      Logger::matrixToString(**mSeriesRes),
                      Logger::matrixToString(**mSeriesInd),
                      Logger::matrixCompToString(impedance),
                      Logger::matrixCompToString(vInitABC),
                      Logger::matrixCompToString(iInit));

  SPDLOG_LOGGER_INFO(
      mSLog,
      "\n--- Initialization from powerflow ---"
      "\nVoltage across: {:s}"
      "\nCurrent: {:s}"
      "\nTerminal 0 voltage: {:s}"
      "\nTerminal 1 voltage: {:s}"
      "\n--- Initialization from powerflow finished ---",
      Logger::matrixToString(**mIntfVoltage),
      Logger::matrixToString(**mIntfCurrent),
      Logger::phasorToString(RMS3PH_TO_PEAK1PH * initialSingleVoltage(0)),
      Logger::phasorToString(RMS3PH_TO_PEAK1PH * initialSingleVoltage(1)));
  mSLog->flush();
}

void EMT::Ph3::PiLine::mnaParentInitialize(Real, Real timeStep,
                                           Attribute<Matrix>::Ptr) {
  if (mUseSubcomponents)
    return;
  if (mHasSeriesResistance)
    mDirectResistanceConductance = (**mSeriesRes).inverse();
  mDirectInductorAdmittance = (timeStep / 2.) * (**mSeriesInd).inverse();
  mDirectCapAdmittance = (**mParallelCap) / timeStep;
}

void EMT::Ph3::PiLine::mnaParentApplySystemMatrixStamp(
    SparseMatrixRow &systemMatrix) {
  if (mUseSubcomponents)
    return;
  updateMatrixNodeIndices();
  const UInt midpoint = mHasSeriesResistance
                            ? mVirtualNodes[0]->matrixNodeIndex(PhaseType::A)
                            : matrixNodeIndex(0);
  if (mHasSeriesResistance)
    MNAStampUtils::stampConductanceMatrix(
        mDirectResistanceConductance, systemMatrix, matrixNodeIndex(0),
        midpoint, terminalNotGrounded(0), true, mSLog);
  MNAStampUtils::stampConductanceMatrix(
      mDirectInductorAdmittance, systemMatrix, midpoint, matrixNodeIndex(1),
      mHasSeriesResistance || terminalNotGrounded(0), terminalNotGrounded(1),
      mSLog);
  const Matrix shunt = mDirectCapAdmittance + (**mParallelCond) / 2.;
  for (UInt terminal = 0; terminal < 2; ++terminal)
    MNAStampUtils::stampConductanceMatrix(
        shunt, systemMatrix, matrixNodeIndex(terminal), 0,
        terminalNotGrounded(terminal), false, mSLog);
}

void EMT::Ph3::PiLine::mnaParentApplyRightSideVectorStamp(Matrix &rightVector) {
  if (mUseSubcomponents)
    return;
  if (mHasSeriesResistance) {
    const UInt midpoint = mVirtualNodes[0]->matrixNodeIndex(PhaseType::A);
    for (UInt phase = 0; phase < 3; ++phase)
      Math::setVectorElement(rightVector, midpoint + phase,
                             mDirectSeriesHistory(phase));
  }
  if (terminalNotGrounded(1))
    for (UInt phase = 0; phase < 3; ++phase)
      Math::setVectorElement(rightVector, matrixNodeIndex(1, phase),
                             -mDirectSeriesHistory(phase) -
                                 mDirectCapHistory1(phase));
  if (terminalNotGrounded(0))
    for (UInt phase = 0; phase < 3; ++phase)
      Math::setVectorElement(
          rightVector, matrixNodeIndex(0, phase),
          -mDirectCapHistory0(phase) +
              (mHasSeriesResistance ? 0. : mDirectSeriesHistory(phase)));
}

Bool EMT::Ph3::PiLine::mnaHasCompactRightSideVectorStamp() const {
  return !mUseSubcomponents ||
         CompositePowerComp<Real>::mnaHasCompactRightSideVectorStamp();
}

void EMT::Ph3::PiLine::mnaAddCompactRightSideVectorStamp(Matrix &rightVector) {
  if (mUseSubcomponents) {
    CompositePowerComp<Real>::mnaAddCompactRightSideVectorStamp(rightVector);
    return;
  }
  if (mHasSeriesResistance) {
    const UInt midpoint = mVirtualNodes[0]->matrixNodeIndex(PhaseType::A);
    for (UInt phase = 0; phase < 3; ++phase)
      Math::addToVectorElement(
          rightVector, midpoint + phase,
          Math::realFromVectorElement(**mRightVector, midpoint + phase));
  }
  for (UInt terminal = 0; terminal < 2; ++terminal) {
    if (!terminalNotGrounded(terminal))
      continue;
    for (UInt phase = 0; phase < 3; ++phase) {
      const UInt idx = matrixNodeIndex(terminal, phase);
      Math::addToVectorElement(
          rightVector, idx, Math::realFromVectorElement(**mRightVector, idx));
    }
  }
}

void EMT::Ph3::PiLine::mnaParentAddPreStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes) {
  if (!mUseSubcomponents)
    prevStepDependencies.push_back(mDirectState);
  else {
    prevStepDependencies.push_back(mIntfCurrent);
    prevStepDependencies.push_back(mIntfVoltage);
  }
  modifiedAttributes.push_back(mRightVector);
}

void EMT::Ph3::PiLine::mnaParentPreStep(Real time, Int timeStepCount) {
  PassiveOperationProfile::Scope operation(mOperationProfile,
                                           PassiveOperationProfile::RhsLocal);
  if (!mUseSubcomponents) {
    mDirectInductorHistory.noalias() =
        mDirectInductorAdmittance * (**mDirectState).block<3, 1>(3, 0);
    mDirectInductorHistory += (**mDirectState).block<3, 1>(0, 0);
    mDirectSeriesHistory = mDirectInductorHistory;
    mDirectCapHistory0 =
        -(**mDirectState).block<3, 1>(9, 0) -
        mDirectCapAdmittance * (**mDirectState).block<3, 1>(6, 0);
    mDirectCapHistory1 =
        -(**mDirectState).block<3, 1>(15, 0) -
        mDirectCapAdmittance * (**mDirectState).block<3, 1>(12, 0);
  }
  mnaCompApplyRightSideVectorStamp(**mRightVector);
}

void EMT::Ph3::PiLine::mnaParentAddPostStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes,
    Attribute<Matrix>::Ptr &leftVector) {
  attributeDependencies.push_back(leftVector);
  modifiedAttributes.push_back(mIntfVoltage);
  modifiedAttributes.push_back(mIntfCurrent);
  modifiedAttributes.push_back(mOutputUpdates);
  modifiedAttributes.push_back(mTerminalCurrent0);
  modifiedAttributes.push_back(mTerminalCurrent1);
  if (!mUseSubcomponents)
    modifiedAttributes.push_back(mDirectState);
}

void EMT::Ph3::PiLine::mnaParentPostStep(Real time, Int timeStepCount,
                                         Attribute<Matrix>::Ptr &leftVector) {
  PassiveOperationProfile::Scope operation(
      mOperationProfile, PassiveOperationProfile::VoltageRecovery);
  ++**mOutputUpdates;
  if (!mUseSubcomponents) {
    const auto voltage = [&](UInt terminal) {
      Eigen::Matrix<Real, 3, 1> result = Eigen::Matrix<Real, 3, 1>::Zero();
      if (terminalNotGrounded(terminal))
        for (UInt phase = 0; phase < 3; ++phase)
          result(phase) = Math::realFromVectorElement(
              **leftVector, matrixNodeIndex(terminal, phase));
      return result;
    };
    const auto v0 = voltage(0), v1 = voltage(1);
    Eigen::Matrix<Real, 3, 1> vmid = v0;
    if (mHasSeriesResistance)
      for (UInt phase = 0; phase < 3; ++phase)
        vmid(phase) = Math::realFromVectorElement(
            **leftVector,
            mVirtualNodes[0]->matrixNodeIndex(static_cast<PhaseType>(phase)));
    **mIntfVoltage = v1 - v0;
    operation.next(PassiveOperationProfile::StateAdvance);
    **mIntfCurrent =
        mDirectInductorAdmittance * (v1 - vmid) + mDirectSeriesHistory;
    (**mDirectState).block<3, 1>(0, 0) = **mIntfCurrent;
    (**mDirectState).block<3, 1>(3, 0) = v1 - vmid;
    (**mDirectState).block<3, 1>(6, 0) = v0;
    (**mDirectState).block<3, 1>(9, 0) =
        mDirectCapAdmittance * v0 + mDirectCapHistory0;
    (**mDirectState).block<3, 1>(12, 0) = v1;
    (**mDirectState).block<3, 1>(15, 0) =
        mDirectCapAdmittance * v1 + mDirectCapHistory1;
    operation.next(PassiveOperationProfile::Output);
    **mTerminalCurrent0 = -**mIntfCurrent + (**mDirectState).block<3, 1>(9, 0) +
                          (**mParallelCond / 2.) * v0;
    **mTerminalCurrent1 = **mIntfCurrent + (**mDirectState).block<3, 1>(15, 0) +
                          (**mParallelCond / 2.) * v1;
    return;
  }
  mnaCompUpdateVoltage(**leftVector);
  operation.next(PassiveOperationProfile::StateAdvance);
  mnaCompUpdateCurrent(**leftVector);
  operation.next(PassiveOperationProfile::Output);
  **mTerminalCurrent0 = -**mIntfCurrent;
  **mTerminalCurrent1 = **mIntfCurrent;
  if (mSubParallelResistor0) {
    **mTerminalCurrent0 += mSubParallelResistor0->intfCurrent();
    **mTerminalCurrent1 += mSubParallelResistor1->intfCurrent();
  }
  if (mSubParallelCapacitor0) {
    **mTerminalCurrent0 += mSubParallelCapacitor0->intfCurrent();
    **mTerminalCurrent1 += mSubParallelCapacitor1->intfCurrent();
  }
}

void EMT::Ph3::PiLine::mnaCompUpdateVoltage(const Matrix &leftVector) {
  // v1 - v0
  (**mIntfVoltage).setZero();
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

void EMT::Ph3::PiLine::mnaCompUpdateCurrent(const Matrix &leftVector) {
  if (!mUseSubcomponents) {
    Eigen::Matrix<Real, 3, 1> midpointVoltage =
        Eigen::Matrix<Real, 3, 1>::Zero();
    if (mHasSeriesResistance || terminalNotGrounded(0))
      for (UInt phase = 0; phase < 3; ++phase)
        midpointVoltage(phase) = Math::realFromVectorElement(
            leftVector, mHasSeriesResistance
                            ? mVirtualNodes[0]->matrixNodeIndex(
                                  static_cast<PhaseType>(phase))
                            : matrixNodeIndex(0, phase));
    Eigen::Matrix<Real, 3, 1> terminal1 = Eigen::Matrix<Real, 3, 1>::Zero();
    if (terminalNotGrounded(1))
      for (UInt phase = 0; phase < 3; ++phase)
        terminal1(phase) =
            Math::realFromVectorElement(leftVector, matrixNodeIndex(1, phase));
    **mIntfCurrent = mDirectInductorAdmittance * (terminal1 - midpointVoltage) +
                     mDirectSeriesHistory;
    return;
  }
  **mIntfCurrent = mSubSeriesInductor->intfCurrent();
}

// #### Tear Methods ####
MNAInterface::List EMT::Ph3::PiLine::mnaTearGroundComponents() {
  MNAInterface::List gndComponents;

  if (mSubParallelResistor0) {
    gndComponents.push_back(mSubParallelResistor0);
    gndComponents.push_back(mSubParallelResistor1);
  }

  if ((**mParallelCap)(0, 0) > 0) {
    gndComponents.push_back(mSubParallelCapacitor0);
    gndComponents.push_back(mSubParallelCapacitor1);
  }

  return gndComponents;
}

void EMT::Ph3::PiLine::mnaTearInitialize(Real omega, Real timeStep) {
  if (mSubSeriesResistor) {
    mSubSeriesResistor->mnaTearSetIdx(mTearIdx);
    mSubSeriesResistor->mnaTearInitialize(omega, timeStep);
  }
  mSubSeriesInductor->mnaTearSetIdx(mTearIdx);
  mSubSeriesInductor->mnaTearInitialize(omega, timeStep);
}

void EMT::Ph3::PiLine::mnaTearApplyMatrixStamp(SparseMatrixRow &tearMatrix) {
  if (mSubSeriesResistor)
    mSubSeriesResistor->mnaTearApplyMatrixStamp(tearMatrix);
  mSubSeriesInductor->mnaTearApplyMatrixStamp(tearMatrix);
}

void EMT::Ph3::PiLine::mnaTearApplyVoltageStamp(Matrix &voltageVector) {
  mSubSeriesInductor->mnaTearApplyVoltageStamp(voltageVector);
}

void EMT::Ph3::PiLine::mnaTearPostStep(MatrixComp voltage, MatrixComp current) {
  mSubSeriesInductor->mnaTearPostStep(voltage - (**mSeriesRes * current),
                                      current);
  (**mIntfCurrent) = mSubSeriesInductor->intfCurrent();
}
