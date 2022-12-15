/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim-models/EMT/EMT_Ph3_NetworkInjection.h>

using namespace CPS;

EMT::Ph3::NetworkInjection::NetworkInjection(String uid, String name, Logger::Level logLevel)
	: CompositePowerComp<Real>(uid, name, true, true, logLevel),
	mVoltageRef(mAttributes->createDynamic<MatrixComp>("V_ref")),
	mSrcFreq(mAttributes->createDynamic<Real>("f_src")) {
	mPhaseType = PhaseType::ABC;
	setVirtualNodeNumber(0);
	setTerminalNumber(1);
	**mIntfVoltage = Matrix::Zero(3, 1);
	**mIntfCurrent = Matrix::Zero(3, 1);

	SPDLOG_LOGGER_INFO(mSLog, "Create {} {}", this->type(), name);

	// Create electrical sub components
	mSubVoltageSource = std::make_shared<EMT::Ph3::VoltageSource>(**mName + "_vs", mLogLevel);
	addMNASubComponent(mSubVoltageSource, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);
	SPDLOG_LOGGER_INFO(mSLog, "Electrical subcomponents: ");
	for (auto subcomp: mSubComponents)
		SPDLOG_LOGGER_INFO(mSLog, "- {}", subcomp->name());

	mSubVoltageSource->mVoltageRef->setReference(mVoltageRef);
	mSubVoltageSource->mSrcFreq->setReference(mSrcFreq);
}

SimPowerComp<Real>::Ptr EMT::Ph3::NetworkInjection::clone(String name) {
	auto copy = NetworkInjection::make(name, mLogLevel);
	copy->setParameters(**mVoltageRef);
	return copy;
}

void EMT::Ph3::NetworkInjection::setParameters(MatrixComp voltageRef, Real srcFreq) {
	mParametersSet = true;

	mSubVoltageSource->setParameters(voltageRef, srcFreq);

	SPDLOG_LOGGER_INFO(mSLog, "\nVoltage Ref={:s} [V]"
				"\nFrequency={:s} [Hz]",
				Logger::matrixCompToString(voltageRef),
				Logger::realToString(srcFreq));
}

void EMT::Ph3::NetworkInjection::setParameters(MatrixComp voltageRef, Real freqStart, Real rocof, Real timeStart, Real duration, bool smoothRamp) {
	mParametersSet = true;

	mSubVoltageSource->setParameters(voltageRef, freqStart, rocof, timeStart, duration, smoothRamp);

	SPDLOG_LOGGER_INFO(mSLog, "\nVoltage Ref={:s} [V]"
				"\nFrequency={:s} [Hz]",
				Logger::matrixCompToString(voltageRef),
				Logger::realToString(freqStart));
}

void EMT::Ph3::NetworkInjection::setParameters(MatrixComp voltageRef, Real modulationFrequency, Real modulationAmplitude, Real baseFrequency /*= 0.0*/, bool zigzag /*= false*/) {
	mParametersSet = true;

	mSubVoltageSource->setParameters(voltageRef, modulationFrequency, modulationAmplitude, baseFrequency, zigzag);

	SPDLOG_LOGGER_INFO(mSLog, "\nVoltage Ref={:s} [V]"
				"\nFrequency={:s} [Hz]",
				Logger::matrixCompToString(voltageRef),
				Logger::realToString(baseFrequency));
}

void EMT::Ph3::NetworkInjection::initializeFromNodesAndTerminals(Real frequency) {
	
	// Connect electrical subcomponents
	mSubVoltageSource->connect({ SimNode::GND, node(0) });

	// Initialize electrical subcomponents
	for (auto subcomp: mSubComponents) {
		subcomp->initialize(mFrequencies);
		subcomp->initializeFromNodesAndTerminals(frequency);
	}
}

// #### MNA functions ####
void EMT::Ph3::NetworkInjection::mnaParentApplyRightSideVectorStamp(Matrix& rightVector) {
	SPDLOG_LOGGER_DEBUG(mSLog, "Right Side Vector: {:s}",
				Logger::matrixToString(rightVector));
}

void EMT::Ph3::NetworkInjection::mnaParentAddPreStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes) {
	prevStepDependencies.push_back(mIntfCurrent);
	prevStepDependencies.push_back(mIntfVoltage);
	modifiedAttributes.push_back(mRightVector);
}

void EMT::Ph3::NetworkInjection::mnaParentPreStep(Real time, Int timeStepCount) {
	mnaCompApplyRightSideVectorStamp(**mRightVector);
}

void EMT::Ph3::NetworkInjection::mnaParentAddPostStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes, Attribute<Matrix>::Ptr &leftVector) {
	attributeDependencies.push_back(leftVector);
	modifiedAttributes.push_back(mIntfVoltage);
	modifiedAttributes.push_back(mIntfCurrent);
}

void EMT::Ph3::NetworkInjection::mnaParentPostStep(Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector) {
	mnaCompUpdateCurrent(**leftVector);
	mnaCompUpdateVoltage(**leftVector);
}

void EMT::Ph3::NetworkInjection::mnaCompUpdateVoltage(const Matrix& leftVector) {
	**mIntfVoltage = **mSubVoltageSource->mIntfVoltage;
}

void EMT::Ph3::NetworkInjection::mnaCompUpdateCurrent(const Matrix& leftVector) {
	**mIntfCurrent = **mSubVoltageSource->mIntfCurrent;
}

// #### DAE functions ####

/// set init value of current, calculate and set initial value of the derivative of the current
void EMT::Ph3::NetworkInjection::setInitialComplexIntfCurrent(Complex initCurrent) {
	mSubVoltageSource->setInitialComplexIntfCurrent(initCurrent);
	**mIntfCurrent = mSubVoltageSource->intfCurrent();
	mIntfDerCurrent = mSubVoltageSource->intfDerCurrent();
}

void EMT::Ph3::NetworkInjection::daeInitialize(double time, double state[],double dstate_dt[], 
	double absoluteTolerances[], double stateVarTypes[], int& offset) {
	// state[offset] = current through component. It is positive when it flows out of the component
	// dstate_dt[offset] = derivative of current through voltage source

	updateMatrixNodeIndices();
	//daePreStep(time);

	state[offset] = (**mIntfCurrent)(0,0);
	dstate_dt[offset] = mIntfDerCurrent(0,0);
	state[offset+1] = (**mIntfCurrent)(1,0);
	dstate_dt[offset+1] = mIntfDerCurrent(1,0);
	state[offset+2] = (**mIntfCurrent)(2,0);
	dstate_dt[offset+2] = mIntfDerCurrent(2,0);

	//set state variable as algebraic variable
	stateVarTypes[offset]   = 0.0;
	stateVarTypes[offset+1] = 0.0;
	stateVarTypes[offset+2] = 0.0;

	//set absolute tolerance
	absoluteTolerances[offset] = mAbsTolerance;
	absoluteTolerances[offset+1] = mAbsTolerance;
	absoluteTolerances[offset+2] = mAbsTolerance;

	mSLog->info(
		"\n--- daeInitialize ---"
		"\nInitial time={:f}s"
		"\nAdded current phase1 of NetworkInjection '{:s}' to state vector, initial value={:f}A"
		"\nAdded current phase2 of NetworkInjection '{:s}' to state vector, initial value={:f}A"
		"\nAdded current phase3 of NetworkInjection '{:s}' to state vector, initial value={:f}A"
		"\nAdded derivative of current phase1 of NetworkInjection '{:s}' to derivative state vector, initial value={:f}"
		"\nAdded derivative of current phase2 of NetworkInjection '{:s}' to derivative state vector, initial value={:f}"
		"\nAdded derivative of current phase3 of NetworkInjection '{:s}' to derivative state vector, initial value={:f}"
		"\nInitial voltage phase1 of NetworkInjection '{:s}' ={:f}V"
		"\nInitial voltage phase2 of NetworkInjection '{:s}' ={:f}V"
		"\nInitial voltage phase3 of NetworkInjection '{:s}' ={:f}V"
		"\nState variable set as differential"
		"\nAbsolute tolerance={:f}"
		"\n--- daeInitialize finished ---",
		time,
		this->name(), state[offset],
		this->name(), state[offset+1],
		this->name(), state[offset+2],
		this->name(), dstate_dt[offset],
		this->name(), dstate_dt[offset+1],
		this->name(), dstate_dt[offset+2],
		this->name(), (**mIntfVoltage)(0,0),
		this->name(), (**mIntfVoltage)(1,0),
		this->name(), (**mIntfVoltage)(2,0),
		absoluteTolerances[offset]
	);
	mSLog->flush();
	offset+=3;
}

void EMT::Ph3::NetworkInjection::daeResidual(double sim_time,
	const double state[], const double dstate_dt[],
	double resid[], std::vector<int>& off) {
	
	// current offset for component
	int c_offset = off[0]+off[1];

	// update voltage
	mSubVoltageSource->updateVoltage(sim_time);
	**mIntfVoltage = mSubVoltageSource->intfVoltage();
	//mIntfVoltage = mSubVoltageSource->attribute<Matrix>("v_intf")->get();

	// residual function of component (v_node - v_component = 0)
	resid[c_offset]   = state[matrixNodeIndex(0, 0)] - (**mIntfVoltage)(0,0);
	resid[c_offset+1] = state[matrixNodeIndex(0, 1)] - (**mIntfVoltage)(1,0);
	resid[c_offset+2] = state[matrixNodeIndex(0, 2)] - (**mIntfVoltage)(2,0);

	// add current to nodal equations
	resid[matrixNodeIndex(0, 0)] += state[c_offset];
	resid[matrixNodeIndex(0, 1)] += state[c_offset+1];
	resid[matrixNodeIndex(0, 2)] += state[c_offset+2];

	mSLog->flush();
	off[1]+=3;
}

void EMT::Ph3::NetworkInjection::daeJacobian(double current_time, const double state[], 
	const double dstate_dt[], SUNMatrix jacobian, double cj, std::vector<int>& off) {

	// current offset for component
	int c_offset = off[0] + off[1]; 

	SM_ELEMENT_D(jacobian, c_offset,   matrixNodeIndex(0, 0)) += 1.0;
	SM_ELEMENT_D(jacobian, c_offset+1, matrixNodeIndex(0, 1)) += 1.0;
	SM_ELEMENT_D(jacobian, c_offset+2, matrixNodeIndex(0, 2)) += 1.0;

	SM_ELEMENT_D(jacobian, matrixNodeIndex(0, 0), c_offset)   += 1.0;
	SM_ELEMENT_D(jacobian, matrixNodeIndex(0, 1), c_offset+1) += 1.0;
	SM_ELEMENT_D(jacobian, matrixNodeIndex(0, 2), c_offset+2) += 1.0;

	off[1] += 3;
}

void EMT::Ph3::NetworkInjection::daePostStep(double Nexttime, const double state[], 
	const double dstate_dt[], int& offset) {

	// update current
	(**mIntfCurrent)(0,0) = state[offset];
	(**mIntfCurrent)(1,0) = state[offset+1];
	(**mIntfCurrent)(2,0) = state[offset+2];

	// update voltage
	(**mIntfVoltage)(0,0) = state[matrixNodeIndex(0, 0)];
	(**mIntfVoltage)(1,0) = state[matrixNodeIndex(0, 1)];
	(**mIntfVoltage)(2,0) = state[matrixNodeIndex(0, 2)];

	// update current derivative
	mIntfDerCurrent(0,0) = dstate_dt[offset];
	mIntfDerCurrent(1,0) = dstate_dt[offset+1];
	mIntfDerCurrent(2,0) = dstate_dt[offset+2];

	offset +=3;
}
