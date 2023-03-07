/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim-models/EMT/EMT_Ph1_Switch.h>

using namespace CPS;

EMT::Ph1::Switch::Switch(String uid, String name, Logger::Level logLevel)
	: SimPowerComp<Real>(uid, name, logLevel) {
	setTerminalNumber(2);
	**mIntfVoltage = Matrix::Zero(1,1);
	**mIntfCurrent = Matrix::Zero(1,1);

	mOpenResistance = Attribute<Real>::create("R_open", mAttributes);
	mClosedResistance = Attribute<Real>::create("R_closed", mAttributes);
	mIsClosed = Attribute<Bool>::create("is_closed", mAttributes);
}

SimPowerComp<Real>::Ptr EMT::Ph1::Switch::clone(String name) {
	auto copy = Switch::make(name, mLogLevel);
	copy->setParameters(**mOpenResistance, **mClosedResistance, **mIsClosed);
	return copy;
}

void EMT::Ph1::Switch::initializeFromNodesAndTerminals(Real frequency) {

	Real resistance = (**mIsClosed) ? **mClosedResistance : **mOpenResistance;

	(**mIntfVoltage)(0,0) = (initialSingleVoltage(1) - initialSingleVoltage(0)).real();
	(**mIntfCurrent)(0,0) = ((**mIntfVoltage)(0,0) / resistance);

	mSLog->info(
		"\n--- Initialization from powerflow ---"
		"\nVoltage across: {:f}"
		"\nCurrent: {:f}"
		"\nTerminal 0 voltage: {:f}"
		"\nTerminal 1 voltage: {:f}"
		"\n--- Initialization from powerflow finished ---",
		(**mIntfVoltage)(0,0),
		(**mIntfCurrent)(0,0),
		initialSingleVoltage(0).real(),
		initialSingleVoltage(1).real());
}

void EMT::Ph1::Switch::mnaInitialize(Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
	MNAInterface::mnaInitialize(omega, timeStep);
	updateMatrixNodeIndices();

	mMnaTasks.push_back(std::make_shared<MnaPostStep>(*this, leftVector));
}

Bool EMT::Ph1::Switch::mnaIsClosed() { return **mIsClosed; }

void EMT::Ph1::Switch::mnaApplySystemMatrixStamp(Matrix& systemMatrix) {
	Matrix conductance = Matrix::Zero(1, 1); 
	conductance(0, 0) = (**mIsClosed) ?
		1./(**mClosedResistance) : 1./(**mOpenResistance);

	// Set diagonal entries
	if (terminalNotGrounded(0)) {
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 0), matrixNodeIndex(0, 0), conductance(0, 0));
	}
	if (terminalNotGrounded(1)) {
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 0), matrixNodeIndex(1, 0), conductance(0, 0));
	}
	if (terminalNotGrounded(0) && terminalNotGrounded(1)) {
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 0), matrixNodeIndex(1, 0), -conductance(0, 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 0), matrixNodeIndex(0, 0), -conductance(0, 0));
	}
	SPDLOG_LOGGER_TRACE(mSLog, 
		"\nConductance matrix: {:s}",
		Logger::matrixToString(conductance));
}

void EMT::Ph1::Switch::mnaApplySwitchSystemMatrixStamp(Bool closed, Matrix& systemMatrix, Int freqIdx) {
	Matrix conductance = Matrix::Zero(1, 1); 
	conductance(0, 0) = (closed) ?
		1./(**mClosedResistance) : 1./(**mOpenResistance);

	// Set diagonal entries
	if (terminalNotGrounded(0)) {
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 0), matrixNodeIndex(0, 0), conductance(0, 0));

	}
	if (terminalNotGrounded(1)) {
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 0), matrixNodeIndex(1, 0), conductance(0, 0));
	}
	// Set off diagonal blocks, 2x3x3 entries
	if (terminalNotGrounded(0) && terminalNotGrounded(1)) {
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 0), matrixNodeIndex(1, 0), -conductance(0, 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 0), matrixNodeIndex(0, 0), -conductance(0, 0));
	}

	SPDLOG_LOGGER_TRACE(mSLog, 
		"\nConductance matrix: {:s}",
		Logger::matrixToString(conductance));
}

void EMT::Ph1::Switch::mnaApplyRightSideVectorStamp(Matrix& rightVector) { }



void EMT::Ph1::Switch::mnaUpdateVoltage(const Matrix& leftVector) {
	// Voltage across component is defined as V1 - V0
	**mIntfVoltage = Matrix::Zero(1, 1);
	if (terminalNotGrounded(1)) {
		(**mIntfVoltage)(0, 0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(1, 0));
	}
	if (terminalNotGrounded(0)) {
		(**mIntfVoltage)(0, 0) = (**mIntfVoltage)(0, 0) - Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 0));
	}
}

void EMT::Ph1::Switch::mnaUpdateCurrent(const Matrix& leftVector) {
	(**mIntfCurrent)(0, 0) = (**mIsClosed) ?
		(**mIntfVoltage)(0, 0)/(**mClosedResistance):
		(**mIntfVoltage)(0, 0)/(**mOpenResistance);
}

void EMT::Ph1::Switch::mnaAddPostStepDependencies(AttributeBase::List &prevStepDependencies,
	AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes,
	Attribute<Matrix>::Ptr &leftVector) {

	attributeDependencies.push_back(leftVector);
	modifiedAttributes.push_back(mIntfVoltage);
	modifiedAttributes.push_back(mIntfCurrent);
}

void EMT::Ph1::Switch::mnaPostStep(Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector) {
	mnaUpdateVoltage(**leftVector);
	mnaUpdateCurrent(**leftVector);
}