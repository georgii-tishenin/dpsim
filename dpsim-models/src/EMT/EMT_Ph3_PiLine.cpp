/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim-models/EMT/EMT_Ph3_PiLine.h>

using namespace CPS;

EMT::Ph3::PiLine::PiLine(String uid, String name, Logger::Level logLevel)
	: Base::Ph3::PiLine(mAttributes), CompositePowerComp<Real>(uid, name, true, true, logLevel) {
	mPhaseType = PhaseType::ABC;
	setVirtualNodeNumber(1);
	setTerminalNumber(2);

	SPDLOG_LOGGER_INFO(mSLog, "Create {} {}", this->type(), name);
	**mIntfVoltage = Matrix::Zero(3, 1);
	**mIntfCurrent = Matrix::Zero(3, 1);

	mSLog->flush();
}

/// DEPRECATED: Delete method
SimPowerComp<Real>::Ptr EMT::Ph3::PiLine::clone(String name) {
	auto copy = PiLine::make(name, mLogLevel);
	copy->setParameters(**mSeriesRes, **mSeriesInd, **mParallelCap, **mParallelCond);
	return copy;
}

void EMT::Ph3::PiLine::initializeFromNodesAndTerminals(Real frequency) {

	// By default there is always a small conductance to ground to
	// avoid problems with floating nodes.
	Matrix defaultParallelCond = Matrix::Zero(3, 3);
	defaultParallelCond <<
		1e-6, 0, 0,
		0, 1e-6, 0,
		0, 0, 1e-6;
	**mParallelCond = ((**mParallelCond)(0, 0) > 0) ? **mParallelCond : defaultParallelCond;

	// Static calculation
	Real omega = 2. * PI * frequency;
	MatrixComp impedance = MatrixComp::Zero(3, 3);
	impedance <<
		Complex((**mSeriesRes)(0, 0), omega * (**mSeriesInd)(0, 0)), Complex((**mSeriesRes)(0, 1), omega * (**mSeriesInd)(0, 1)), Complex((**mSeriesRes)(0, 2), omega * (**mSeriesInd)(0, 2)),
		Complex((**mSeriesRes)(1, 0), omega * (**mSeriesInd)(1, 0)), Complex((**mSeriesRes)(1, 1), omega * (**mSeriesInd)(1, 1)), Complex((**mSeriesRes)(1, 2), omega * (**mSeriesInd)(1, 2)),
		Complex((**mSeriesRes)(2, 0), omega * (**mSeriesInd)(2, 0)), Complex((**mSeriesRes)(2, 1), omega * (**mSeriesInd)(2, 1)), Complex((**mSeriesRes)(2, 2), omega * (**mSeriesInd)(2, 2));

	MatrixComp vInitABC = MatrixComp::Zero(3, 1);
	vInitABC(0, 0) = RMS3PH_TO_PEAK1PH * initialSingleVoltage(1) - RMS3PH_TO_PEAK1PH * initialSingleVoltage(0);
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

	mVirtualNodes[0]->setInitialVoltage(PEAK1PH_TO_RMS3PH*(vInitTerm0 + **mSeriesRes * iInit));

	// Create series sub components
	mSubSeriesResistor = std::make_shared<EMT::Ph3::Resistor>(**mName + "_res", mLogLevel);
	mSubSeriesResistor->setParameters(**mSeriesRes);
	mSubSeriesResistor->connect({ mTerminals[0]->node(), mVirtualNodes[0] });
	mSubSeriesResistor->initialize(mFrequencies);
	mSubSeriesResistor->initializeFromNodesAndTerminals(frequency);
	addMNASubComponent(mSubSeriesResistor, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, false);

	mSubSeriesInductor = std::make_shared<EMT::Ph3::Inductor>(**mName + "_ind", mLogLevel);
	mSubSeriesInductor->setParameters(**mSeriesInd);
	mSubSeriesInductor->connect({ mVirtualNodes[0], mTerminals[1]->node() });
	mSubSeriesInductor->initialize(mFrequencies);
	mSubSeriesInductor->initializeFromNodesAndTerminals(frequency);
	addMNASubComponent(mSubSeriesInductor, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);

	// Create parallel sub components
	mSubParallelResistor0 = std::make_shared<EMT::Ph3::Resistor>(**mName + "_con0", mLogLevel);
	mSubParallelResistor0->setParameters(2. * (**mParallelCond).inverse());
	mSubParallelResistor0->connect(SimNode::List{ SimNode::GND, mTerminals[0]->node() });
	mSubParallelResistor0->initialize(mFrequencies);
	mSubParallelResistor0->initializeFromNodesAndTerminals(frequency);
	addMNASubComponent(mSubParallelResistor0, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, false);

	mSubParallelResistor1 = std::make_shared<EMT::Ph3::Resistor>(**mName + "_con1", mLogLevel);
	mSubParallelResistor1->setParameters(2. * (**mParallelCond).inverse());
	mSubParallelResistor1->connect(SimNode::List{ SimNode::GND, mTerminals[1]->node() });
	mSubParallelResistor1->initialize(mFrequencies);
	mSubParallelResistor1->initializeFromNodesAndTerminals(frequency);
	addMNASubComponent(mSubParallelResistor1, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, false);

	if ((**mParallelCap)(0,0) > 0) {
		mSubParallelCapacitor0 = std::make_shared<EMT::Ph3::Capacitor>(**mName + "_cap0", mLogLevel);
		mSubParallelCapacitor0->setParameters(**mParallelCap / 2.);
		mSubParallelCapacitor0->connect(SimNode::List{ SimNode::GND, mTerminals[0]->node() });
		mSubParallelCapacitor0->initialize(mFrequencies);
		mSubParallelCapacitor0->initializeFromNodesAndTerminals(frequency);
		addMNASubComponent(mSubParallelCapacitor0, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);

		mSubParallelCapacitor1 = std::make_shared<EMT::Ph3::Capacitor>(**mName + "_cap1", mLogLevel);
		mSubParallelCapacitor1->setParameters(**mParallelCap / 2.);
		mSubParallelCapacitor1->connect(SimNode::List{ SimNode::GND, mTerminals[1]->node() });
		mSubParallelCapacitor1->initialize(mFrequencies);
		mSubParallelCapacitor1->initializeFromNodesAndTerminals(frequency);
		addMNASubComponent(mSubParallelCapacitor1, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);
	}

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

	SPDLOG_LOGGER_INFO(mSLog, 
		"\n--- Initialization from powerflow ---"
		"\nVoltage across: {:s}"
		"\nCurrent: {:s}"
		"\nTerminal 0 voltage: {:s}"
		"\nTerminal 1 voltage: {:s}"
		"\nVirtual Node 1 voltage: {:s}"
		"\n--- Initialization from powerflow finished ---",
		Logger::matrixToString(**mIntfVoltage),
		Logger::matrixToString(**mIntfCurrent),
		Logger::phasorToString(RMS3PH_TO_PEAK1PH * initialSingleVoltage(0)),
		Logger::phasorToString(RMS3PH_TO_PEAK1PH * initialSingleVoltage(1)),
		Logger::phasorToString(mVirtualNodes[0]->initialSingleVoltage()));
	mSLog->flush();
}

void EMT::Ph3::PiLine::mnaParentAddPreStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes){
	prevStepDependencies.push_back(mIntfCurrent);
	prevStepDependencies.push_back(mIntfVoltage);
	modifiedAttributes.push_back(mRightVector);
}

void EMT::Ph3::PiLine::mnaParentPreStep(Real time, Int timeStepCount) {
	mnaCompApplyRightSideVectorStamp(**mRightVector);
}

void EMT::Ph3::PiLine::mnaParentAddPostStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes, Attribute<Matrix>::Ptr &leftVector) {
	attributeDependencies.push_back(leftVector);
	modifiedAttributes.push_back(mIntfVoltage);
	modifiedAttributes.push_back(mIntfCurrent);
}

void EMT::Ph3::PiLine::mnaParentPostStep(Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector) {
	mnaCompUpdateVoltage(**leftVector);
	mnaCompUpdateCurrent(**leftVector);
}

void EMT::Ph3::PiLine::mnaCompUpdateVoltage(const Matrix& leftVector) {
	// v1 - v0
	**mIntfVoltage = Matrix::Zero(3, 1);
	if (terminalNotGrounded(1)) {
		(**mIntfVoltage)(0, 0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(1, 0));
		(**mIntfVoltage)(1, 0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(1, 1));
		(**mIntfVoltage)(2, 0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(1, 2));
	}
	if (terminalNotGrounded(0)) {
		(**mIntfVoltage)(0, 0) = (**mIntfVoltage)(0, 0) - Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 0));
		(**mIntfVoltage)(1, 0) = (**mIntfVoltage)(1, 0) - Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 1));
		(**mIntfVoltage)(2, 0) = (**mIntfVoltage)(2, 0) - Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 2));
	}
}

void EMT::Ph3::PiLine::mnaCompUpdateCurrent(const Matrix& leftVector) {
	**mIntfCurrent = mSubSeriesInductor->intfCurrent();
}

// #### DAE functions ####

void EMT::Ph3::PiLine::daeInitialize(double time, double state[], double dstate_dt[], 
	double absoluteTolerances[], double stateVarTypes[], int& offset) {
	// state variables are: 3xinductor_current

	updateMatrixNodeIndices();
	int index_node00 = matrixNodeIndex(0, 0);
	int index_node01 = matrixNodeIndex(0, 1);
	int index_node02 = matrixNodeIndex(0, 2);
	int index_node10 = matrixNodeIndex(1, 0);
	int index_node11 = matrixNodeIndex(1, 1);
	int index_node12 = matrixNodeIndex(1, 2);

	// initialize inductor variables
	Matrix inductorVoltage = mSubSeriesInductor->intfVoltage();
	state[offset] = (**mIntfCurrent)(0,0);
	dstate_dt[offset]   = inductorVoltage(0,0) / (**mSeriesInd)(0,0);
	state[offset+1] = (**mIntfCurrent)(1,0);
	dstate_dt[offset+1] = inductorVoltage(1,0) / (**mSeriesInd)(1,1);
	state[offset+2] = (**mIntfCurrent)(2,0);
	dstate_dt[offset+2] = inductorVoltage(2,0) / (**mSeriesInd)(2, 2);

	//set inductor state variables as differential variable
	stateVarTypes[offset+0] = 0.0;
	stateVarTypes[offset+1] = 0.0;
	stateVarTypes[offset+2] = 0.0;

	//set absolute tolerance =
	absoluteTolerances[offset]   = mAbsTolerance;
	absoluteTolerances[offset+1] = mAbsTolerance;
	absoluteTolerances[offset+2] = mAbsTolerance;

	mSLog->info(
		"\n--- daeInitialize ---"
		"\nAdded current-phase1 through the inductor of PiLine '{:s}' to state vector, initial value={:f}A"
		"\nAdded current-phase2 through the inductor of PiLine '{:s}' to state vector, initial value={:f}A"
		"\nAdded current-phase3 through the inductor of PiLine '{:s}' to state vector, initial value={:f}A"
		"\nAdded derivative of current-phase1 through the inductor of PiLine '{:s}' to derivative state vector, initial value={:f}"
		"\nAdded derivative of current-phase2 through the inductor of PiLine '{:s}' to derivative state vector, initial value={:f}"
		"\nAdded derivative of current-phase3 through the inductor of PiLine '{:s}' to derivative state vector, initial value={:f}"
		"\nState variables of inductors set as differential"
		"\nAbsolute tolerances={:f}"	
		"\n--- daeInitialize finished ---",

		this->name(), state[offset],
		this->name(), state[offset+1],
		this->name(), state[offset+2],
		this->name(), dstate_dt[offset],
		this->name(), dstate_dt[offset+1],
		this->name(), dstate_dt[offset+2],
		absoluteTolerances[offset]
	);
	
	mSLog->flush();
	offset+=3;
}

void EMT::Ph3::PiLine::daeResidual(double sim_time,
	const double state[], const double dstate_dt[],
	double resid[], std::vector<int>& off) {
	// offset+3, offset+4, offset+5 --> cap0 --> node1	(left node)
	// offset+3, offset+4, offset+5 --> cap1 --> node0  (right node)

	int c_offset = off[0]+off[1]; //current offset for component
	int pos_node00 = matrixNodeIndex(0, 0);
	int pos_node01 = matrixNodeIndex(0, 1);
	int pos_node02 = matrixNodeIndex(0, 2);
	int pos_node10 = matrixNodeIndex(1, 0);
	int pos_node11 = matrixNodeIndex(1, 1);
	int pos_node12 = matrixNodeIndex(1, 2);

	// residual function of inductors: v1-(v0+vr) - L*di(t)/dt = 0
	resid[c_offset]    = -((**mSeriesInd)(0, 0) * dstate_dt[c_offset+0] + state[c_offset+0] * (**mSeriesRes)(0, 0));
	resid[c_offset+1]  = -((**mSeriesInd)(1, 1) * dstate_dt[c_offset+1] + state[c_offset+1] * (**mSeriesRes)(1, 1));
	resid[c_offset+2]  = -((**mSeriesInd)(2, 2) * dstate_dt[c_offset+2] + state[c_offset+2] * (**mSeriesRes)(2, 2));
	if (terminalNotGrounded(0)) {
		resid[c_offset]   -= state[pos_node00];
		resid[c_offset+1] -= state[pos_node01];
		resid[c_offset+2] -= state[pos_node02];

		// update residual equations of nodes
		resid[pos_node00] -= state[c_offset];
		resid[pos_node01] -= state[c_offset+1];
		resid[pos_node02] -= state[c_offset+2];
	}
	if (terminalNotGrounded(1)) {
		resid[c_offset]   += state[pos_node10];
		resid[c_offset+1] += state[pos_node11];
		resid[c_offset+2] += state[pos_node12];

		// update residual equations of nodes
		resid[pos_node10] += state[c_offset];
		resid[pos_node11] += state[c_offset+1];
		resid[pos_node12] += state[c_offset+2];
	}

	//add cap and cond currents to nodal equations
	Real conductance = (**mParallelCond)(0,0)/2;
	Real capacitance = (**mParallelCap)(0,0)/2;
	resid[pos_node10] += dstate_dt[pos_node10] * capacitance + state[pos_node10] * conductance;
	resid[pos_node11] += dstate_dt[pos_node11] * capacitance + state[pos_node11] * conductance;
	resid[pos_node12] += dstate_dt[pos_node12] * capacitance + state[pos_node12] * conductance;
	resid[pos_node00] += dstate_dt[pos_node00] * capacitance + state[pos_node00] * conductance;
	resid[pos_node01] += dstate_dt[pos_node01] * capacitance + state[pos_node01] * conductance;
	resid[pos_node02] += dstate_dt[pos_node02] * capacitance + state[pos_node02] * conductance;
	
	mSLog->debug(
		"\n\n--- daeResidual PiLoad - name: {:s} SimStep = {:f} ---"
		"\nresid[c_offset]   = v1 - (v0+vr) - v_l(t) = state[pos_node10] - (state[pos_node00]+state[c_offset+0]*mSeriesRes(0, 0)) -mSeriesInd(0, 0)*dstate_dt[c_offset+0] = {:f} - ({:f}+{:f}*{:f}) - {:f}*{:f} = {:f}"
		"\nresid[c_offset]   = v1 - (v0+vr) - v_l(t) = state[pos_node11] - (state[pos_node01]+state[c_offset+1]*mSeriesRes(1, 1)) -mSeriesInd(1, 1)*dstate_dt[c_offset+1] = {:f} - ({:f}+{:f}*{:f}) - {:f}*{:f} = {:f}"
		"\nresid[c_offset]   = v1 - (v0+vr) - v_l(t) = state[pos_node12] - (state[pos_node02]+state[c_offset+2]*mSeriesRes(2, 2)) -mSeriesInd(2, 2)*dstate_dt[c_offset+2] = {:f} - ({:f}+{:f}*{:f}) - {:f}*{:f} = {:f}"

		"\nupdate nodal equations:"
		"\nresid[pos_node10] += state[c_offset]   + state[pos_node10]*mConductance(0,0) + dstate_dt[pos_node10]*capacitance = {:f} + {:f}*{:f} + {:f}*{:f} = {:f}"
		"\nresid[pos_node11] += state[c_offset+1] + state[pos_node11]*mConductance(1,1) + dstate_dt[pos_node11]*capacitance = {:f} + {:f}*{:f} + {:f}*{:f} = {:f}"
		"\nresid[pos_node12] += state[c_offset+2] + state[pos_node12]*mConductance(2,2) + dstate_dt[pos_node12]*capacitance = {:f} + {:f}*{:f} + {:f}*{:f} = {:f}"
		"\n*** state[pos_node10]*mConductance(0,0)={:f}"
		"\n*** state[pos_node10]*capacitance(0,0)={:f}",

		this->name(), sim_time,
		state[pos_node10], state[pos_node00], state[c_offset],   (**mSeriesRes)(0, 0), (**mSeriesInd)(0,0), dstate_dt[c_offset],  resid[c_offset],
		state[pos_node11], state[pos_node01], state[c_offset+1], (**mSeriesRes)(1, 1), (**mSeriesInd)(1,1), dstate_dt[c_offset+1], resid[c_offset+1],
		state[pos_node12], state[pos_node02], state[c_offset+2], (**mSeriesRes)(2, 2), (**mSeriesInd)(2,2), dstate_dt[c_offset+2], resid[c_offset+2],

		state[c_offset],   state[pos_node10], conductance, dstate_dt[pos_node10], capacitance, resid[pos_node10],
		state[c_offset+1], state[pos_node11], conductance, dstate_dt[pos_node11], capacitance, resid[pos_node11],
		state[c_offset+2], state[pos_node12], conductance, dstate_dt[pos_node12], capacitance, resid[pos_node12],
		conductance * state[pos_node10],
		capacitance * dstate_dt[pos_node10]
	);
	mSLog->flush();

	off[1] += 3;

}

void EMT::Ph3::PiLine::daePostStep(double Nexttime, const double state[], const double dstate_dt[], int& offset) {
	(**mIntfVoltage)(0,0) = state[matrixNodeIndex(1, 0)] - state[matrixNodeIndex(0, 0)];
	(**mIntfVoltage)(1,0) = state[matrixNodeIndex(1, 1)] - state[matrixNodeIndex(0, 1)];
	(**mIntfVoltage)(2,0) = state[matrixNodeIndex(1, 2)] - state[matrixNodeIndex(0, 2)];
	(**mIntfCurrent)(0, 0) = state[offset];
	(**mIntfCurrent)(1, 0) = state[offset+1];
	(**mIntfCurrent)(2, 0) = state[offset+2];

	offset+=3;
}


