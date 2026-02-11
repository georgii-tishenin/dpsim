/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/
#pragma once

#include <dpsim-models/CompositePowerComp.h>
#include <dpsim-models/Solver/MNAInterface.h>
#include <dpsim-models/Definitions.h>
#include <dpsim-models/EMT/EMT_Ph3_Resistor.h>
#include <dpsim-models/EMT/EMT_Ph3_Inductor.h>
#include <dpsim-models/EMT/EMT_Ph3_Capacitor.h>
#include <dpsim-models/EMT/EMT_Ph3_VoltageSource.h>
#include <dpsim-models/EMT/EMT_Ph3_Transformer.h>
#include <dpsim-models/Base/Base_AvVoltageSourceInverterDQ.h>
#include <dpsim-models/Signal/VCO.h>
#include <dpsim-models/Signal/VoltageControllerVSI.h>

namespace CPS {
namespace EMT {
namespace Ph3 {

class VSIVoltageControlDQ :
	public CompositePowerComp<Real>,
	public Base::AvVoltageSourceInverterDQ,
	public SharedFactory<VSIVoltageControlDQ> {

protected:
	// ### General Parameters ###
	/// Nominal voltage (unused, kept for compatibility)
	Real mVnom = 0;
	/// Simulation step
	Real mTimeStep = 0;

	// ### Control Subcomponents ###
	/// VCO (free running at omegaNom; no droop)
	std::shared_ptr<Signal::VCO> mVCO;
	/// Voltage Controller
	std::shared_ptr<Signal::VoltageControllerVSI> mVoltageControllerVSI;

	// ### Electrical Subcomponents ###
	/// Controlled voltage source
	std::shared_ptr<EMT::Ph3::VoltageSource> mSubCtrledVoltageSource;
	/// Resistor Rf as part of LCL filter
	std::shared_ptr<EMT::Ph3::Resistor> mSubResistorF;
	/// Capacitor Cf as part of LCL filter
	std::shared_ptr<EMT::Ph3::Capacitor> mSubCapacitorF;
	/// Inductor Lf as part of LCL filter
	std::shared_ptr<EMT::Ph3::Inductor> mSubInductorF;
	/// Resistor Rc as part of LCL filter
	std::shared_ptr<EMT::Ph3::Resistor> mSubResistorC;
	/// Optional connection transformer
	std::shared_ptr<EMT::Ph3::Transformer> mConnectionTransformer;

	/// Flag for connection transformer usage
	Bool mWithConnectionTransformer = false;
	/// Flag for controller usage
	Bool mWithControl = true;

	// #### solver ####
	std::vector<const Matrix*> mRightVectorStamps;

public:
	// ### General Parameters ###
	/// Nominal frequency
	const Attribute<Real>::Ptr mOmegaN;
	/// Voltage d reference
	const Attribute<Real>::Ptr mVdRef;
	/// Voltage q reference
	const Attribute<Real>::Ptr mVqRef;
	/// Active power reference
	const Attribute<Real>::Ptr mPRef;

	// ### Inverter Interfacing Variables ###
	// Control inputs
	const Attribute<Real>::Ptr mVcd;
	const Attribute<Real>::Ptr mVcq;
	const Attribute<Real>::Ptr mIrcd;
	const Attribute<Real>::Ptr mIrcq;

	// Electrical power logging
	const Attribute<Real>::Ptr mElecActivePower;
	const Attribute<Real>::Ptr mElecPassivePower;

	// Control outputs
	const Attribute<Matrix>::Ptr mVsref;

	// Sub voltage source logging
	const Attribute<Matrix>::Ptr mVs;

	// VCO output logging
	const Attribute<Real>::Ptr mVCOOutput;

	// input, state and output vector for logging
	const Attribute<Matrix>::Ptr mVoltagectrlInputs;
	const Attribute<Matrix>::Ptr mVoltagectrlStates;
	const Attribute<Matrix>::Ptr mVoltagectrlOutputs;

	VSIVoltageControlDQ(String name, Logger::Level logLevel = Logger::Level::off)
		: VSIVoltageControlDQ(name, name, logLevel) {}

	VSIVoltageControlDQ(String uid, String name, Logger::Level logLevel = Logger::Level::off, Bool withTrafo = false);

	// #### General ####
	void initializeFromNodesAndTerminals(Real frequency);
	void setParameters(Real Omega, Real VdRef, Real VqRef, Real Pref);

	// Keep the old signature for compatibility; droop parameters are ignored.
	void setControllerParameters(Real Kp_voltageCtrl, Real Ki_voltageCtrl,
	                             Real Kp_currCtrl, Real Ki_currCtrl,
	                             Real Omega, Real tau_p, Real tau_i, Real m_p);

	void setTransformerParameters(Real nomVoltageEnd1, Real nomVoltageEnd2, Real ratedPower,
	                              Real ratioAbs, Real ratioPhase, Real resistance, Real inductance, Real omega);

	void setFilterParameters(Real Lf, Real Cf, Real Rf, Real Rc);

	void setInitialStateValues(Real phi_dInit, Real phi_qInit, Real gamma_dInit, Real gamma_qInit);

	void withControl(Bool controlOn) { mWithControl = controlOn; };

	// #### Mathematical Matrix Transforms ####
	Matrix getParkTransformMatrixPowerInvariant(Real theta);
	Matrix parkTransformPowerInvariant(Real theta, const Matrix &fabc);
	Matrix getInverseParkTransformMatrixPowerInvariant(Real theta);
	Matrix inverseParkTransformPowerInvariant(Real theta, const Matrix &fdq);

	// #### MNA section ####
	void mnaParentInitialize(Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) override;
	void mnaCompUpdateCurrent(const Matrix& leftVector) override;
	void mnaCompUpdateVoltage(const Matrix& leftVector) override;
	void mnaParentPreStep(Real time, Int timeStepCount) override;
	void mnaParentPostStep(Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector) override;
	void mnaParentAddPreStepDependencies(AttributeBase::List &prevStepDependencies,
	                                     AttributeBase::List &attributeDependencies,
	                                     AttributeBase::List &modifiedAttributes) override;
	void mnaParentAddPostStepDependencies(AttributeBase::List &prevStepDependencies,
	                                      AttributeBase::List &attributeDependencies,
	                                      AttributeBase::List &modifiedAttributes,
	                                      Attribute<Matrix>::Ptr &leftVector) override;

	// #### Control section ####
	void controlPreStep(Real time, Int timeStepCount);
	void controlStep(Real time, Int timeStepCount);
	void addControlPreStepDependencies(AttributeBase::List &prevStepDependencies,
	                                   AttributeBase::List &attributeDependencies,
	                                   AttributeBase::List &modifiedAttributes);
	void addControlStepDependencies(AttributeBase::List &prevStepDependencies,
	                                AttributeBase::List &attributeDependencies,
	                                AttributeBase::List &modifiedAttributes);

	class ControlPreStep : public CPS::Task {
	public:
		ControlPreStep(VSIVoltageControlDQ& comp) :
			Task(**comp.mName + ".ControlPreStep"), mComp(comp) {
			mComp.addControlPreStepDependencies(mPrevStepDependencies, mAttributeDependencies, mModifiedAttributes);
		}
		void execute(Real time, Int timeStepCount) override { mComp.controlPreStep(time, timeStepCount); }
	private:
		VSIVoltageControlDQ& mComp;
	};

	class ControlStep : public CPS::Task {
	public:
		ControlStep(VSIVoltageControlDQ& comp) :
			Task(**comp.mName + ".ControlStep"), mComp(comp) {
			mComp.addControlStepDependencies(mPrevStepDependencies, mAttributeDependencies, mModifiedAttributes);
		}
		void execute(Real time, Int timeStepCount) override { mComp.controlStep(time, timeStepCount); }
	private:
		VSIVoltageControlDQ& mComp;
	};
};

} // namespace Ph3
} // namespace EMT
} // namespace CPS
