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
#include <dpsim-models/Solver/DAEInterface.h>
#include <dpsim-models/EMT/EMT_Ph3_VoltageSource.h>
#include <dpsim-models/Solver/DAEInterface.h>

namespace CPS {
	namespace EMT {
		namespace Ph3 {
			/// \brief Network injection model
			///
			/// This model represents network injections by an ideal voltage source.
			class NetworkInjection :
				public CompositePowerComp<Real>,
				public DAEInterface,
				public SharedFactory<NetworkInjection> {
			private:
				// ### Electrical Subcomponents ###
				/// Voltage source
				std::shared_ptr<EMT::Ph3::VoltageSource> mSubVoltageSource;

				// #### solver ####
				/// Vector to collect subcomponent right vector stamps
				std::vector<const Matrix*> mRightVectorStamps;
			public:
				const CPS::Attribute<MatrixComp>::Ptr mVoltageRef;
				const CPS::Attribute<Real>::Ptr mSrcFreq;

				/// Defines UID, name, component parameters and logging level
				NetworkInjection(String uid, String name, Logger::Level loglevel = Logger::Level::off);
				/// Defines UID, name, component parameters and logging level
				NetworkInjection(String name, Logger::Level logLevel = Logger::Level::off)
					: NetworkInjection(name, name, logLevel) { }
				/// Defines name, component parameters and logging level
				NetworkInjection(String name,
					Complex voltage, Logger::Level logLevel = Logger::Level::off);
				///
				SimPowerComp<Real>::Ptr clone(String name);

				// #### General ####
				/// Initializes component from power flow data
				void initializeFromNodesAndTerminals(Real frequency);
				/// Setter for reference voltage parameters
				void setParameters(MatrixComp voltageRef, Real srcFreq = 50.0);
				/// Setter for reference signal of type frequency ramp
				void setParameters(MatrixComp voltageRef, Real freqStart, Real rocof, Real timeStart, Real duration, bool smoothRamp = true);
				/// Setter for reference signal of type cosine frequency modulation
				void setParameters(MatrixComp voltageRef, Real modulationFrequency, Real modulationAmplitude, Real baseFrequency = 50.0, bool zigzag = false);

				// #### MNA Section ####
				/// Stamps right side (source) vector
				void mnaParentApplyRightSideVectorStamp(Matrix& rightVector) override;
				/// Returns current through the component
				void mnaCompUpdateCurrent(const Matrix& leftVector) override;
				/// Updates voltage across component
				void mnaCompUpdateVoltage(const Matrix& leftVector) override;
				/// MNA pre step operations
				void mnaParentPreStep(Real time, Int timeStepCount) override;
				/// MNA post step operations
				void mnaParentPostStep(Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector) override;
				/// Add MNA pre step dependencies
				void mnaParentAddPreStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes) override;
				/// Add MNA post step dependencies
				void mnaParentAddPostStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes, Attribute<Matrix>::Ptr &leftVector) override;

				// #### DAE Section ####
				
				/// Derivative of the current
				MatrixVar<Real> mIntfDerCurrent;
				void setInitialComplexIntfCurrent(Complex initCurrent);
				///
				void daeInitialize(double time, double state[], double dstate_dt[],
					double absoluteTolerances[], double stateVarTypes[], int& offset) override;
				/// Residual function for DAE Solver
				void daeResidual(double time, const double state[], const double dstate_dt[], double resid[], 
					 std::vector<int>& off) override;
				/// Calculation of jacobian
				void daeJacobian(double current_time, const double state[], const double dstate_dt[], 
					SUNMatrix jacobian, double cj, std::vector<int>& off) override;
				///
				void daePostStep(double Nexttime, const double state[], const double dstate_dt[], int& offset) override;
				///
				int getNumberOfStateVariables() override {return 3;}
			};
		}
	}
}
