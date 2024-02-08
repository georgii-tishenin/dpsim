/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#pragma once

#include <dpsim-models/Base/Base_ReducedOrderSynchronGenerator.h>
#include <dpsim-models/Solver/MNANonlinearVariableCompInterface.h>

namespace CPS{
namespace EMT{
namespace Ph3{
class SynchronGenerator4OrderSSN:  
    public Base::ReducedOrderSynchronGenerator<Real>,
    public MNANonlinearVariableCompInterface {

public:
	/// Defines UID, name, component parameters and logging level
	SynchronGenerator4OrderSSN(String uid, String name, Logger::Level logLevel = Logger::Level::off);
	/// Defines name, component parameters and logging level
	SynchronGenerator4OrderSSN(String name, Logger::Level logLevel = Logger::Level::off)
		: SynchronGenerator4OrderSSN(name, name, logLevel) { }

    virtual void calculateNonlinearFunctionResult() override;
    virtual void mnaCompApplySystemMatrixStamp(Matrix& systemMatrix) override;
    virtual void mnaCompPostStep(Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector) override;
    virtual void mnaCompUpdateVoltage(const Matrix& leftVector) override;
    virtual void mnaCompUpdateCurrent(const Matrix& leftVector) override;
    virtual void mnaCompAddPostStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes, Attribute<Matrix>::Ptr &leftVector) override;
    virtual void mnaCompPreStep(Real time, Int timeStepCount) override;
    virtual void mnaCompAddPreStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes) override;
    virtual void mnaCompApplyRightSideVectorStamp(Matrix& rightVector) override;
    virtual void mnaCompInitialize(Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) override;
    void updateJacobian();
    void updateStates();
    virtual void iterationUpdate(const Matrix& leftVector) override;

protected:
    Matrix Jacobian = Matrix::Zero(4,4);

    //inputs
    double P_mech;
    double P_mech_old;
    double Ef;
    double Ef_old;

    //states
    double theta;
    double theta_old;
    double Ed;
    double Ed_old;
    double Eq;
    double Eq_old;
    double omega;
    double omega_old;

private:

}
}
}
}

