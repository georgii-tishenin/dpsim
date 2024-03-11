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
    public MNANonlinearVariableCompInterface,
    public SharedFactory<SynchronGenerator4OrderSSN> {

public:
	/// Defines UID, name, component parameters and logging level
	SynchronGenerator4OrderSSN(String uid, String name, Logger::Level logLevel = Logger::Level::off);
	/// Defines name, component parameters and logging level
	SynchronGenerator4OrderSSN(String name, Logger::Level logLevel = Logger::Level::off)
		: SynchronGenerator4OrderSSN(name, name, logLevel) { }

    virtual void calculateNonlinearFunctionResult() override;
    virtual void mnaCompApplySystemMatrixStamp(SparseMatrixRow& systemMatrix) override;
    virtual void mnaCompPostStep(const Matrix &leftVector) override;
    virtual void mnaCompUpdateVoltage(const Matrix& leftVector) override;
    virtual void mnaCompUpdateCurrent(const Matrix& leftVector) override;
    virtual void mnaCompAddPostStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes, Attribute<Matrix>::Ptr &leftVector) override;
    virtual void mnaCompPreStep(Real time, Int timeStepCount) override;
    virtual void mnaCompAddPreStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes) override;
    virtual void mnaCompApplyRightSideVectorStamp(Matrix& rightVector) override;
    virtual void mnaCompInitialize(Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) override;
    void updateJacobian();
    void updateCurrentStates();
    void updateOldStates();
    void updateImplicitStates(const Matrix& leftVector);
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

    //helpers
    double Vd = 0;
    double Vq = 0;
    double Vd_old = 0;
    double Vq_old = 0;

private:
    //constants
    const double C_d = (mTimeStep*mLd_t)/(2.*mTd0_t*mLd_t+mTimeStep*mLd);
    const double C_dd = (mTimeStep*(mLd-mLd_t))/(2.*mTd0_t*mLd_t+mTimeStep*mLd);
    const double C_0dd = (2.*mTd0_t*mLd_t-mTimeStep*mLd)/(2.*mTd0_t*mLd_t+mTimeStep*mLd);
    const double C_qq = (mTimeStep*(mLq-mLq_t))/(2.*mTq0_t*mLq_t+mTimeStep*mLq);
    const double C_0qq = (2.*mTq0_t*mLq_t-mTimeStep*mLq)/(2.*mTq0_t*mLq_t+mTimeStep*mLq);
    const double C_wbq = (mTimeStep*mTimeStep*mBase_OmElec)/(4.*mH*mLq_t);
    const double C_wbd = (mTimeStep*mTimeStep*mBase_OmElec)/(4.*mH*mLd_t);
    const double C_wb = (mTimeStep*mTimeStep*mBase_OmElec)/(8.*mH);
    const double C_h = (mTimeStep)/(4.*mH);
};
}
}
}

