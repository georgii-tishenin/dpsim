#include "EMT_Ph3_SynchronGenerator4OrderSSN.h"

CPS::EMT::Ph3::SynchronGenerator4OrderSSN::SynchronGenerator4OrderSSN(String uid, String name, Logger::Level logLevel)
    	: MNASimPowerComp<Real>(uid, name, false, false, logLevel)
{
    mPhaseType = PhaseType::ABC;
    setTerminalNumber(2);
    setVirtualNodeNumber(1);
    **mIntfCurrent = Matrix::Zero(3,1);
    **mIntfVoltage = Matrix::Zero(3,1);
}

void CPS::EMT::Ph3::SynchronGenerator4OrderSSN::calculateNonlinearFunctionResult()
{

    double Vd = mIntfVoltage(0,0)*cos(theta)+mIntfVoltage(1,0)*cos(theta-(2*M_PI/3))+mIntfVoltage(2,0)*cos(theta+(2*M_PI/3));
    double Vq = -(mIntfVoltage(0,0)*sin(theta)+mIntfVoltage(1,0)*sin(theta-(2*M_PI/3))+mIntfVoltage(2,0)*sin(theta+(2*M_PI/3)));

    **intfCurrent(0,0) = ((Eq-Vq)/Xd)*cos(theta) - ((Vd-Ed)/Xq)*sin(theta);
    **intfCurrent(1,0) = ((Eq-Vq)/Xd)*cos(theta - (2*M_PI/3)) - ((Vd-Ed)/Xq)*sin(theta - (2*M_PI/3));
    **intfCurrent(2,0) = ((Eq-Vq)/Xd)*cos(theta + (2*M_PI/3)) - ((Vd-Ed)/Xq)*sin(theta + (2*M_PI/3));

	theta = ((mTimeStep*mTimeStep*omega_base)/(8.*H))*P_mech - ((mTimeStep*mTimeStep*omega_base)/(8.*H))*((Vd*Vd)/Xq_dash)
			+((mTimeStep*mTimeStep*omega_base)/(8.*H))*((Vd)/Xq_dash)*(
			((mTimeStep*(Xq-Xq_dash))/(2*Tq0_dash*Xq_dash+mTimeStep*Xq))*Vd+((mTimeStep*(Xq-Xq_dash))/(2.*Tq0_dash*Xq_dash+mTimeStep*Xq))*Vd_old
			+((2*Tq0_dash*Xq_dash-Xq*mTimeStep)/(2.*Tq0_dash*Xq_dash+mTimeStep*Xq))*Ed_dash_old)
			-((mTimeStep*mTimeStep*omega_base)/(8.*H))*(Vq/Xd_dash)*(
			((mTimeStep*(Xd-Xd_dash))/(2*Td0_dash*Xd_dash+mTimeStep*Xd))*Vq+((mTimeStep*Xd_dash)/(2.*Td0_dash*Xd_dash+mTimeStep*Xd))*Ef
			+((2*Td0_dash*Xd_dash-mTimeStep*Xd)/(2.*Td0_dash*Xd_dash+mTimeStep*Xd))*Eq_dash_old
			+((mTimeStep*(Xd-Xd_dash))/(2.*Td0_dash*Xd_dash+mTimeStep*Xd))*Vq_old
			+((mTimeStep*Xd_dash)/(2.*Td0_dash*Xd_dash+mTimeStep*Xd))*Ef_old)
			+((mTimeStep*mTimeStep*omega_base)/(8.*H))*((Vq*Vq)/Xd_dash)
			+((mTimeStep*mTimeStep*omega_base)/(8.*H))*(P_mech_old-((Vd_old*Vd_old)/Xq_dash)+((Vd_old*Ed_dash_old)/Xq_dash)-((Vq_old*Eq_dash_old)/Xd_dash)+(Vq_old*Vq_old)/Xd_dash)
			+((mTimeStep*omega_base)/2.)*omega_old+((mTimeStep*omega_base*(omega_old-2.))/2.)+theta_old;

    if (terminalNotGrounded(0)) {
		Math::setVectorElement(mNonlinearFunctionStamp, matrixNodeIndex(0, 0), -(**intfCurrent)(0,0));
		Math::setVectorElement(mNonlinearFunctionStamp, matrixNodeIndex(0, 1), -(**intfCurrent)(1,0));
		Math::setVectorElement(mNonlinearFunctionStamp, matrixNodeIndex(0, 2), -(**intfCurrent)(2,0));
	}
	if (terminalNotGrounded(1)) {
		Math::setVectorElement(mNonlinearFunctionStamp, matrixNodeIndex(1, 0), (**intfCurrent)(0,0));
		Math::setVectorElement(mNonlinearFunctionStamp, matrixNodeIndex(1, 1), (**intfCurrent)(1,0));
		Math::setVectorElement(mNonlinearFunctionStamp, matrixNodeIndex(1, 2), (**intfCurrent)(2,0));
	}
	Math::setVectorElement(mNonlinearFunctionStamp, mVirtualNodes[0]->matrixNodeIndex(PhaseType::Single), theta);
}

void CPS::EMT::Ph3::SynchronGenerator4OrderSSN::mnaCompApplySystemMatrixStamp(Matrix& systemMatrix) {
	if (terminalNotGrounded(0)) {
		// set upper left block, 3x3 entries
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 0), matrixNodeIndex(0, 0), Jacobian(0, 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 0), matrixNodeIndex(0, 1), Jacobian(0, 1));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 0), matrixNodeIndex(0, 2), Jacobian(0, 2));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 0), mVirtualNodes[0]->matrixNodeIndex(PhaseType::Single), Jacobian(0, 3));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 1), matrixNodeIndex(0, 0), Jacobian(1, 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 1), matrixNodeIndex(0, 1), Jacobian(1, 1));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 1), matrixNodeIndex(0, 2), Jacobian(1, 2));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 1), mVirtualNodes[0]->matrixNodeIndex(PhaseType::Single), Jacobian(1, 3));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 2), matrixNodeIndex(0, 0), Jacobian(2, 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 2), matrixNodeIndex(0, 1), Jacobian(2, 1));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 2), matrixNodeIndex(0, 2), Jacobian(2, 2));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 2), mVirtualNodes[0]->matrixNodeIndex(PhaseType::Single), -Jacobian(2, 3));
	}
	if (terminalNotGrounded(1)) {
		// set bottom right block, 3x3 entries
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 0), matrixNodeIndex(1, 0), Jacobian(0, 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 0), matrixNodeIndex(1, 1), Jacobian(0, 1));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 0), matrixNodeIndex(1, 2), Jacobian(0, 2));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 0), mVirtualNodes[0]->matrixNodeIndex(PhaseType::Single), Jacobian(0, 3));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 1), matrixNodeIndex(1, 0), Jacobian(1, 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 1), matrixNodeIndex(1, 1), Jacobian(1, 1));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 1), matrixNodeIndex(1, 2), Jacobian(1, 2));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 1), mVirtualNodes[0]->matrixNodeIndex(PhaseType::Single), Jacobian(1, 3));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 2), matrixNodeIndex(1, 0), Jacobian(2, 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 2), matrixNodeIndex(1, 1), Jacobian(2, 1));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 2), matrixNodeIndex(1, 2), Jacobian(2, 2));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 2), mVirtualNodes[0]->matrixNodeIndex(PhaseType::Single), Jacobian(2, 3));
	}
	// Set off diagonal blocks, 2x3x3 entries
	if (terminalNotGrounded(0) && terminalNotGrounded(1)) {
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 0), matrixNodeIndex(1, 0), -Jacobian(0, 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 0), matrixNodeIndex(1, 1), -Jacobian(0, 1));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 0), matrixNodeIndex(1, 2), -Jacobian(0, 2));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 1), matrixNodeIndex(1, 0), -Jacobian(1, 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 1), matrixNodeIndex(1, 1), -Jacobian(1, 1));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 1), matrixNodeIndex(1, 2), -Jacobian(1, 2));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 2), matrixNodeIndex(1, 0), -Jacobian(2, 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 2), matrixNodeIndex(1, 1), -Jacobian(2, 1));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 2), matrixNodeIndex(1, 2), -Jacobian(2, 2));


		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 0), matrixNodeIndex(0, 0), -Jacobian(0, 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 0), matrixNodeIndex(0, 1), -Jacobian(0, 1));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 0), matrixNodeIndex(0, 2), -Jacobian(0, 2));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 1), matrixNodeIndex(0, 0), -Jacobian(1, 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 1), matrixNodeIndex(0, 1), -Jacobian(1, 1));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 1), matrixNodeIndex(0, 2), -Jacobian(1, 2));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 2), matrixNodeIndex(0, 0), -Jacobian(2, 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 2), matrixNodeIndex(0, 1), -Jacobian(2, 1));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 2), matrixNodeIndex(0, 2), -Jacobian(2, 2));
	}
		//Internal Voltages are v1-v0: Positive for Terminal 1, negative for terminal 0 (Thus also the signs for Jacobian(2, 3) above).
        //Theta is not a difference: We only have one "virtual terminal" for it.

		Math::addToMatrixElement(systemMatrix, mVirtualNodes[0]->matrixNodeIndex(PhaseType::Single), matrixNodeIndex(0, 0) ,-Jacobian(3, 0));
		Math::addToMatrixElement(systemMatrix, mVirtualNodes[0]->matrixNodeIndex(PhaseType::Single), matrixNodeIndex(0, 1) ,-Jacobian(3, 1));
		Math::addToMatrixElement(systemMatrix, mVirtualNodes[0]->matrixNodeIndex(PhaseType::Single), matrixNodeIndex(0, 2) ,-Jacobian(3, 2));

		Math::addToMatrixElement(systemMatrix, mVirtualNodes[0]->matrixNodeIndex(PhaseType::Single), matrixNodeIndex(1, 0) ,Jacobian(3, 0));
		Math::addToMatrixElement(systemMatrix, mVirtualNodes[0]->matrixNodeIndex(PhaseType::Single), matrixNodeIndex(1, 1) ,Jacobian(3, 1));
		Math::addToMatrixElement(systemMatrix, mVirtualNodes[0]->matrixNodeIndex(PhaseType::Single), matrixNodeIndex(1, 2) ,Jacobian(3, 2));

		Math::addToMatrixElement(systemMatrix, mVirtualNodes[0]->matrixNodeIndex(PhaseType::Single), mVirtualNodes[0]->matrixNodeIndex(PhaseType::Single), Jacobian(3, 3));

		//Math::addToMatrixElement(systemMatrix, mVirtualNodes[0]->matrixNodeIndex(PhaseType::B), mVirtualNodes[0]->matrixNodeIndex(PhaseType::B), 1.);
		//Math::addToMatrixElement(systemMatrix, mVirtualNodes[0]->matrixNodeIndex(PhaseType::C), mVirtualNodes[0]->matrixNodeIndex(PhaseType::C), 1.);
}

void CPS::EMT::Ph3::SynchronGenerator4OrderSSN::mnaCompInitialize(Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
		updateMatrixNodeIndices();


}

void CPS::EMT::Ph3::SynchronGenerator4OrderSSN::mnaCompApplyRightSideVectorStamp(Matrix& rightVector) {

}

void CPS::EMT::Ph3::SynchronGenerator4OrderSSN::mnaCompAddPreStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes) {
	// actually depends on C, but then we'd have to modify the system matrix anyway
	modifiedAttributes.push_back(mRightVector);
	prevStepDependencies.push_back(mIntfCurrent);
	prevStepDependencies.push_back(mIntfVoltage);
}


void CPS::EMT::Ph3::SynchronGenerator4OrderSSN::mnaCompPreStep(Real time, Int timeStepCount) {
	mnaCompApplyRightSideVectorStamp(**mRightVector);
}

void CPS::EMT::Ph3::SynchronGenerator4OrderSSN::mnaCompAddPostStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes, Attribute<Matrix>::Ptr &leftVector) {
	attributeDependencies.push_back(leftVector);
	modifiedAttributes.push_back(mIntfVoltage);
	modifiedAttributes.push_back(mIntfCurrent);
}

void CPS::EMT::Ph3::SynchronGenerator4OrderSSN::mnaCompPostStep(Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector) {
    updateStates();
	mnaCompUpdateVoltage(**leftVector);
	mnaCompUpdateCurrent(**leftVector);
}

void CPS::EMT::Ph3::SynchronGenerator4OrderSSN::mnaCompUpdateVoltage(const Matrix& leftVector) {
	// v1 - v0
	(**mIntfVoltage) = Matrix::Zero(3,1);
	if (terminalNotGrounded(1))
		(**mIntfVoltage)(0,0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(1,0));
		(**mIntfVoltage)(1,0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(1,1));
		(**mIntfVoltage)(2,0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(1,2));
	if (terminalNotGrounded(0))
		(**mIntfVoltage)(0,0) = (**mIntfVoltage)(0,0) - Math::realFromVectorElement(leftVector, matrixNodeIndex(0,0));
        (**mIntfVoltage)(1,0) = (**mIntfVoltage)(1,0) - Math::realFromVectorElement(leftVector, matrixNodeIndex(0,1));
        (**mIntfVoltage)(2,0) = (**mIntfVoltage)(2,0) - Math::realFromVectorElement(leftVector, matrixNodeIndex(0,2));
}

void CPS::EMT::Ph3::SynchronGenerator4OrderSSN::mnaCompUpdateCurrent(const Matrix& leftVector) {
	calculateNonlinearFunctionResult();
}


void CPS::EMT::Ph3::SynchronGenerator4OrderSSN::iterationUpdate(const Matrix& leftVector)
{
    mnaCompUpdateVoltage(leftVector);
	calculateNonlinearFunctionResult();
	updateJacobian();
}

void CPS::EMT::Ph3::SynchronGenerator4OrderSSN::updateJacobian()
{
}

void CPS::EMT::Ph3::SynchronGenerator4OrderSSN::updateStates()
{
}
