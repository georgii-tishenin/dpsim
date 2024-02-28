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

	///FIXME: f_theta instead of theta
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
		//Internal Voltages are v1-v0: Positive for Terminal 1, negative for terminal 0
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
	//Phase current output equations do not contain constant history terms

	//Math::setVectorElement(rightVector, matrixNodeIndex(0, PhaseType::A), 0.);
	//Math::setVectorElement(rightVector, matrixNodeIndex(0, PhaseType::B), 0.);
	//Math::setVectorElement(rightVector, matrixNodeIndex(0, PhaseType::C), 0.);

	//Math::setVectorElement(rightVector, matrixNodeIndex(0, PhaseType::A), 0.);
	//Math::setVectorElement(rightVector, matrixNodeIndex(0, PhaseType::B), 0.);
	//Math::setVectorElement(rightVector, matrixNodeIndex(0, PhaseType::C), 0.);

	//Equation for theta does contain constant history terms:

	double linear_theta_hist = -C_wb*(P_mech_old-(V_d_old*V_d_old/X_Q)+(V_d_old*E_d_old/X_Q)-(V_q_old*E_q_old/X_D)+(V_q*V_q/X_D))-(0.5*timeStep*omega_old)-(0.5*timeStep*omega_base(omega_old-2.))-omega_old-1.;
	Math::setVectorElement(rightVector, mVirtualNodes[0]->matrixNodeIndex(PhaseType::Single), linear_theta_hist);
}

void CPS::EMT::Ph3::SynchronGenerator4OrderSSN::mnaCompAddPreStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes) {
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
	Jacobian(0,0)=(cos(theta)*sin(theta)*(1.-C_dd)/X_D)-(cos(theta)*sin(theta)*(1.-C_qq)/X_Q);
	Jacobian(0,1)=(cos(theta)*sin(theta-(2.*M_PI/3.))*(1.-C_dd)/X_D)-(cos(theta-(2.*M_PI/3.))*sin(theta)*(1.-C_qq)/X_Q);
	Jacobian(0,2)=(cos(theta)*sin(theta+(2.*M_PI/3.))*(1.-C_dd)/X_D)-(cos(theta+(2.*M_PI/3.))*sin(theta)*(1.-C_qq)/X_Q);
	Jacobian(0,3)=(cos(theta)/X_D)*(V_a*cos(theta)+V_b*cos(theta-(2.*M_PI/3.))+V_c*cos(theta+(2.*M_PI/3.)))*(1.-C_dd)-(sin(theta)/X_D)*(E_q-V_q)+(sin(theta)/X_Q)*(V_a*sin(theta)+V_b*sin(theta-(2.*M_PI/3.))+V_c*sin(theta+(2.*M_PI/3.)))*(1.-C_qq)+(cos(theta)/X_Q)*(E_d-V_d);
	Jacobian(1,0)=(cos(theta-(2.*M_PI/3.))*sin(theta)*(1.-C_dd)/X_D)-(cos(theta)*sin(theta-(2.*M_PI/3.))*(1.-C_qq)/X_Q);
	Jacobian(1,1)=(cos(theta-(2.*M_PI/3.))*sin(theta-(2.*M_PI/3.))*(1.-C_dd)/X_D)-(cos(theta-(2.*M_PI/3.))*sin(theta-(2.*M_PI/3.))*(1.-C_qq)/X_Q);
	Jacobian(1,2)=(cos(theta-(2.*M_PI/3.))*sin(theta+(2.*M_PI/3.))*(1.-C_dd)/X_D)-(cos(theta+(2.*M_PI/3.))*sin(theta-(2.*M_PI/3.))*(1.-C_qq)/X_Q);
	Jacobian(1,3)=(cos(theta-(2.*M_PI/3.))/X_D)*(V_a*cos(theta)+V_b*cos(theta-(2.*M_PI/3.))+V_c*cos(theta+(2.*M_PI/3.)))*(1.-C_dd)-(sin(theta-(2.*M_PI/3.))/X_D)*(E_q-V_q)+(sin(theta-(2.*M_PI/3.))/X_Q)*(V_a*sin(theta)+V_b*sin(theta-(2.*M_PI/3.))+V_c*sin(theta+(2.*M_PI/3.)))*(1.-C_qq)+(cos(theta-(2.*M_PI/3.))/X_Q)*(E_d-V_d);
	Jacobian(2,0)=(cos(theta+(2.*M_PI/3.))*sin(theta)*(1.-C_dd)/X_D)-(cos(theta)*sin(theta+(2.*M_PI/3.))*(1.-C_qq)/X_Q);
	Jacobian(2,1)=(cos(theta+(2.*M_PI/3.))*sin(theta-(2.*M_PI/3.))*(1.-C_dd)/X_D)-(cos(theta-(2.*M_PI/3.))*sin(theta+(2.*M_PI/3.))*(1.-C_qq)/X_Q);
	Jacobian(2,2)=(cos(theta+(2.*M_PI/3.))*sin(theta+(2.*M_PI/3.))*(1.-C_dd)/X_D)-(cos(theta+(2.*M_PI/3.))*sin(theta+(2.*M_PI/3.))*(1.-C_qq)/X_Q);
	Jacobian(2,3)=(cos(theta+(2.*M_PI/3.))/X_D)*(V_a*cos(theta)+V_b*cos(theta-(2.*M_PI/3.))+V_c*cos(theta+(2.*M_PI/3.)))*(1.-C_dd)-(sin(theta+(2.*M_PI/3.))/X_D)*(E_q-V_q)+(sin(theta+(2.*M_PI/3.))/X_Q)*(V_a*sin(theta)+V_b*sin(theta-(2.*M_PI/3.))+V_c*sin(theta+(2.*M_PI/3.)))*(1.-C_qq)+(cos(theta+(2.*M_PI/3.))/X_Q)*(E_d-V_d);
	Jacobian(3,0)=C_wbq*V_d*cos(theta)*(1.-C_qq)-0.5*C_wbq*cos(theta)*(C_qq*V_d_old+C_0qq*E_d_old)-C_wbd*V_q*sin(theta)*(C_dd-1.)-0.5*C_wbd*sin(theta)*(C_d*(E_f+E_f_old)+C_0dd*E_q_old+C_dd*V_q_old);
	Jacobian(3,1)=C_wbq*V_d*cos(theta-(2.*M_PI/3.))*(1.-C_qq)-0.5*C_wbq*cos(theta-(2.*M_PI/3.))*(C_qq*V_d_old+C_0qq*E_d_old)-C_wbd*V_q*sin(theta-(2.*M_PI/3.))*(C_dd-1.)-0.5*C_wbd*sin(theta-(2.*M_PI/3.))*(C_d*(E_f+E_f_old)+C_0dd*E_q_old+C_dd*V_q_old);
	Jacobian(3,2)=C_wbq*V_d*cos(theta+(2.*M_PI/3.))*(1.-C_qq)-0.5*C_wbq*cos(theta+(2.*M_PI/3.))*(C_qq*V_d_old+C_0qq*E_d_old)-C_wbd*V_q*sin(theta+(2.*M_PI/3.))*(C_dd-1.)-0.5*C_wbd*sin(theta+(2.*M_PI/3.))*(C_d*(E_f+E_f_old)+C_0dd*E_q_old+C_dd*V_q_old);
	Jacobian(3,3)=C_wbq*V_d*(V_a*sin(theta)+V_b*sin(theta-(2.*M_PI/3.))+V_c*sin(theta+(2.*M_PI/3.)))*(1.-C_qq)-0.5*C_wbq*(V_a*sin(theta)+V_b*sin(theta-(2.*M_PI/3.))+V_c*sin(theta+(2.*M_PI/3.)))*(C_qq*V_d_old+C_0qq*E_d_old)+C_wbd*V_q*(V_a*cos(theta)+V_b*cos(theta-(2.*M_PI/3.))+V_c*cos(theta+(2.*M_PI/3.)))*(C_dd-1.)+0.5*C_wbd*(V_a*cos(theta)+V_b*cos(theta-(2.*M_PI/3.))+V_c*cos(theta+(2.*M_PI/3.)))*(C_d*(E_f+E_f_old)+C_0dd*E_q_old+C_dd*V_q_old)-1.;
}

void CPS::EMT::Ph3::SynchronGenerator4OrderSSN::updateStates()
{
	E_d		= C_qq*(V_d+V_d_old)+C_0qq*E_d_old;
	E_q		= C_dd*(V_q+V_q_old)+C_d*(E_f+E_f_old)+C_0dd*E_q_old;
	omega	= C_h*(P_mech-(V_d*V_d/X_Q)+(V_d*E_d/X_Q)-(V_q*E_q/X_d)+(V_q*V_q/X_D))+C_h*(P_mech_old-(V_d_old*V_d_old/X_Q)+(V_d_old*E_d_old/X_Q)-(V_q_old*E_q_old/X_D)+(V_q_old*V_q_old/X_D))+omega_old;
}
