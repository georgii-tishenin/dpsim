
#include <dpsim-models/EMT/EMT_Ph1_ParkTransformer.h>
#include <dpsim-models/EMT/EMT_Ph1_InertiaMoment.h>


using namespace CPS;


EMT::Ph1::ParkTransformer::ParkTransformer(String uid, String name,
    Logger::Level logLevel)
: MNASimPowerComp<Real>(uid, name, true, true, logLevel), mOmega(std::make_shared<Real>()),
  mTheta_init(std::make_shared<Real>()), mTheta(std::make_shared<Real>()), mTheta_atr(mAttributes->create<Real>("theta")){

setVirtualNodeNumber(3);

setTerminalNumber(6);

//mTheta_atr = mAttributes->createDynamic<Real>("theta");

SPDLOG_LOGGER_INFO(mSLog, "Create {} {}", this->type(), name);
**mIntfVoltage = Matrix::Zero(3, 1);
**mIntfCurrent = Matrix::Zero(3, 1);
}


void EMT::Ph1::ParkTransformer::setParameters(Real omega, Real theta_initial) {

    *mOmega = omega;

    if(theta_initial < 0){
        theta_initial = theta_initial + 2*M_PI;
    }

    *mTheta_init = theta_initial;

   
    SPDLOG_LOGGER_INFO(mSLog, "The rotational frequency={} [ ] ",
                        std::abs(omega));
    SPDLOG_LOGGER_INFO(mSLog, "The initial electrical angle={} [ ] ",
                        std::abs(theta_initial));
   
     mParametersSet = true;
     }





     void EMT::Ph1::ParkTransformer::initializeFromNodesAndTerminals(Real frequency) {

 
        // Set initial voltage of virtual nodes
      mVirtualNodes[0]->setInitialVoltage(initialSingleVoltage(1));
      mVirtualNodes[1]->setInitialVoltage(initialSingleVoltage(2));
      mVirtualNodes[2]->setInitialVoltage(initialSingleVoltage(3));


      *mTheta = *mTheta_init;
    
    
      // Log initialization results
      SPDLOG_LOGGER_INFO(
          mSLog,
          "--- Initialization ---\n"
          "Terminal 0 voltage: {:s}\n"
          "Terminal 1 voltage: {:s}\n"
          "Terminal 2 voltage: {:s}\n"
          "Terminal 3 voltage: {:s}\n"
          "Terminal 4 voltage: {:s}\n"
          "Terminal 5 voltage: {:s}\n"
          "Virtual Node 1 voltage: {:s}\n"
          "Virtual Node 2 voltage: {:s}\n"
          "Virtual Node 3 voltage: {:s}\n"
          "--- Initialization finished ---",
          Logger::phasorToString(initialSingleVoltage(0)),
          Logger::phasorToString(initialSingleVoltage(1)),
          Logger::phasorToString(initialSingleVoltage(2)),
          Logger::phasorToString(initialSingleVoltage(3)),
          Logger::phasorToString(initialSingleVoltage(4)),
          Logger::phasorToString(initialSingleVoltage(5)),
          Logger::phasorToString(mVirtualNodes[0]->initialSingleVoltage()),
          Logger::phasorToString(mVirtualNodes[1]->initialSingleVoltage()),
          Logger::phasorToString(mVirtualNodes[2]->initialSingleVoltage()));
    }
    
    
    
void EMT::Ph1::ParkTransformer::mnaCompInitialize(Real omega, Real timeStep,
                                               Attribute<Matrix>::Ptr leftVector) {
      updateMatrixNodeIndices();
      //**mRightVector = Matrix::Zero(0, 0);
}
    
void EMT::Ph1::ParkTransformer::mnaCompPreStep(Real time,
        Int timeStepCount) {

    if (mIsOmegaConstant == false) {

        *mTheta = *mTheta + ( mTimeStep / 2 ) * (*mOmega + mOmega_prev);
        *mTheta = std::fmod(*mTheta, 2*M_PI);
        mTheta_atr->set(*mTheta);
        mOmega_prev = *mOmega;

    } else {
    
    *mTheta = std::fmod(*mTheta_init + *mOmega * time, 2*M_PI);
    mTheta_atr->set(*mTheta);

  }
    mnaCompApplyRightSideVectorStamp(**mRightVector);
}





void EMT::Ph1::ParkTransformer::mnaCompApplySystemMatrixStamp(
    SparseMatrixRow &systemMatrix) {
  
  // Park transformation equations

    // Entries related to virtual node 1
  if (terminalNotGrounded(0)) {

    Math::setMatrixElement(systemMatrix, mVirtualNodes[0]->matrixNodeIndex(), matrixNodeIndex(0), 1); // virtual node 1 / node a
    Math::setMatrixElement(systemMatrix, matrixNodeIndex(0), mVirtualNodes[0]->matrixNodeIndex(), 1); // node a / virtual node 1

    Math::setMatrixElement(systemMatrix, mVirtualNodes[0]->matrixNodeIndex(), matrixNodeIndex(3), -sqrt(2.0/3.0)*cos(*mTheta)); //  virtual node 1 / node d
    Math::setMatrixElement(systemMatrix, matrixNodeIndex(3), mVirtualNodes[0]->matrixNodeIndex(), -sqrt(2.0/3.0)*cos(*mTheta)); // node d / virtual node 1

    Math::setMatrixElement(systemMatrix, mVirtualNodes[0]->matrixNodeIndex(), matrixNodeIndex(4), +sqrt(2.0/3.0)*sin(*mTheta)); // virtual node 1 / node q
    Math::setMatrixElement(systemMatrix, matrixNodeIndex(4), mVirtualNodes[0]->matrixNodeIndex(), +sqrt(2.0/3.0)*sin(*mTheta)); // node q / virtual node 1


    Math::setMatrixElement(systemMatrix, mVirtualNodes[0]->matrixNodeIndex(), matrixNodeIndex(5), -sqrt(1.0/3.0));  // virtual node 1 / node 0    
    Math::setMatrixElement(systemMatrix, matrixNodeIndex(5), mVirtualNodes[0]->matrixNodeIndex(), -sqrt(1.0/3.0));  // node 0 / virtual node 1
  }

    // Entries related to virtual node 2
  if (terminalNotGrounded(1)) {
    Math::setMatrixElement(systemMatrix, mVirtualNodes[1]->matrixNodeIndex(), matrixNodeIndex(1), 1);  // node b / virtual node 2
    Math::setMatrixElement(systemMatrix, matrixNodeIndex(1), mVirtualNodes[1]->matrixNodeIndex(), 1);  // virtual node 2 / node b

    Math::setMatrixElement(systemMatrix, mVirtualNodes[1]->matrixNodeIndex(), matrixNodeIndex(3), -sqrt(2.0/3.0)*cos(*mTheta - 2*M_PI/3));  // virtual node 2 / node d
    Math::setMatrixElement(systemMatrix, matrixNodeIndex(3), mVirtualNodes[1]->matrixNodeIndex(), -sqrt(2.0/3.0)*cos(*mTheta - 2*M_PI/3));  // node d / virtual node 2

    Math::setMatrixElement(systemMatrix, mVirtualNodes[1]->matrixNodeIndex(), matrixNodeIndex(4), +sqrt(2.0/3.0)*sin(*mTheta - 2*M_PI/3));  // virtual node 2 / node q
    Math::setMatrixElement(systemMatrix, matrixNodeIndex(4), mVirtualNodes[1]->matrixNodeIndex(), +sqrt(2.0/3.0)*sin(*mTheta - 2*M_PI/3));  // node q / virtual node 2

    Math::setMatrixElement(systemMatrix, mVirtualNodes[1]->matrixNodeIndex(), matrixNodeIndex(5), -sqrt(1.0/3.0));  // virtual node 2 / node 0
    Math::setMatrixElement(systemMatrix, matrixNodeIndex(5), mVirtualNodes[1]->matrixNodeIndex(), -sqrt(1.0/3.0));  // node 0 / virtual node 2
  }

    // Entries related to virtual node 3
  if (terminalNotGrounded(2)) {
    Math::setMatrixElement(systemMatrix, mVirtualNodes[2]->matrixNodeIndex(), matrixNodeIndex(2), 1);  // virtual node 3 / node c
    Math::setMatrixElement(systemMatrix, matrixNodeIndex(2), mVirtualNodes[2]->matrixNodeIndex(), 1);  // node c / virtual node 3

    Math::setMatrixElement(systemMatrix, mVirtualNodes[2]->matrixNodeIndex(), matrixNodeIndex(3), -sqrt(2.0/3.0)*cos(*mTheta + 2*M_PI/3));  // virtual node 3 / node d
    Math::setMatrixElement(systemMatrix, matrixNodeIndex(3), mVirtualNodes[2]->matrixNodeIndex(), -sqrt(2.0/3.0)*cos(*mTheta + 2*M_PI/3));  // node d / virtual node 3

    Math::setMatrixElement(systemMatrix, mVirtualNodes[2]->matrixNodeIndex(), matrixNodeIndex(4), +sqrt(2.0/3.0)*sin(*mTheta + 2*M_PI/3));  // virtual node 3 / node q
    Math::setMatrixElement(systemMatrix, matrixNodeIndex(4), mVirtualNodes[2]->matrixNodeIndex(), +sqrt(2.0/3.0)*sin(*mTheta + 2*M_PI/3));  // node q / virtual node 3

    Math::setMatrixElement(systemMatrix, mVirtualNodes[2]->matrixNodeIndex(), matrixNodeIndex(5), -sqrt(1.0/3.0));  // virtual node 3 / node 0
    Math::setMatrixElement(systemMatrix, matrixNodeIndex(5), mVirtualNodes[2]->matrixNodeIndex(), -sqrt(1.0/3.0));  // node 0 / virtual node 3
  }

}




void EMT::Ph1::ParkTransformer::mnaCompAddPostStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes,
    Attribute<Matrix>::Ptr &leftVector) {
  attributeDependencies.push_back(leftVector);
  modifiedAttributes.push_back(mIntfVoltage);
  modifiedAttributes.push_back(mIntfCurrent);
}



void EMT::Ph1::ParkTransformer::mnaCompPostStep(Real time, Int timeStepCount,
    Attribute<Matrix>::Ptr &leftVector) {
mnaCompUpdateVoltage(**leftVector);
mnaCompUpdateCurrent(**leftVector);
//*mTheta = std::fmod(*mTheta_init + *mOmega * time, 2*M_PI);
//mTheta_atr->set(*mTheta);

if(mIsOmegaConstant == false) {
  updateOmega();
}
}


void EMT::Ph1::ParkTransformer::mnaCompUpdateVoltage(const Matrix &leftVector) {


  /*

    (**mIntfVoltage)(0, 0) = 0;
  if (terminalNotGrounded(3))
    (**mIntfVoltage)(0, 0) =
        Math::realFromVectorElement(leftVector, matrixNodeIndex(3));


  (**mIntfVoltage)(1, 0) = 0;
  if (terminalNotGrounded(4))
    (**mIntfVoltage)(1, 0) =
        Math::realFromVectorElement(leftVector, matrixNodeIndex(4));


  (**mIntfVoltage)(2, 0) = 0;
  if (terminalNotGrounded(5))
    (**mIntfVoltage)(2, 0) =
        Math::realFromVectorElement(leftVector, matrixNodeIndex(5));

  */

  /*     
  (**mIntfVoltage)(0, 0) =
      Math::realFromVectorElement(leftVector, matrixNodeIndex(3));

  (**mIntfVoltage)(1, 0) =
      Math::realFromVectorElement(leftVector, matrixNodeIndex(4));

  (**mIntfVoltage)(2, 0) =
      Math::realFromVectorElement(leftVector, matrixNodeIndex(5));
  */ 

 
}



void EMT::Ph1::ParkTransformer::mnaCompUpdateCurrent(const Matrix &leftVector) {
  /*
  (**mIntfCurrent)(0, 0) = Math::realFromVectorElement(
      **mRightVector, matrixNodeIndex(3));
  (**mIntfCurrent)(1, 0) = Math::realFromVectorElement(
      **mRightVector, matrixNodeIndex(4));
  (**mIntfCurrent)(2, 0) = Math::realFromVectorElement(
      **mRightVector, matrixNodeIndex(5));
  */
}



void EMT::Ph1::ParkTransformer::stampBranchNodeIncidenceMatrix(UInt branchIdx,
    Matrix &branchNodeIncidenceMatrix) {

      

UInt virtual_node0__branchIdx_1 = branchIdx - 10;
UInt virtual_node0__branchIdx_2 = branchIdx - 9;
UInt virtual_node0__branchIdx_3 = branchIdx - 8;
UInt virtual_node1__branchIdx_1 = branchIdx - 7;
UInt virtual_node1__branchIdx_2 = branchIdx - 6;
UInt virtual_node1__branchIdx_3 = branchIdx - 5;
UInt virtual_node1__branchIdx_4 = branchIdx - 4;
UInt virtual_node2__branchIdx_1 = branchIdx - 3;
UInt virtual_node2__branchIdx_2 = branchIdx - 2;
UInt virtual_node2__branchIdx_3 = branchIdx - 1;
UInt virtual_node2__branchIdx_4 = branchIdx;


UInt nodeIdx_phaseA = matrixNodeIndex(0);
UInt nodeIdx_phaseB = matrixNodeIndex(1);
UInt nodeIdx_phaseC = matrixNodeIndex(2);
UInt nodeIdx_alpha = matrixNodeIndex(3);
UInt nodeIdx_beta = matrixNodeIndex(4);
UInt nodeIdx_zero = matrixNodeIndex(5);

UInt virtual_node0_Idx = mVirtualNodes[0]->matrixNodeIndex();
UInt virtual_node1_Idx = mVirtualNodes[1]->matrixNodeIndex();
UInt virtual_node2_Idx = mVirtualNodes[2]->matrixNodeIndex();

//// Virtual node 0

branchNodeIncidenceMatrix(virtual_node0__branchIdx_1, nodeIdx_phaseA) = 1.0;
branchNodeIncidenceMatrix(virtual_node0__branchIdx_1, virtual_node0_Idx) = -1.0; ////////////// duhet -

branchNodeIncidenceMatrix(virtual_node0__branchIdx_2, nodeIdx_alpha) = -1.0; ////////////// duhet -
branchNodeIncidenceMatrix(virtual_node0__branchIdx_2, virtual_node0_Idx) = 1.0;

branchNodeIncidenceMatrix(virtual_node0__branchIdx_3, nodeIdx_zero) = -1.0; ////////////// duhet -
branchNodeIncidenceMatrix(virtual_node0__branchIdx_3, virtual_node0_Idx) = 1.0;

//// Virtual node 1

branchNodeIncidenceMatrix(virtual_node1__branchIdx_1, nodeIdx_phaseB) = 1.0;
branchNodeIncidenceMatrix(virtual_node1__branchIdx_1, virtual_node1_Idx) = -1.0; ////////////// duhet -

branchNodeIncidenceMatrix(virtual_node1__branchIdx_2, nodeIdx_alpha) = 1.0;
branchNodeIncidenceMatrix(virtual_node1__branchIdx_2, virtual_node1_Idx) = -1.0; ////////////// duhet -

branchNodeIncidenceMatrix(virtual_node1__branchIdx_3, nodeIdx_beta) = -1.0; ////////////// duhet -
branchNodeIncidenceMatrix(virtual_node1__branchIdx_3, virtual_node1_Idx) = 1.0;

branchNodeIncidenceMatrix(virtual_node1__branchIdx_4, nodeIdx_zero) = -1.0; ////////////// duhet -
branchNodeIncidenceMatrix(virtual_node1__branchIdx_4, virtual_node1_Idx) = 1.0;

//// Virtual node 2

branchNodeIncidenceMatrix(virtual_node2__branchIdx_1, nodeIdx_phaseC) = 1.0;
branchNodeIncidenceMatrix(virtual_node2__branchIdx_1, virtual_node2_Idx) = -1.0;  ////////////// duhet -

branchNodeIncidenceMatrix(virtual_node2__branchIdx_2, nodeIdx_alpha) = 1.0;
branchNodeIncidenceMatrix(virtual_node2__branchIdx_2, virtual_node2_Idx) = -1.0;  ////////////// duhet -

branchNodeIncidenceMatrix(virtual_node2__branchIdx_3, nodeIdx_beta) = 1.0;
branchNodeIncidenceMatrix(virtual_node2__branchIdx_3, virtual_node2_Idx) = -1.0;  ////////////// duhet -

branchNodeIncidenceMatrix(virtual_node2__branchIdx_4, nodeIdx_zero) = -1.0; ////////////// duhet -
branchNodeIncidenceMatrix(virtual_node2__branchIdx_4, virtual_node2_Idx) = 1.0;


}   




double EMT::Ph1::ParkTransformer::getNumberOfBranches() {
    return 11;
  }




void EMT::Ph1::ParkTransformer::mnaCompAddPreStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes) {
  // add pre-step dependencies of component itself
  prevStepDependencies.push_back(mIntfCurrent);
  prevStepDependencies.push_back(mIntfVoltage);
  modifiedAttributes.push_back(mRightVector);
}



void EMT::Ph1::ParkTransformer::isOmegaConstant(bool isOmegaConstant) {
   mIsOmegaConstant = isOmegaConstant;
}


void EMT::Ph1::ParkTransformer::setTimeStep(Real timeStep) {
  mTimeStep = timeStep;
}


void EMT::Ph1::ParkTransformer::updateOmega() {

  mOmega_prev = *mOmega;
  // Get the new omega from the inertia moment
  Real newOmega = (**(mInertiaMoment->mIntfVoltage))(0, 0);
  *mTheta_init = 0.0;
  setParameters(newOmega, *mTheta_init);

}