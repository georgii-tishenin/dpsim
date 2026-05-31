
#pragma once

#include <iostream>

#include <dpsim-models/Attribute.h>
#include <dpsim-models/AttributeList.h>
#include <dpsim-models/Definitions.h>
#include <dpsim-models/MNASimPowerComp.h>
#include <dpsim-models/Solver/EigenvalueCompInterface.h>
#include <dpsim-models/Solver/MNAInterface.h>
#include <dpsim-models/Solver/MNAVariableCompInterface.h>

namespace CPS {
template <typename VarType> class SimNode;
}

namespace CPS {
namespace EMT {
namespace Ph1 {

class ThetaControlledVoltageSource
    : public MNASimPowerComp<Real>,
      public MNAVariableCompInterface,
      public SharedFactory<ThetaControlledVoltageSource>,
      public EigenvalueCompInterface {
public:

  /// flux is turns ratio (v = flux * omega, torque = flux * i)
  const typename Attribute<Real>::Ptr mFlux;
  
  std::shared_ptr<CPS::SimNode<Real>> mVoltageReferenceNode;

  std::shared_ptr<CPS::SimNode<Real>> mPhaseAReferenceNode;

  std::shared_ptr<CPS::SimNode<Real>> mPhaseBReferenceNode;

  std::shared_ptr<CPS::SimNode<Real>> mPhaseCReferenceNode;

  Real mCoefficient;

  // In order to get the discrete integration of a voltage (flux)
  Real mTimeStep = 0.0;

  Real mV_a = 0.0;

  Real mV_b = 0.0;

  Real mV_c = 0.0;

  /// Defines UID, name and logging level
  ThetaControlledVoltageSource(String uid, String name,
                                Logger::Level logLevel = Logger::Level::off);

  /// Defines name and logging level
  ThetaControlledVoltageSource(String name,
                                Logger::Level logLevel = Logger::Level::off)
      : ThetaControlledVoltageSource(name, name, logLevel) {}

  // #### General ####
  /// Sets coefficient
  void setCoefficient(Real coefficient) { mCoefficient = coefficient; };

  void setInitialFlux(Real flux);

  void setVoltageReferenceNode(const std::shared_ptr<CPS::SimNode<Real>> &pt) {
    mVoltageReferenceNode = pt;
  }

  void setPhaseABCReferenceNode( 
      const std::shared_ptr<CPS::SimNode<Real>> &ptA,
      const std::shared_ptr<CPS::SimNode<Real>> &ptB,
      const std::shared_ptr<CPS::SimNode<Real>> &ptC) {
    mPhaseAReferenceNode = ptA;
    mPhaseBReferenceNode = ptB;
    mPhaseCReferenceNode = ptC;
  }

  // #### MNA section ####
  /// Initializes internal variables of the component
  void mnaCompInitialize(Real omega, Real timeStep,
                         Attribute<Matrix>::Ptr leftSideVector) override;
  /// Stamps system matrix
  void mnaCompApplySystemMatrixStamp(SparseMatrixRow &systemMatrix) override;

  void mnaCompPostStep(Real time, Int timeStepCount,
                       Attribute<Matrix>::Ptr &leftVector) override;
  /// Add MNA post step dependencies
  void
  mnaCompAddPostStepDependencies(AttributeBase::List &prevStepDependencies,
                                 AttributeBase::List &attributeDependencies,
                                 AttributeBase::List &modifiedAttributes,
                                 Attribute<Matrix>::Ptr &leftVector) override;

  // #### Implementation of eigenvalue component interface ####
  void stampBranchNodeIncidenceMatrix(UInt branchIdx,
                                      Matrix &branchNodeIncidenceMatrix) final;

  // Mark that parameter changes so that system matrix is updated
  Bool hasParameterChanged() override { return true; }

};
} // namespace Ph1
} // namespace EMT
} // namespace CPS