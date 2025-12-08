
#pragma once

#include <iostream>

#include <dpsim-models/Attribute.h>
#include <dpsim-models/AttributeList.h>
#include <dpsim-models/Definitions.h>
#include <dpsim-models/MNASimPowerComp.h>
#include <dpsim-models/Solver/EigenvalueCompInterface.h>
#include <dpsim-models/Solver/MNAInterface.h>
#include <dpsim-models/Solver/MNAVariableCompInterface.h>
#include <dpsim-models/EMT/EMT_Ph1_Inductor.h>

namespace CPS {
namespace EMT {
namespace Ph1 {

class CurrentControlledTorqueSource
    : public MNASimPowerComp<Real>,
      public MNAVariableCompInterface,
      public SharedFactory<CurrentControlledTorqueSource>,
      public EigenvalueCompInterface {
public:
  std::shared_ptr<CPS::EMT::Ph1::Inductor> mInductor;
  const typename Attribute<Real>::Ptr mFlux;
  Real mInductance = 0.0;

  /// Defines UID, name and logging level
  CurrentControlledTorqueSource(String uid, String name,
                                Logger::Level logLevel = Logger::Level::off);

  /// Defines name and logging level
  CurrentControlledTorqueSource(String name,
                                Logger::Level logLevel = Logger::Level::off)
      : CurrentControlledTorqueSource(name, name, logLevel) {}

  // #### General ####
  void setInductor(const std::shared_ptr<CPS::EMT::Ph1::Inductor> &pt) {
    mInductor = pt;
  }

  void setInductance(Real inductance) { mInductance = inductance; }

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