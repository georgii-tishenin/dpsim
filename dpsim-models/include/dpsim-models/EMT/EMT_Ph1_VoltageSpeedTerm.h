
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

class VoltageSpeedTerm : public MNASimPowerComp<Real>,
                         public MNAVariableCompInterface,
                         public SharedFactory<VoltageSpeedTerm>,
                         public EigenvalueCompInterface {
public:
  const typename Attribute<Real>::Ptr mOmega;

  std::shared_ptr<CPS::SimNode<Real>> mOmegaReferenceNode;

  Real mOmegaOffset = 0.0;

  Real mInductance;

  bool mIsConstantSpeed = false;

  bool mIsNegative = false;

  /// Defines UID, name and logging level
  VoltageSpeedTerm(String uid, String name,
                   Logger::Level logLevel = Logger::Level::off);

  /// Defines name and logging level
  VoltageSpeedTerm(String name, Logger::Level logLevel = Logger::Level::off)
      : VoltageSpeedTerm(name, name, logLevel) {}

  // #### General ####
  /// Sets initial omega
  void setInitialOmega(Real omega);

  void setOmegaOffset(Real omegaOffset) { mOmegaOffset = omegaOffset; }

  void setIsConstantSpeed(bool isConstantSpeed) {
    mIsConstantSpeed = isConstantSpeed;
  }

  void setIsNegative(bool isNegative) { mIsNegative = isNegative; }

  /// Sets inductance
  void setInductance(Real inductance);

  void setOmegaReferenceNode(const std::shared_ptr<CPS::SimNode<Real>> &pt) {
    mOmegaReferenceNode = pt;
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