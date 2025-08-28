
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
namespace EMT {
namespace Ph1 {

class VoltageSpeedTerm : public MNASimPowerComp<Real>,
                         public MNAVariableCompInterface,
                         public SharedFactory<VoltageSpeedTerm>,
                         public EigenvalueCompInterface {
public:
  const typename Attribute<Real>::Ptr mOmega;

  Real mInductance;

  bool mIsNegative = false;

  void setIsNegative(bool isNegative) { mIsNegative = isNegative; }

  /// Defines UID, name and logging level
  VoltageSpeedTerm(String uid, String name,
                   Logger::Level logLevel = Logger::Level::off);

  /// Defines name and logging level
  VoltageSpeedTerm(String name, Logger::Level logLevel = Logger::Level::off)
      : VoltageSpeedTerm(name, name, logLevel) {}

  // #### General ####
  /// Sets initial flux
  void setInitialOmega(Real omega);

  /// Sets inductance
  void setInductance(Real inductance);

  // #### MNA section ####
  /// Initializes internal variables of the component
  void mnaCompInitialize(Real omega, Real timeStep,
                         Attribute<Matrix>::Ptr leftSideVector) override;
  /// Stamps system matrix
  void mnaCompApplySystemMatrixStamp(SparseMatrixRow &systemMatrix) override;

  // #### Implementation of eigenvalue component interface ####
  void stampBranchNodeIncidenceMatrix(UInt branchIdx,
                                      Matrix &branchNodeIncidenceMatrix) final;

  // Mark that parameter changes so that system matrix is updated
  Bool hasParameterChanged() override { return true; }
};
} // namespace Ph1
} // namespace EMT
} // namespace CPS