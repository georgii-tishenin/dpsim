
#pragma once

#include <iostream>

#include <dpsim-models/Attribute.h>
#include <dpsim-models/AttributeList.h>
#include <dpsim-models/Definitions.h>
#include <dpsim-models/EMT/EMT_Ph1_Inductor.h>
#include <dpsim-models/MNASimPowerComp.h>
#include <dpsim-models/Solver/EigenvalueCompInterface.h>
#include <dpsim-models/Solver/MNAInterface.h>
#include <dpsim-models/Solver/MNAVariableCompInterface.h>

namespace CPS {
namespace EMT {
namespace Ph1 {

class CurrentControlledTorqueSource
    : public MNASimPowerComp<Real>,
      public MNAVariableCompInterface,
      public SharedFactory<CurrentControlledTorqueSource>,
      public EigenvalueCompInterface {
public:
  const typename Attribute<Real>::Ptr mCoefficientCurrent;

  std::shared_ptr<CPS::EMT::Ph1::Inductor> mInductor;

  Real mCoefficient;

  bool mIsNegative = false;

  /// Defines UID, name and logging level
  CurrentControlledTorqueSource(String uid, String name,
                                Logger::Level logLevel = Logger::Level::off);

  /// Defines name and logging level
  CurrentControlledTorqueSource(String name,
                                Logger::Level logLevel = Logger::Level::off)
      : CurrentControlledTorqueSource(name, name, logLevel) {}

  // #### General ####
  /// Sets coefficient
  void setCoefficient(Real coefficient) { mCoefficient = coefficient; };

  void setIsNegative(bool isNegative) { mIsNegative = isNegative; }

  void setInductorForCurrentRerefence(
      const std::shared_ptr<CPS::EMT::Ph1::Inductor> &pt) {
    mInductor = pt;
  }

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

private:
  double calculateAlpha();

  void updateCoefficientCurrent();
};
} // namespace Ph1
} // namespace EMT
} // namespace CPS