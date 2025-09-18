#pragma once

#include <dpsim-models/MNASimPowerComp.h>
#include <dpsim-models/Solver/EigenvalueCompInterface.h>
#include <dpsim-models/Solver/MNAInterface.h>

namespace CPS {
namespace EMT {
namespace Ph1 {
class TimeLaggingTorqueSource : public MNASimPowerComp<Real>,
                                 public SharedFactory<TimeLaggingTorqueSource>,
                                 public EigenvalueCompInterface{
public:
  const Attribute<Real>::Ptr mCurrent;
  int mM;
  int mN;
  int mP;
  int mQ;
  int mImn;
  
  TimeLaggingTorqueSource(String uid, String name,
                           Logger::Level logLevel = Logger::Level::off);
  ///
  TimeLaggingTorqueSource(String name,
                           Logger::Level logLevel = Logger::Level::off)
      : TimeLaggingTorqueSource(name, name, logLevel) {}

  // #### MNA section ####
  /// Initializes internal variables of the component
  void mnaCompInitialize(Real omega, Real timeStep,
                         Attribute<Matrix>::Ptr leftVector) override;
  /// Stamps system matrix
  void mnaCompApplySystemMatrixStamp(SparseMatrixRow &systemMatrix) override;
  void mnaCompPreStep(Real time, Int timeStepCount) override;
  void mnaCompPostStep(Real time, Int timeStepCount,
                       Attribute<Matrix>::Ptr &leftVector) override;

  /// Add MNA pre step dependencies
  void mnaCompAddPreStepDependencies(
      AttributeBase::List &prevStepDependencies,
      AttributeBase::List &attributeDependencies,
      AttributeBase::List &modifiedAttributes) override;

  /// Add MNA post step dependencies
  void
  mnaCompAddPostStepDependencies(AttributeBase::List &prevStepDependencies,
                                 AttributeBase::List &attributeDependencies,
                                 AttributeBase::List &modifiedAttributes,
                                 Attribute<Matrix>::Ptr &leftVector) override;

  // #### Implementation of eigenvalue component interface ####
  void stampBranchNodeIncidenceMatrix(UInt branchIdx,
                                      Matrix &branchNodeIncidenceMatrix) final;
};
} // namespace Ph1
} // namespace EMT
} // namespace CPS
