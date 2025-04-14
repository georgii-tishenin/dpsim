
#pragma once

#include <iostream>

#include <dpsim-models/AttributeList.h>
#include <dpsim-models/Definitions.h>
#include <dpsim-models/MNASimPowerComp.h>
#include <dpsim-models/Solver/EigenvalueCompInterface.h>
#include <dpsim-models/Solver/MNAInterface.h>


namespace CPS {
namespace EMT {
namespace Ph1 {
class IdealTransformer: public MNASimPowerComp<Real>,
                        public SharedFactory<IdealTransformer>,
                        public EigenvalueCompInterface {
public:
  /// Turns ratio (v1 = v2 * ratio)
  const std::shared_ptr<Real> mN;

  /// Defines UID, name and logging level
  IdealTransformer(String uid, String name,
              Logger::Level logLevel = Logger::Level::off);

/// Defines name and logging level
  IdealTransformer(String name, Logger::Level logLevel = Logger::Level::off)
      : IdealTransformer(name, name, logLevel) {}



// #### General ####
  /// Defines component parameters
  void setParameters(Real N);

  /// Initializes component from power flow data
  void initializeFromNodesAndTerminals(Real frequency) override;

  // #### MNA section ####
  /// Initializes internal variables of the component
  void mnaCompInitialize(Real omega, Real timeStep,
                         Attribute<Matrix>::Ptr leftSideVector) override;
  /// Stamps system matrix
  void mnaCompApplySystemMatrixStamp(SparseMatrixRow &systemMatrix) override;
  /// Stamps right side (source) vector
  void mnaCompApplyRightSideVectorStamp(Matrix &rightVector) override {}
  /// Update interface voltage from MNA system result
  void mnaCompUpdateVoltage(const Matrix &leftVector) override;
  /// Update interface current from MNA system result
  void mnaCompUpdateCurrent(const Matrix &leftVector) override;
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

  double getNumberOfBranches() final;

};
} // namespace Ph1
} // namespace Base
} // namespace CPS
