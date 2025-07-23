
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

class InertiaMoment; // Forward declaration

class ParkTransformer : public MNASimPowerComp<Real>,
                        public MNAVariableCompInterface,
                        public SharedFactory<ParkTransformer>,
                        public EigenvalueCompInterface {

public:
  // The rotating frame frequency
  const std::shared_ptr<Real> mOmega;
  const std::shared_ptr<Real> mTheta_init;
  const std::shared_ptr<Real> mTheta;

  const typename Attribute<Real>::Ptr mTheta_atr;

  Real mOmega_prev = 0.0;

  Real mTimeStep = 0.0;

  bool mIsOmegaConstant = false;

  std::shared_ptr<CPS::EMT::Ph1::InertiaMoment> mInertiaMoment;

  void
  setInertiaMoment(const std::shared_ptr<CPS::EMT::Ph1::InertiaMoment> &pt) {
    mInertiaMoment = pt;
  }

  void updateOmega();

  /// Defines UID, name and logging level
  ParkTransformer(String uid, String name,
                  Logger::Level logLevel = Logger::Level::off);

  /// Defines name and logging level
  ParkTransformer(String name, Logger::Level logLevel = Logger::Level::off)
      : ParkTransformer(name, name, logLevel) {}

  // #### General ####
  /// Defines component parameters
  void setParameters(Real omega, Real theta_initial);

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
  void
  mnaCompUpdateVoltage(const Matrix &leftVector) override; // Idk if I need this
  /// Update interface current from MNA system result
  void mnaCompUpdateCurrent(
      const Matrix &leftVector) override; // Idk if I need this as well
  void mnaCompPostStep(Real time, Int timeStepCount,
                       Attribute<Matrix>::Ptr &leftVector) override;

  /// MNA pre step operations
  void mnaCompPreStep(Real time, Int timeStepCount) override;

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
  /// Mark that parameter changes so that system matrix is updated
  Bool hasParameterChanged() override { return true; }

  // #### Implementation of eigenvalue component interface ####
  void stampBranchNodeIncidenceMatrix(UInt branchIdx,
                                      Matrix &branchNodeIncidenceMatrix) final;

  void isOmegaConstant(bool isOmegaConstant);

  void setTimeStep(Real timeStep);
};
} // namespace Ph1
} // namespace EMT
} // namespace CPS