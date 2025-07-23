
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

class Inductor;

class ElectroMechanicalConverter
    : public MNASimPowerComp<Real>,
      public MNAVariableCompInterface,
      public SharedFactory<ElectroMechanicalConverter>,
      public EigenvalueCompInterface {
public:
  /// Turns ratio (v1 = v2 * ratio)
  const std::shared_ptr<Real> mN;

  const typename Attribute<Real>::Ptr mRatio;

  // In order to get the discrete integration of a voltage (flux) as the turns ratio of the ideal transformer

  Real mTimeStep = 0.0;

  Real mOldVoltage = 0.0;

  Real mVoltage = 0.0;

  std::shared_ptr<CPS::EMT::Ph1::Inductor> mInductor;

  void setStatorInductor(const std::shared_ptr<CPS::EMT::Ph1::Inductor> &pt) {
    mInductor = pt;
  }

  /// Defines UID, name and logging level
  ElectroMechanicalConverter(String uid, String name,
                                Logger::Level logLevel = Logger::Level::off);

  /// Defines name and logging level
  ElectroMechanicalConverter(String name,
                                Logger::Level logLevel = Logger::Level::off)
      : ElectroMechanicalConverter(name, name, logLevel) {}

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

  /// Add MNA pre step dependencies
  void mnaCompAddPreStepDependencies(
      AttributeBase::List &prevStepDependencies,
      AttributeBase::List &attributeDependencies,
      AttributeBase::List &modifiedAttributes) override;

  void mnaCompPreStep(Real time, Int timeStepCount) override;

  // #### Implementation of eigenvalue component interface ####
  void stampBranchNodeIncidenceMatrix(UInt branchIdx,
                                      Matrix &branchNodeIncidenceMatrix) final;

  /// Mark that parameter changes so that system matrix is updated
  Bool hasParameterChanged() override { return true; }
};
} // namespace Ph1
} // namespace EMT
} // namespace CPS