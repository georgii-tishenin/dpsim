
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

class ElectroMechanicalConverter
    : public MNASimPowerComp<Real>,
      public MNAVariableCompInterface,
      public SharedFactory<ElectroMechanicalConverter>,
      public EigenvalueCompInterface {
public:
  /// flux is turns ratio (v = flux * omega, torque = flux * i)
  const typename Attribute<Real>::Ptr mFlux;

  // In order to get the discrete integration of a voltage (flux) as the turns ratio of the ideal transformer

  Real mTimeStep = 0.0;

  Real mOldVoltage = 0.0;

  Real mVoltage = 0.0;

  std::shared_ptr<CPS::SimNode<Real>> mVoltageReferenceNode;

  void setVoltageReferenceNode(const std::shared_ptr<CPS::SimNode<Real>> &pt) {
    mVoltageReferenceNode = pt;
  }

  /// Defines UID, name and logging level
  ElectroMechanicalConverter(String uid, String name,
                             Logger::Level logLevel = Logger::Level::off);

  /// Defines name and logging level
  ElectroMechanicalConverter(String name,
                             Logger::Level logLevel = Logger::Level::off)
      : ElectroMechanicalConverter(name, name, logLevel) {}

  // #### General ####
  /// Sets initial flux
  void setInitialFlux(Real flux);

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

  /// Mark that parameter changes so that system matrix is updated
  Bool hasParameterChanged() override { return true; }
};
} // namespace Ph1
} // namespace EMT
} // namespace CPS