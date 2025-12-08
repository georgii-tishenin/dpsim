#pragma once

#include <dpsim-models/MNASimPowerComp.h>
#include <dpsim-models/Solver/EigenvalueCompInterface.h>
#include <dpsim-models/Solver/MNAInterface.h>

namespace CPS {
template <typename VarType> class SimNode;
}

namespace CPS {
namespace EMT {
namespace Ph1 {
class TimeLaggingVoltageSource : public MNASimPowerComp<Real>,
                                 public SharedFactory<TimeLaggingVoltageSource>,
                                 public EigenvalueCompInterface{
public:
  const Attribute<Real>::Ptr mVoltage;
  int mP;
  int mQ;
  int mIpq;
  std::shared_ptr<CPS::SimNode<Real>> mNodeM;
  std::shared_ptr<CPS::SimNode<Real>> mNodeN;

  TimeLaggingVoltageSource(String uid, String name,
                           Logger::Level logLevel = Logger::Level::off);
  ///
  TimeLaggingVoltageSource(String name,
                           Logger::Level logLevel = Logger::Level::off)
      : TimeLaggingVoltageSource(name, name, logLevel) {}

  void
  setVoltageReferenceNodes(const std::shared_ptr<CPS::SimNode<Real>> &ptM,
                           const std::shared_ptr<CPS::SimNode<Real>> &ptN) {
    mNodeM = ptM;
    mNodeN = ptN;
  }

  // #### MNA section ####
  /// Initializes internal variables of the component
  void mnaCompInitialize(Real omega, Real timeStep,
                         Attribute<Matrix>::Ptr leftVector) override;
  /// Stamps system matrix
  void mnaCompApplySystemMatrixStamp(SparseMatrixRow &systemMatrix) override;
  void mnaCompPreStep(Real time, Int timeStepCount) override;

  /// Add MNA pre step dependencies
  void mnaCompAddPreStepDependencies(
      AttributeBase::List &prevStepDependencies,
      AttributeBase::List &attributeDependencies,
      AttributeBase::List &modifiedAttributes) override;

  // #### Implementation of eigenvalue component interface ####
  void stampBranchNodeIncidenceMatrix(UInt branchIdx,
                                      Matrix &branchNodeIncidenceMatrix) final;
};
} // namespace Ph1
} // namespace EMT
} // namespace CPS
