// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include <dpsim-models/DP/DP_Ph1_MixedVTypeVariableSSNComp.h>
#include <dpsim-models/Solver/MNAIterativeCompInterface.h>

namespace CPS {
namespace DP {
namespace Ph1 {

/// Averaged grid-following VSI SSN port of EMT::Ph3::AvVoltSourceInverterStateSpace: 8 real control states plus 2 complex envelope states (Vc, If).
class AvVoltSourceInverterStateSpace final
    : public MixedVTypeVariableSSNComp,
      public MNAIterativeCompInterface,
      public SharedFactory<AvVoltSourceInverterStateSpace> {
private:
  static constexpr Real mIterationRelativeTolerance = 1e-9;
  static constexpr Real mIterationAbsoluteTolerance = 1e-9;

  enum StateIndex : Int {
    Psi = 0,
    PhiPLL = 1,
    PFiltered = 2,
    QFiltered = 3,
    PhiD = 4,
    PhiQ = 5,
    GammaD = 6,
    GammaQ = 7,
    VcRe = 8,
    VcIm = 9,
    IfRe = 10,
    IfIm = 11
  };

  Real mLf;
  Real mCf;
  Real mRf;
  Real mRc;

  Real mOmegaN;
  Real mKpPLL;
  Real mKiPLL;

  Real mOmegaCutoff;
  Real mPRef;
  Real mQRef;
  Real mKpPowerCtrl;
  Real mKiPowerCtrl;
  Real mKpCurrCtrl;
  Real mKiCurrCtrl;

  const Attribute<Real>::Ptr mVcD;
  const Attribute<Real>::Ptr mVcQ;
  const Attribute<Real>::Ptr mIrcD;
  const Attribute<Real>::Ptr mIrcQ;
  const Attribute<Real>::Ptr mPInst;
  const Attribute<Real>::Ptr mQInst;
  const Attribute<Real>::Ptr mOmegaPLL;

  Bool mIterativeSolutionEnabled = false;
  Matrix mIterationState;
  Matrix mIterationInput;
  const Attribute<Int>::Ptr mIterationCount;
  const Attribute<Real>::Ptr mIterationStateResidual;
  const Attribute<Real>::Ptr mIterationInputResidual;

  /// Builds the affine real model (A,B,C,D,E,F) around (x,u): RHS + analytic Jacobian, E = f(x,u) - A*x - B*u.
  void buildStateSpaceModel(const Matrix &x, const Matrix &u, Matrix &A,
                            Matrix &B, Matrix &C, Matrix &D, Matrix &E,
                            Matrix &F) const;

  Bool updateComponentParameters(const Matrix &state, const Matrix &input);
  Real iterationResidual(const Matrix &value,
                         const Matrix &previousValue) const;

protected:
  Bool updateComponentParameters() override final;
  void updateLogAttributes(const Matrix &u) const override final;

public:
  using SharedFactory<AvVoltSourceInverterStateSpace>::make;

  AvVoltSourceInverterStateSpace(String uid, String name,
                                 Logger::Level logLevel = Logger::Level::off);
  AvVoltSourceInverterStateSpace(String name,
                                 Logger::Level logLevel = Logger::Level::off)
      : AvVoltSourceInverterStateSpace(name, name, logLevel) {}

  void setParameters(Real lf, Real cf, Real rf, Real rc, Real omegaN,
                     Real kpPLL, Real kiPLL, Real omegaCutoff, Real pRef,
                     Real qRef, Real kpPowerCtrl, Real kiPowerCtrl,
                     Real kpCurrCtrl, Real kiCurrCtrl);

  /// Enable repeated MNA solutions within each time step. Disabled by
  /// default to preserve the existing execution behavior.
  void setIterativeSolution(Bool enabled) {
    mIterativeSolutionEnabled = enabled;
  }

  Bool iterativeSolutionEnabled() const { return mIterativeSolutionEnabled; }

  void initializeFromNodesAndTerminals(Real frequency) override;

  void mnaInitializeIteration(Real time, Int timeStepCount) override final;
  MNAIterationUpdate
  mnaUpdateIteration(const Matrix &leftVector) override final;
  void mnaFinalizeIteration() override final;
};

} // namespace Ph1
} // namespace DP
} // namespace CPS
