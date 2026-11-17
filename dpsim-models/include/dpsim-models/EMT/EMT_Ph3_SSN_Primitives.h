// SPDX-License-Identifier: MPL-2.0
#pragma once

#include <dpsim-models/EMT/EMT_Ph3_Capacitor.h>
#include <dpsim-models/EMT/EMT_Ph3_Inductor.h>
#include <dpsim-models/EMT/EMT_Ph3_Resistor.h>

#include <array>
#include <stdexcept>
#include <type_traits>

namespace CPS::EMT::Ph3::SSN {

/// Inspection interface for the elementary discrete-time Norton SSNs.
class ElementarySSNInterface {
public:
  virtual ~ElementarySSNInterface() = default;
  /// Stored companion memories, not a minimal continuous-time state order.
  virtual UInt getStateCount() const = 0;
};

/// An elementary R, L or C represented by precomputed discrete Norton maps.
///
/// For L and C, x=[v_previous; i_previous], h=H*x, i=D*v+h, followed by
/// x=[v; i]. Trapezoidal integration gives
///   L: D=dt/2*L^-1, H=[ D,  I];
///   C: D=2/dt*C,    H=[-D, -I].
/// A resistor has no state/history and D=R^-1. Each primitive retains the
/// physical terminals of the native companion circuit. In particular, C
/// needs no extra current unknown.
///
/// The native base supplies parameter attributes, phasor initialization and
/// the established public component interface. Its MNA companion algorithms
/// are all replaced below; existing equipment wiring can therefore be reused.
template <class Native>
class NortonPrimitive final : public Native, public ElementarySSNInterface {
  static constexpr Bool IsResistor = std::is_same_v<Native, Ph3::Resistor>;
  static constexpr Bool IsInductor = std::is_same_v<Native, Ph3::Inductor>;
  static constexpr Bool IsCapacitor = std::is_same_v<Native, Ph3::Capacitor>;
  static_assert(IsResistor || IsInductor || IsCapacitor);

  Attribute<Matrix>::Ptr mState;
  Matrix *mStateValue = nullptr, *mVoltageValue = nullptr,
         *mCurrentValue = nullptr;
  Matrix mD = Matrix::Zero(3, 3), mH;
  Matrix mHistory = Matrix::Zero(3, 1);
  std::array<UInt, 6> mRows{};
  std::array<Real, 6> mSigns{};
  UInt mRowCount = 0;

  const Matrix &parameter() const {
    if constexpr (IsResistor)
      return **this->mResistance;
    else if constexpr (IsInductor)
      return **this->mInductance;
    else
      return **this->mCapacitance;
  }

public:
  NortonPrimitive(String uid, String name,
                  Logger::Level logLevel = Logger::Level::off)
      : Native(uid, name, logLevel) {
    mState = this->mAttributes->template create<Matrix>("ssn_x");
    **mState = Matrix::Zero(getStateCount(), 1);
    mH = Matrix::Zero(3, getStateCount());
  }
  NortonPrimitive(String name, Logger::Level logLevel = Logger::Level::off)
      : NortonPrimitive(name, name, logLevel) {}

  UInt getStateCount() const override { return IsResistor ? 0 : 6; }
  const Matrix &getNortonD() const { return mD; }
  const Matrix &getHistoryMap() const { return mH; }

  SimPowerComp<Real>::Ptr clone(String name) override {
    auto copy = std::make_shared<NortonPrimitive>(name, this->mLogLevel);
    copy->setParameters(parameter());
    return copy;
  }

  void mnaCompInitialize(Real, Real timeStep, Attribute<Matrix>::Ptr) override {
    if (timeStep <= 0. || parameter().rows() != 3 || parameter().cols() != 3)
      throw std::invalid_argument(
          "Elementary SSN requires dt>0 and a 3x3 parameter.");
    this->updateMatrixNodeIndices();
    mStateValue = &mState->get();
    mVoltageValue = &this->mIntfVoltage->get();
    mCurrentValue = &this->mIntfCurrent->get();
    if constexpr (IsResistor) {
      mD = parameter().inverse();
      **this->mRightVector = Matrix::Zero(0, 0);
    } else {
      if constexpr (IsInductor)
        mD = (timeStep / 2.) * parameter().inverse();
      else
        mD = (2. * parameter()) / timeStep;
      const Real sign = IsInductor ? 1. : -1.;
      mH.leftCols(3) = sign * mD;
      mH.rightCols(3) = sign * Matrix::Identity(3, 3);
      mStateValue->topRows(3) = *mVoltageValue;
      mStateValue->bottomRows(3) = *mCurrentValue;
      mHistory.noalias() = mH * *mStateValue;
      // Keep the native public conductance getter meaningful for inspection.
      this->mEquivCond = mD;
    }
    if (!mD.allFinite())
      throw std::invalid_argument(
          "Elementary SSN parameter produces a non-finite Norton map.");
    mRowCount = 0;
    for (UInt terminal = 0; terminal < 2; ++terminal)
      if (this->terminalNotGrounded(terminal))
        for (UInt phase = 0; phase < 3; ++phase) {
          mRows[mRowCount] = this->matrixNodeIndex(terminal, phase);
          mSigns[mRowCount++] = terminal == 0 ? 1. : -1.;
        }
  }

  void mnaCompApplySystemMatrixStamp(SparseMatrixRow &matrix) override {
    // Use the native stamp; the MNA solver prunes exact static zeros for all forms.
    MNAStampUtils::stampConductanceMatrix(
        mD, matrix, this->matrixNodeIndex(0), this->matrixNodeIndex(1),
        this->terminalNotGrounded(0), this->terminalNotGrounded(1),
        this->mSLog);
  }

  void mnaCompApplyRightSideVectorStamp(Matrix &rightVector) override {
    if constexpr (!IsResistor) {
      PassiveOperationProfile::Scope operation(
          this->mOperationProfile, PassiveOperationProfile::History);
      mHistory.noalias() = mH * *mStateValue;
      operation.next(PassiveOperationProfile::RhsLocal);
      for (UInt row = 0; row < mRowCount; ++row)
        rightVector(mRows[row], 0) = mSigns[row] * mHistory(row % 3, 0);
    }
  }
  Bool mnaHasCompactRightSideVectorStamp() const override { return true; }
  void mnaAddCompactRightSideVectorStamp(Matrix &rightVector) override {
    if constexpr (!IsResistor)
      for (UInt row = 0; row < mRowCount; ++row)
        rightVector(mRows[row], 0) += (*this->mRightVectorValue)(mRows[row], 0);
  }
  void mnaCompPreStep(Real, Int) override {
    mnaCompApplyRightSideVectorStamp(*this->mRightVectorValue);
  }

  void mnaCompUpdateVoltage(const Matrix &leftVector) override {
    for (UInt phase = 0; phase < 3; ++phase) {
      Real voltage = 0.;
      if (this->terminalNotGrounded(1))
        voltage += leftVector(this->matrixNodeIndex(1, phase), 0);
      if (this->terminalNotGrounded(0))
        voltage -= leftVector(this->matrixNodeIndex(0, phase), 0);
      (*mVoltageValue)(phase, 0) = voltage;
    }
  }
  void mnaCompUpdateCurrent(const Matrix &) override {
    mCurrentValue->noalias() = mD * *mVoltageValue;
    if constexpr (!IsResistor)
      *mCurrentValue += mHistory;
  }
  void mnaCompPostStep(Real, Int, Attribute<Matrix>::Ptr &leftVector) override {
    PassiveOperationProfile::Scope operation(
        this->mOperationProfile, PassiveOperationProfile::VoltageRecovery);
    mnaCompUpdateVoltage(**leftVector);
    operation.next(PassiveOperationProfile::StateAdvance);
    mnaCompUpdateCurrent(**leftVector);
    if constexpr (!IsResistor) {
      mStateValue->topRows(3) = *mVoltageValue;
      mStateValue->bottomRows(3) = *mCurrentValue;
    }
  }
  void mnaCompAddPreStepDependencies(AttributeBase::List &previous,
                                     AttributeBase::List &,
                                     AttributeBase::List &modified) override {
    if constexpr (!IsResistor) {
      previous.push_back(mState);
      modified.push_back(this->mRightVector);
    }
  }
  void
  mnaCompAddPostStepDependencies(AttributeBase::List &,
                                 AttributeBase::List &dependencies,
                                 AttributeBase::List &modified,
                                 Attribute<Matrix>::Ptr &leftVector) override {
    dependencies.push_back(leftVector);
    modified.push_back(this->mIntfVoltage);
    modified.push_back(this->mIntfCurrent);
    if constexpr (!IsResistor)
      modified.push_back(mState);
  }

  // The new representation currently targets MNA. Never fall through to the
  // inherited tear-specific native history algorithm with a separate state.
  void mnaTearInitialize(Real, Real) {
    throw std::logic_error(
        "Elementary Norton SSNs currently support MNA only.");
  }
  void mnaTearApplyMatrixStamp(SparseMatrixRow &) {
    throw std::logic_error(
        "Elementary Norton SSNs currently support MNA only.");
  }
};

using ElementaryResistor = NortonPrimitive<Ph3::Resistor>;
using ElementaryInductor = NortonPrimitive<Ph3::Inductor>;
using ElementaryCapacitor = NortonPrimitive<Ph3::Capacitor>;

} // namespace CPS::EMT::Ph3::SSN

namespace CPS::EMT::Ph3 {

/// Select only the numerical primitive; the equipment owns the physical wiring.
template <class Native>
std::shared_ptr<Native> makePassivePrimitive(Bool elementarySSN, String name,
                                             Logger::Level logLevel) {
  if (elementarySSN)
    return std::make_shared<SSN::NortonPrimitive<Native>>(name, logLevel);
  return std::make_shared<Native>(name, logLevel);
}

} // namespace CPS::EMT::Ph3
