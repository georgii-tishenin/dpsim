/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#pragma once

#include <stdexcept>

#include <dpsim-models/EMT/EMT_PassiveOperationProfile.h>

#include <dpsim-models/CompositePowerComp.h>
#include <dpsim-models/EMT/EMT_Ph3_Capacitor.h>
#include <dpsim-models/EMT/EMT_Ph3_Inductor.h>
#include <dpsim-models/EMT/EMT_Ph3_Resistor.h>
#include <dpsim-models/Solver/MNAInterface.h>

namespace CPS {
namespace EMT {
namespace Ph3 {

class Shunt : public CompositePowerComp<Real>, public SharedFactory<Shunt> {
private:
  /// Three-phase current from the terminal node into the shunt.
  Attribute<UInt>::Ptr mOutputUpdates;
  Attribute<Matrix>::Ptr mTerminalCurrent0;
  /// Conductance [S]
  Matrix mConductance;
  /// Susceptance [S]
  Matrix mSusceptance;

  /// Capacitor between terminal and ground
  std::shared_ptr<Capacitor> mSubCapacitor;
  /// Inductor between terminal and ground
  std::shared_ptr<Inductor> mSubInductor;
  /// Resistor between terminal and ground
  std::shared_ptr<Resistor> mSubResistor;
  Bool mUseSubcomponents = true;
  Bool mUseElementarySSN = false;
  Matrix mCapEquivalentConductance = Matrix::Zero(3, 3);
  Matrix mIndEquivalentConductance = Matrix::Zero(3, 3);
  Matrix mCapEquivalentCurrent = Matrix::Zero(3, 1);
  Matrix mIndEquivalentCurrent = Matrix::Zero(3, 1);
  Matrix mCapBranchCurrent = Matrix::Zero(3, 1);
  Matrix mIndBranchCurrent = Matrix::Zero(3, 1);

  PassiveOperationProfile mOperationProfile;

public:
  /// Select discrete Norton SSNs for the existing R/L/C circuit children.
  /// Must be called before initialization; physical nodes remain unchanged.
  void setUseElementarySSN(Bool enabled) {
    if (mSubCompCreated || !mnaSubComponents().empty())
      throw std::logic_error("Select elementary SSNs before initialization.");
    if (enabled && !mUseSubcomponents)
      throw std::invalid_argument(
          "Elementary SSNs require the subcomponent circuit.");
    mUseElementarySSN = enabled;
  }
  /// Number of actual R/L/C children, excluding ideal transformer constraints.
  UInt primitiveCount() const {
    return static_cast<UInt>(mnaSubComponents().size());
  }

  void mnaSetDetailedProfiling(Bool enabled) override {
    mOperationProfile.enabled = enabled;
    for (const auto &component : mnaSubComponents())
      component->mnaSetDetailedProfiling(enabled);
  }
  void mnaWriteProfile(std::ostream &profile) const override {
    mOperationProfile.write(profile, **mName);
    for (const auto &component : mnaSubComponents())
      component->mnaWriteProfile(profile);
  }

  /// Defines UID, name, component parameters and logging level
  Shunt(String uid, String name, Logger::Level logLevel = Logger::Level::off);

  /// Defines name and logging level
  Shunt(String name, Logger::Level logLevel = Logger::Level::off)
      : Shunt(name, name, logLevel) {}
  /// Direct MNA implementation that preserves the legacy constructor by default.
  Shunt(String uid, String name, Bool useSubcomponents,
        Logger::Level logLevel = Logger::Level::off);
  Shunt(String name, Bool useSubcomponents,
        Logger::Level logLevel = Logger::Level::off)
      : Shunt(name, name, useSubcomponents, logLevel) {}

  // #### General ####
  /// Set shunt specific parameters
  void setParameters(Real conductance, Real susceptance);
  void setParameters(Matrix conductance, Matrix susceptance);

  // #### MNA section ####
  /// Initializes component from power flow data
  void initializeParentFromNodesAndTerminals(Real frequency) final;
  void mnaParentInitialize(Real omega, Real timeStep,
                           Attribute<Matrix>::Ptr leftVector) override;
  void mnaParentApplySystemMatrixStamp(SparseMatrixRow &systemMatrix) override;
  void mnaParentApplyRightSideVectorStamp(Matrix &rightVector) override;
  /// Add MNA pre step dependencies
  void mnaParentAddPreStepDependencies(
      AttributeBase::List &prevStepDependencies,
      AttributeBase::List &attributeDependencies,
      AttributeBase::List &modifiedAttributes) final;
  /// MNA pre step operations
  void mnaParentPreStep(Real time, Int timeStepCount) final;
  /// Updates internal current variable of the component
  void mnaCompUpdateCurrent(const Matrix &leftVector) final;
  /// Updates internal voltage variable of the component
  void mnaCompUpdateVoltage(const Matrix &leftVector) final;
  /// MNA post-step operations
  void mnaParentPostStep(Real time, Int timeStepCount,
                         Attribute<Matrix>::Ptr &leftVector) final;
  /// add MNA post-step dependencies
  void
  mnaParentAddPostStepDependencies(AttributeBase::List &prevStepDependencies,
                                   AttributeBase::List &attributeDependencies,
                                   AttributeBase::List &modifiedAttributes,
                                   Attribute<Matrix>::Ptr &leftVector) final;
};
} // namespace Ph3
} // namespace EMT
} // namespace CPS
