/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#pragma once

#include <dpsim-models/Base/Base_Ph1_Inductor.h>
#include <dpsim-models/MNASimPowerComp.h>
#include <dpsim-models/Solver/MNATearInterface.h>

namespace CPS {
namespace SP {
namespace Ph1 {
/// Static phasor inductor model
class Inductor : public MNASimPowerComp<Complex>,
                 public Base::Ph1::Inductor,
                 public MNATearInterface,
                 public SharedFactory<Inductor> {
protected:
  /// susceptance [S]
  Complex mSusceptance;

public:
  /// Defines UID, name and log level
  Inductor(String uid, String name,
           Logger::Level logLevel = Logger::Level::off);
  /// Defines name and log level
  Inductor(String name, Logger::Level logLevel = Logger::Level::off)
      : Inductor(name, name, logLevel) {}

  SimPowerComp<Complex>::Ptr clone(String name) override;

  // #### General ####
  /// Initializes component from power flow data
  void initializeFromNodesAndTerminals(Real frequency) override;
  // #### MNA section ####
  /// Initializes internal variables of the component
  void mnaCompInitialize(Real omega, Real timeStep,
                         Attribute<Matrix>::Ptr leftVector) override;
  /// Stamps system matrix
  void mnaCompApplySystemMatrixStamp(SparseMatrixRow &systemMatrix) override;
  /// Update interface voltage from MNA system results
  void mnaCompUpdateVoltage(const Matrix &leftVector) override;
  /// Update interface current from MNA system results
  void mnaCompUpdateCurrent(const Matrix &leftVector) override;
  /// MNA post step operations
  void mnaCompPostStep(Real time, Int timeStepCount,
                       Attribute<Matrix>::Ptr &leftVector) override;
  /// Add MNA post step dependencies
  void
  mnaCompAddPostStepDependencies(AttributeBase::List &prevStepDependencies,
                                 AttributeBase::List &attributeDependencies,
                                 AttributeBase::List &modifiedAttributes,
                                 Attribute<Matrix>::Ptr &leftVector) override;
  //
  void mnaTearApplyMatrixStamp(SparseMatrixRow &tearMatrix) override;
};
} // namespace Ph1
} // namespace SP
} // namespace CPS
