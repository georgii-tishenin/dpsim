/* Copyright 2017-2020 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/
#pragma once

#include <dpsim-models/MNASimPowerComp.h>
#include <dpsim-models/EMT/EMT_Ph1_PieceWiseNonLinearCharacteristic.h>
#include <dpsim-models/Solver/MNAInterface.h>
#include <dpsim-models/Attribute.h>
#include <dpsim-models/AttributeList.h>

namespace CPS {
namespace EMT {
namespace Ph3 {
/// \brief 3-phase non-linear inductor model
///
/// This model uses modified nodal analysis to represent a 3-phase non-linear inductor model as current sources that inject a current every time step
/// This involves the stamping of the current to the right side vector.
class NonLinearInductor : public MNASimPowerComp<Real>,
                      public SharedFactory<NonLinearInductor> {

protected:
  void updateCurrent(Real time);

public:

  Matrix mFlux = Matrix::Zero(3, 1); ///< Flux linkage vector
  Matrix mOldVoltages = Matrix::Zero(3, 1); ///< Old voltage vector
  Matrix mNewVoltages = Matrix::Zero(3, 1); ///< New voltage vector
  Real mTimeStep = 0.0; ///< Time step for trapezoidal


  std::shared_ptr<CPS::EMT::Ph1::PieceWiseNonLinearCharacteristic> mPieceWiseCharacteristic;

  void setPieceWiseCharacteristic(const std::shared_ptr<CPS::EMT::Ph1::PieceWiseNonLinearCharacteristic>& pt) {
    mPieceWiseCharacteristic = pt;
    }

  /// Defines UID, name and logging level
  NonLinearInductor(String uid, String name,
                Logger::Level logLevel = Logger::Level::off);
  ///
  NonLinearInductor(String name, Logger::Level logLevel = Logger::Level::off)
      : NonLinearInductor(name, name, logLevel) {}

  // #### General ####
  /// Initializes component from power flow data
  void initializeFromNodesAndTerminals(Real frequency) override;

  // #### MNA section ####
  /// Initializes internal variables of the component
  void mnaCompInitialize(Real omega, Real timeStep,
                         Attribute<Matrix>::Ptr leftVector) override;
  /// Stamps right side (source) vector
  void mnaCompApplyRightSideVectorStamp(Matrix &rightVector) override;
  /// Returns voltage through the component
  void mnaCompUpdateVoltage(const Matrix &leftVector) override;
  /// MNA pre step operations
  void mnaCompPreStep(Real time, Int timeStepCount) override;
  /// MNA post step operations
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


  void setTimeStep(Real timeStep);

};



} // namespace Ph3
} // namespace EMT
} // namespace CPS
