/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#pragma once

#include <dpsim-models/Base/Base_Ph1_CurrentSource.h>
#include <dpsim-models/MNASimPowerComp.h>
#include <dpsim-models/Solver/MNAInterface.h>
#include <dpsim-models/EMT/EMT_Ph1_PieceWiseNonLinearCharacteristic.h>
#include <dpsim-models/Attribute.h>
#include <dpsim-models/AttributeList.h>

namespace CPS {
namespace EMT {
namespace Ph1 {
/// \brief Ideal current source model
///
/// A positive current is flowing out of
/// node1 and into node2.
class NonLinearInductor : public MNASimPowerComp<Real>,
                      public SharedFactory<NonLinearInductor> {
public:
  const Attribute<Complex>::Ptr mCurrentRef;
  const Attribute<Real>::Ptr mSrcFreq;

  const typename Attribute<Real>::Ptr mFluxLinkage;
  const typename Attribute<Real>::Ptr mInductanceValue;


  Real mFlux = 0.0; // Magnetic flux
  Real mOldVoltage = 0.0;
  Real mNewVoltage = 0.0;
  Real mTimeStep = 0.0; // Needed for the trapezoidal integration
  Real mInductance = 0.0; // Inductance value at the current time step



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

  SimPowerComp<Real>::Ptr clone(String name) override;

  void setParameters(Complex currentRef, Real srcFreq = -1);
  // #### General ####
  /// Initializes component from power flow data
  void initializeFromNodesAndTerminals(Real frequency) override {}

  // #### MNA section ####
  /// Initializes internal variables of the component
  void mnaCompInitialize(Real omega, Real timeStep,
                         Attribute<Matrix>::Ptr leftVector) override;
  /// Stamps system matrix
  void mnaCompApplySystemMatrixStamp(SparseMatrixRow &systemMatrix) override {}
  /// Stamps right side (source) vector
  void mnaCompApplyRightSideVectorStamp(Matrix &rightVector) override;
  ///
  void mnaCompUpdateVoltage(const Matrix &leftVector) override;

  void updateState(Real time);

  void mnaCompPreStep(Real time, Int timeStepCount) override;
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
} // namespace Ph1
} // namespace EMT
} // namespace CPS
