/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#pragma once

#include <dpsim-models/MNASimPowerComp.h>
#include <dpsim-models/Solver/EigenvalueCompInterface.h>
#include <dpsim-models/Solver/MNAInterface.h>


namespace CPS {
namespace EMT {
namespace Ph1 {


  
    class Inductor;  // Forward declaration so that we can get the current from the resistor
    class InertiaMoment;  // Forward declaration so that we can get the angular speed from the inertia moment



/// \brief Cross-Coupling source model
///
/// This model uses modified nodal analysis to represent a cross-coupling term (speed term) as an ideal voltage source.
/// For a voltage source between nodes j and k, a new variable (current across the voltage source)
/// is added to the left side vector
/// as unkown and it is taken into account for the equation of node j as positve and for the equation
/// of node k as negative. Moreover
/// a new equation ej - ek = V is added to the problem.
class CrossCoupling : public MNASimPowerComp<Real>,
                      public SharedFactory<CrossCoupling>,
                      public EigenvalueCompInterface {
private:
  Real mTimeStep;
  Real mConstantSpeed = 0;  // Default value for constant speed

protected:
  void updateVoltage(Real time);

public:
  const Attribute<Real>::Ptr mAngularSpeedRef;
  /// Defines UID, name and logging level
  CrossCoupling(String uid, String name,
                Logger::Level logLevel = Logger::Level::off);
  ///
  CrossCoupling(String name, Logger::Level logLevel = Logger::Level::off)
      : CrossCoupling(name, name, logLevel) {}


    /// Pointer to the resistor component from which the speed term shall get its current
  std::shared_ptr<CPS::EMT::Ph1::Inductor> mInductor;

  /// Pointer to the inertia moment component from which the angular speed is derived
  std::shared_ptr<CPS::EMT::Ph1::InertiaMoment> mInertiaMoment;

void setInertiaMoment(const std::shared_ptr<CPS::EMT::Ph1::InertiaMoment>& pt) {
    mInertiaMoment = pt;
  }

void setResistor(const std::shared_ptr<CPS::EMT::Ph1::Inductor>& pt) {
    mInductor = pt;
}



  Real mOldVoltage = 0.0;

  Real mVoltage = 0.0;

  Real mFluxLinkage = 0.0;

  bool mNegativeSpeedTermVoltageFlag = false;

  void setVoltageNegative(bool flag) { 
    mNegativeSpeedTermVoltageFlag = flag; 
  }







  void setParameters(Real omega);

  SimPowerComp<Real>::Ptr clone(String name) override;
  // #### General ####
  /// Initializes component from power flow data
  void initializeFromNodesAndTerminals(Real frequency) override {}

  // #### MNA section ####
  /// Initializes internal variables of the component
  void mnaCompInitialize(Real omega, Real timeStep,
                         Attribute<Matrix>::Ptr leftVector) override;
  /// Stamps system matrix
  void mnaCompApplySystemMatrixStamp(SparseMatrixRow &systemMatrix) override;
  /// Stamps right side (source) vector
  void mnaCompApplyRightSideVectorStamp(Matrix &rightVector) override;
  /// Returns current through the component
  void mnaCompUpdateCurrent(const Matrix &leftVector) override;

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
  // #### Implementation of eigenvalue component interface ####
  void stampBranchNodeIncidenceMatrix(UInt branchIdx,
                                      Matrix &branchNodeIncidenceMatrix) final;
};
} // namespace Ph1
} // namespace EMT
} // namespace CPS
