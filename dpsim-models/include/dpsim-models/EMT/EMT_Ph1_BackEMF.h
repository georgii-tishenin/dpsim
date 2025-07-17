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


  
    class InertiaMoment;  // Forward declaration so that we can get the angular speed (omega) from the inertia moment



/// \brief Back EMF source model
///
/// This model uses modified nodal analysis to represent a back EMF term (psi_F*omega) as an ideal voltage source.
/// For a voltage source between nodes j and k, a new variable (current across the voltage source)
/// is added to the left side vector
/// as unkown and it is taken into account for the equation of node j as positve and for the equation
/// of node k as negative. Moreover
/// a new equation ej - ek = V is added to the problem.
class BackEMF : public MNASimPowerComp<Real>,
                      public SharedFactory<BackEMF>,
                      public EigenvalueCompInterface {
private:
  Real mTimeStep;

protected:
  void updateVoltage(Real time);

public:
  const Attribute<Real>::Ptr mFieldFluxLinkageRef;
  //const Attribute<Real>::Ptr mAngularSpeedRef;
  /// Defines UID, name and logging level
  BackEMF(String uid, String name,
                Logger::Level logLevel = Logger::Level::off);
  ///
  BackEMF(String name, Logger::Level logLevel = Logger::Level::off)
      : BackEMF(name, name, logLevel) {}

  /// Pointer to the inertia moment component from which the angular speed is derived
  std::shared_ptr<CPS::EMT::Ph1::InertiaMoment> mInertiaMoment;

  void setInertiaMoment(const std::shared_ptr<CPS::EMT::Ph1::InertiaMoment>& pt) {
    mInertiaMoment = pt;
  }


  void setParameters(Real field_flux_linkage);

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
