/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#pragma once

#include <dpsim-models/Logger.h>
#include <dpsim-models/SimSignalComp.h>

namespace CPS {
namespace Signal {

/// Turbine governor model to be used with synchronous generator
class TurbineGovernor : public SimSignalComp,
                        public SharedFactory<TurbineGovernor> {
protected:
  // ### Steam Turbine Parameters ####
  /// Time constant of main inlet volume and steam chest
  Real mTa;
  /// Time constant of reheater
  Real mTb;
  /// Time constant of cross over piping and LP inlet volumes
  Real mTc;
  /// Fraction of total turbine power generated by HP section
  Real mFa;
  /// Fraction of total turbine power generated by IP section
  Real mFb;
  /// Fraction of total turbine power generated by LP section
  Real mFc;

  // ### Regulator ####
  /// Main droop
  Real mK;

  // ### Speed Relay and servo motor ####
  /// Time constant of speed relay
  Real mTsr;
  /// Time constant of servo motor
  Real mTsm;
  /// Opening rate limit
  Real mLc1 = 0.1;
  /// Closing rate limit
  Real mLc2 = -1;

  // ### Vaariables ###
  /// Mechanical Torque in pu (output of steam turbine)
  Real mTm = 0;
  /// Valve position
  Real mVcv = 0;
  /// Valve position changing rate
  Real mpVcv = 0;
  /// Boiler pressure
  Real mPb = 1;
  /// Speed reference
  Real mOmRef = 1;
  /// Speed
  Real mOm = 0;
  /// Power Reference
  Real mPmRef;
  /// Input of speed realy
  Real Psr_in = 0;
  /// Input of servor motor
  Real Psm_in = 0;

  Real AuxVar = 0;

  Real T1;
  Real T2;

public:
  ///
  TurbineGovernor(String name) : SimSignalComp(name, name) {}

  /// Constructor with log level
  TurbineGovernor(String name, CPS::Logger::Level logLevel)
      : SimSignalComp(name, name, logLevel) {}

  /// Initializes exciter parameters
  void setParameters(Real Ta, Real Tb, Real Tc, Real Fa, Real Fb, Real Fc,
                     Real K, Real Tsr, Real Tsm);

  ///
  void initialize(Real PmRef, Real Tm_init);

  /// Performs an step to update field voltage value
  Real step(Real mOm, Real mOmRef, Real PmRef, Real dt);
};
} // namespace Signal
} // namespace CPS
