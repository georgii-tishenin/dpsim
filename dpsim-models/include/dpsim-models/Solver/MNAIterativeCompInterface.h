// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include <memory>
#include <vector>

#include <dpsim-models/Definitions.h>

namespace CPS {

/// Describes the changes made by a component during an MNA iteration.
struct MNAIterationUpdate {
  Bool requiresIteration = false;
  Bool systemMatrixChanged = false;
  Bool rightSideVectorChanged = false;
};

/// Interface for components that require repeated MNA solutions within one
/// simulation time step.
class MNAIterativeCompInterface {
public:
  using Ptr = std::shared_ptr<MNAIterativeCompInterface>;
  using List = std::vector<Ptr>;

  virtual ~MNAIterativeCompInterface() = default;

  /// Prepare component-local data that must remain fixed during iteration.
  virtual void mnaInitializeIteration(Real time, Int timeStepCount) = 0;

  /// Update the component linearization from the latest global MNA solution.
  /// The component must update its matrix and right-vector stamps before
  /// returning the corresponding change flags.
  virtual MNAIterationUpdate mnaUpdateIteration(const Matrix &leftVector) = 0;

  /// Finish the iteration without advancing the component state. State
  /// advancement remains part of the normal MNA post-step task.
  virtual void mnaFinalizeIteration() = 0;
};

} // namespace CPS
