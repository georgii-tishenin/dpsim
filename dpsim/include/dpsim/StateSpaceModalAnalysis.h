// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0
#pragma once

#include <dpsim/Definitions.h>
#include <dpsim/MNAStateSpaceExtractor.h>

#include <stdexcept>
#include <vector>

namespace DPsim {

enum class StateSpaceAnalysisFrame {
  Native,
  GlobalDQ0,
};

enum class StateSpacePoleMapping {
  /// Inverse trapezoidal (bilinear) mapping. Use when the complete model was
  /// discretized with the trapezoidal rule and contains no explicit delays.
  Bilinear,

  /// Sampled-data mapping lambda = log(z) / dt. Use for models containing
  /// explicit discrete delays or other non-trapezoidal update operations.
  Logarithmic,
};

/// Performs modal analysis of an extracted discrete-time state-space model.
///
/// The analysis uses the state matrix provided by MNAStateSpaceExtractor and
/// maps discrete-time eigenvalues to continuous-time equivalent eigenvalues
/// using a selectable pole mapping.
class StateSpaceModalAnalysis {
public:
  explicit StateSpaceModalAnalysis(const MNAStateSpaceExtractor &extractor);

  /// Select the coordinate frame used for modal analysis.
  ///
  /// Native uses the extracted state matrix directly.
  /// GlobalDQ0 transforms registered abc state blocks to a global dq0 frame.
  void setAnalysisFrame(StateSpaceAnalysisFrame frame) {
    mAnalysisFrame = frame;
  }

  /// Set the global dq0 frame.
  ///
  /// omega is the constant synchronous angular speed.
  /// theta0 is the global frame angle at t = 0.
  void setGlobalDq0Frame(Real omega, Real theta0 = 0.0) {
    mGlobalOmega = omega;
    mGlobalTheta0 = theta0;
  }

  void setPoleMapping(StateSpacePoleMapping mapping) { mPoleMapping = mapping; }

  /// Eliminate metadata-registered auxiliary embedding coordinates on the
  /// reachable invariant manifold before modal quantities are evaluated.
  ///
  /// The raw extracted matrix is not modified. This is intended for
  /// nonminimal exact-step realizations whose auxiliary coordinates are not
  /// independent physical states and must therefore not receive physical
  /// participation factors.
  void setReduceAuxiliaryStates(Bool reduce) {
    mReduceAuxiliaryStates = reduce;
  }

  /// Exclude metadata-registered zero-sequence coordinates from a balanced
  /// GlobalDQ0 modal analysis.
  ///
  /// The reduction is accepted only when the discarded zero-sequence
  /// coordinates are decoupled from the retained states within the configured
  /// relative tolerance. It must not be used for unbalanced studies where the
  /// zero-sequence subsystem is part of the physical dynamics of interest.
  void setExcludeDecoupledZeroSequenceStates(Bool exclude) {
    mExcludeDecoupledZeroSequenceStates = exclude;
  }

  void setZeroSequenceCouplingTolerance(Real tolerance) {
    if (tolerance < 0.0)
      throw std::invalid_argument(
          "Zero-sequence coupling tolerance must be nonnegative.");
    mZeroSequenceCouplingTolerance = tolerance;
  }

  /// Update modal quantities from the current extracted state matrix.
  void update();

  /// Eigenvalues of the extracted discrete-time state matrix in the selected analysis frame.
  const CPS::VectorComp &getDiscreteEigenvalues() const {
    return mDiscreteEigenvalues;
  }

  /// Continuous-time equivalent eigenvalues reconstructed with the selected
  /// pole mapping.
  const CPS::VectorComp &getContinuousEigenvalues() const {
    return mContinuousEigenvalues;
  }

  /// Right eigenvectors of the selected discrete analysis state matrix.
  ///
  /// Columns correspond to modes in the same order as getDiscreteEigenvalues().
  /// These eigenvectors represent mode shapes in the selected analysis frame.
  const CPS::MatrixComp &getRightEigenvectors() const {
    return mRightEigenvectors;
  }

  /// Left eigenvectors of the selected discrete analysis state matrix.
  ///
  /// Rows correspond to modes and are normalized such that
  /// getLeftEigenvectors() * getRightEigenvectors() = I.
  const CPS::MatrixComp &getLeftEigenvectors() const {
    return mLeftEigenvectors;
  }

  /// Participation factors of the selected discrete analysis state matrix.
  ///
  /// P(state, mode) = rightEigenvectors(state, mode)
  ///                  * leftEigenvectors(mode, state)
  ///
  /// Rows correspond to states, columns correspond to modes. Columns follow the
  /// same mode order as getDiscreteEigenvalues() and getContinuousEigenvalues().
  const CPS::MatrixComp &getParticipationFactors() const {
    return mParticipationFactors;
  }

  /// State names in the selected analysis frame.
  ///
  /// In Native frame, registered abc states keep their native abc names.
  /// In GlobalDQ0 frame, registered abc state blocks are labelled as dq0 states.
  const std::vector<String> &getStateNames() const { return mStateNames; }

  /// Relative invariance residual ||A E - E A_r|| / ||A|| of the most recent
  /// auxiliary-state reduction; zero when no reduction was requested.
  Real getAuxiliaryReductionResidual() const {
    return mAuxiliaryReductionResidual;
  }

  /// Maximum nearest-neighbour difference between each reduced pole and the
  /// unreduced augmented spectrum in the most recent update.
  Real getAuxiliaryReductionPoleError() const {
    return mAuxiliaryReductionPoleError;
  }

  /// Relative off-diagonal coupling of the discarded zero-sequence block in
  /// the most recent update; zero when no such reduction was requested.
  Real getZeroSequenceCouplingResidual() const {
    return mZeroSequenceCouplingResidual;
  }

private:
  Matrix buildDiscreteStateMatrixInAnalysisFrame() const;

  Matrix reduceToReachablePhysicalStateMatrix(const Matrix &matrix,
                                               std::vector<UInt> &physicalIndices);

  Matrix excludeDecoupledZeroSequenceStates(
      const Matrix &matrix, std::vector<UInt> &originalStateIndices);

  Matrix buildGlobalDq0Transformation(Real theta) const;

  std::vector<String> buildStateNamesInAnalysisFrame() const;

  Complex mapDiscreteToContinuous(const Complex &z) const;

  const MNAStateSpaceExtractor &mExtractor;

  StateSpaceAnalysisFrame mAnalysisFrame = StateSpaceAnalysisFrame::Native;

  StateSpacePoleMapping mPoleMapping = StateSpacePoleMapping::Bilinear;

  Bool mReduceAuxiliaryStates = false;

  Bool mExcludeDecoupledZeroSequenceStates = false;

  Real mZeroSequenceCouplingTolerance = 1e-8;

  Real mZeroSequenceCouplingResidual = 0.0;

  Real mAuxiliaryReductionResidual = 0.0;

  Real mAuxiliaryReductionPoleError = 0.0;

  Real mGlobalOmega = 0.0;

  Real mGlobalTheta0 = 0.0;

  CPS::VectorComp mDiscreteEigenvalues;

  CPS::VectorComp mContinuousEigenvalues;

  CPS::MatrixComp mRightEigenvectors;
  CPS::MatrixComp mLeftEigenvectors;
  CPS::MatrixComp mParticipationFactors;

  std::vector<String> mStateNames;
};

} // namespace DPsim
