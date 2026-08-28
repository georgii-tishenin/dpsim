// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <Eigen/QR>

#include <dpsim-models/MathUtils.h>
#include <dpsim/StateSpaceModalAnalysis.h>

#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>

namespace DPsim {

namespace {

Matrix parkTransformDQ0(Real theta) {
  Matrix transform(3, 3);

  const Real k = std::sqrt(2.0 / 3.0);
  const Real k0 = 1.0 / std::sqrt(3.0);

  transform.row(0) << k * std::cos(theta), k * std::cos(theta - 2.0 * PI / 3.0),
      k * std::cos(theta + 2.0 * PI / 3.0);

  transform.row(1) << -k * std::sin(theta),
      -k * std::sin(theta - 2.0 * PI / 3.0),
      -k * std::sin(theta + 2.0 * PI / 3.0);

  transform.row(2) << k0, k0, k0;

  return transform;
}

String fallbackStateName(UInt index) { return "x" + std::to_string(index); }

} // namespace

StateSpaceModalAnalysis::StateSpaceModalAnalysis(
    const MNAStateSpaceExtractor &extractor)
    : mExtractor(extractor) {}

void StateSpaceModalAnalysis::update() {
  if (!mExtractor.isInitialized())
    throw std::logic_error("StateSpaceModalAnalysis requires an initialized "
                           "MNAStateSpaceExtractor.");

  Matrix Ad = buildDiscreteStateMatrixInAnalysisFrame();

  if (Ad.rows() == 0) {
    mDiscreteEigenvalues.resize(0);
    mContinuousEigenvalues.resize(0);

    mRightEigenvectors.resize(0, 0);
    mLeftEigenvectors.resize(0, 0);
    mParticipationFactors.resize(0, 0);

    mStateNames.clear();

    return;
  }

  if (Ad.rows() != Ad.cols())
    throw std::logic_error(
        "StateSpaceModalAnalysis requires a square state matrix.");

  mStateNames = buildStateNamesInAnalysisFrame();
  mAuxiliaryReductionResidual = 0.0;
  mAuxiliaryReductionPoleError = 0.0;
  mZeroSequenceCouplingResidual = 0.0;

  std::vector<UInt> originalStateIndices(static_cast<UInt>(Ad.rows()));
  std::iota(originalStateIndices.begin(), originalStateIndices.end(), 0U);
  if (mReduceAuxiliaryStates &&
      !mExtractor.getMetadata().auxiliaryStateIndices.empty()) {
    std::vector<UInt> physicalIndices;
    Ad = reduceToReachablePhysicalStateMatrix(Ad, physicalIndices);
    std::vector<String> physicalStateNames;
    physicalStateNames.reserve(physicalIndices.size());
    for (const UInt idx : physicalIndices)
      physicalStateNames.push_back(mStateNames[idx]);
    mStateNames = std::move(physicalStateNames);
    originalStateIndices = std::move(physicalIndices);
  }

  if (mExcludeDecoupledZeroSequenceStates) {
    Ad = excludeDecoupledZeroSequenceStates(Ad, originalStateIndices);
  }

  Eigen::EigenSolver<Matrix> eigenSolver(Ad, true);

  if (eigenSolver.info() != Eigen::Success)
    throw std::runtime_error(
        "StateSpaceModalAnalysis: eigenvalue computation failed.");

  mDiscreteEigenvalues = eigenSolver.eigenvalues();

  mContinuousEigenvalues.resize(mDiscreteEigenvalues.rows());

  for (Eigen::Index idx = 0; idx < mDiscreteEigenvalues.rows(); ++idx)
    mContinuousEigenvalues(idx) =
        mapDiscreteToContinuous(mDiscreteEigenvalues(idx));

  mRightEigenvectors = eigenSolver.eigenvectors();

  Eigen::FullPivLU<CPS::MatrixComp> eigenvectorLu(mRightEigenvectors);

  if (!eigenvectorLu.isInvertible())
    throw std::runtime_error(
        "StateSpaceModalAnalysis: cannot compute participation factors because "
        "the eigenvector matrix is singular.");

  mLeftEigenvectors = eigenvectorLu.inverse();

  mParticipationFactors = CPS::Math::elementwiseProduct(
      mRightEigenvectors, mLeftEigenvectors.transpose());
}

Matrix StateSpaceModalAnalysis::reduceToReachablePhysicalStateMatrix(
    const Matrix &matrix, std::vector<UInt> &physicalIndices) {
  const UInt fullStateCount = static_cast<UInt>(matrix.rows());
  const auto &auxiliaryIndices =
      mExtractor.getMetadata().auxiliaryStateIndices;
  std::vector<Bool> isAuxiliary(fullStateCount, false);
  for (const UInt idx : auxiliaryIndices) {
    if (idx >= fullStateCount)
      throw std::logic_error(
          "Auxiliary state index lies outside the modal state matrix.");
    if (isAuxiliary[idx])
      throw std::logic_error("Duplicate auxiliary state index.");
    isAuxiliary[idx] = true;
  }

  physicalIndices.clear();
  physicalIndices.reserve(fullStateCount - auxiliaryIndices.size());
  for (UInt idx = 0; idx < fullStateCount; ++idx) {
    if (!isAuxiliary[idx])
      physicalIndices.push_back(idx);
  }
  if (physicalIndices.empty())
    throw std::logic_error(
        "Auxiliary-state reduction requires at least one physical state.");

  const UInt physicalStateCount = physicalIndices.size();
  const UInt auxiliaryStateCount = auxiliaryIndices.size();
  Matrix physicalAdvance(physicalStateCount, fullStateCount);
  Matrix auxiliaryAdvance(auxiliaryStateCount, fullStateCount);
  for (UInt row = 0; row < physicalStateCount; ++row)
    physicalAdvance.row(row) = matrix.row(physicalIndices[row]);
  for (UInt row = 0; row < auxiliaryStateCount; ++row)
    auxiliaryAdvance.row(row) = matrix.row(auxiliaryIndices[row]);

  Eigen::FullPivLU<Matrix> physicalAdvanceLu(physicalAdvance);
  if (static_cast<UInt>(physicalAdvanceLu.rank()) != physicalStateCount) {
    throw std::runtime_error(
        "Cannot reduce auxiliary states because the retained physical "
        "coordinates do not parameterize the reachable subspace.");
  }

  // A reached augmented state lies on u_aux = H x_physical. Solve
  // A_aux = H A_physical, then express the invariant dynamics through the
  // embedding E x_physical = [x_physical, H x_physical].
  const Matrix graph =
      physicalAdvance.transpose()
          .colPivHouseholderQr()
          .solve(auxiliaryAdvance.transpose())
          .transpose();
  Matrix embedding = Matrix::Zero(fullStateCount, physicalStateCount);
  for (UInt col = 0; col < physicalStateCount; ++col)
    embedding(physicalIndices[col], col) = 1.0;
  for (UInt row = 0; row < auxiliaryStateCount; ++row)
    embedding.row(auxiliaryIndices[row]) = graph.row(row);

  const Matrix reduced = physicalAdvance * embedding;
  mAuxiliaryReductionResidual =
      (matrix * embedding - embedding * reduced).norm() /
      std::max(matrix.norm(), std::numeric_limits<Real>::epsilon());
  if (mAuxiliaryReductionResidual > 1e-8) {
    throw std::runtime_error(
        "Auxiliary-state reduction failed the invariant-subspace check.");
  }

  Eigen::EigenSolver<Matrix> augmentedSolver(matrix, false);
  Eigen::EigenSolver<Matrix> reducedSolver(reduced, false);
  if (augmentedSolver.info() != Eigen::Success ||
      reducedSolver.info() != Eigen::Success) {
    throw std::runtime_error(
        "Auxiliary-state reduction pole-preservation check failed.");
  }
  const CPS::VectorComp augmentedPoles = augmentedSolver.eigenvalues();
  const CPS::VectorComp reducedPoles = reducedSolver.eigenvalues();
  for (Eigen::Index reducedIdx = 0; reducedIdx < reducedPoles.rows();
       ++reducedIdx) {
    Real nearest = std::numeric_limits<Real>::infinity();
    for (Eigen::Index augmentedIdx = 0;
         augmentedIdx < augmentedPoles.rows(); ++augmentedIdx) {
      nearest = std::min(
          nearest,
          std::abs(reducedPoles(reducedIdx) - augmentedPoles(augmentedIdx)));
    }
    mAuxiliaryReductionPoleError =
        std::max(mAuxiliaryReductionPoleError, nearest);
  }
  return reduced;
}

Matrix StateSpaceModalAnalysis::excludeDecoupledZeroSequenceStates(
    const Matrix &matrix, std::vector<UInt> &originalStateIndices) {
  if (mAnalysisFrame != StateSpaceAnalysisFrame::GlobalDQ0) {
    throw std::logic_error(
        "Zero-sequence exclusion requires the GlobalDQ0 analysis frame.");
  }
  if (originalStateIndices.size() != static_cast<std::size_t>(matrix.rows())) {
    throw std::logic_error(
        "Modal state-index mapping is inconsistent with the state matrix.");
  }

  const UInt fullStateCount = mExtractor.getStateCount();
  std::vector<Bool> isZeroSequence(fullStateCount, false);
  for (const auto &abcBlock : mExtractor.getMetadata().abcStateBlocks) {
    const UInt zeroIndex = abcBlock.indices[2];
    if (zeroIndex >= fullStateCount)
      throw std::logic_error(
          "Zero-sequence state index lies outside the modal state matrix.");
    isZeroSequence[zeroIndex] = true;
  }

  std::vector<UInt> retainedPositions;
  std::vector<UInt> discardedPositions;
  retainedPositions.reserve(originalStateIndices.size());
  discardedPositions.reserve(originalStateIndices.size());
  for (UInt position = 0; position < originalStateIndices.size(); ++position) {
    const UInt originalIndex = originalStateIndices[position];
    if (originalIndex >= fullStateCount)
      throw std::logic_error(
          "Modal state-index mapping lies outside the extracted state set.");
    (isZeroSequence[originalIndex] ? discardedPositions : retainedPositions)
        .push_back(position);
  }

  if (discardedPositions.empty())
    return matrix;
  if (retainedPositions.empty())
    throw std::logic_error(
        "Zero-sequence exclusion would discard every modal state.");

  Real couplingSquared = 0.0;
  for (const UInt retained : retainedPositions) {
    for (const UInt discarded : discardedPositions) {
      const Real retainedToDiscarded = matrix(discarded, retained);
      const Real discardedToRetained = matrix(retained, discarded);
      couplingSquared += retainedToDiscarded * retainedToDiscarded;
      couplingSquared += discardedToRetained * discardedToRetained;
    }
  }
  mZeroSequenceCouplingResidual =
      std::sqrt(couplingSquared) /
      std::max(matrix.norm(), std::numeric_limits<Real>::epsilon());
  if (mZeroSequenceCouplingResidual > mZeroSequenceCouplingTolerance) {
    throw std::runtime_error(
        "Cannot exclude zero-sequence states because they are coupled to the "
        "retained modal subsystem.");
  }

  const UInt retainedCount = retainedPositions.size();
  Matrix reduced(retainedCount, retainedCount);
  std::vector<String> retainedNames;
  std::vector<UInt> retainedOriginalIndices;
  retainedNames.reserve(retainedCount);
  retainedOriginalIndices.reserve(retainedCount);
  for (UInt row = 0; row < retainedCount; ++row) {
    retainedNames.push_back(mStateNames[retainedPositions[row]]);
    retainedOriginalIndices.push_back(originalStateIndices[retainedPositions[row]]);
    for (UInt col = 0; col < retainedCount; ++col) {
      reduced(row, col) =
          matrix(retainedPositions[row], retainedPositions[col]);
    }
  }
  mStateNames = std::move(retainedNames);
  originalStateIndices = std::move(retainedOriginalIndices);
  return reduced;
}

Matrix
StateSpaceModalAnalysis::buildDiscreteStateMatrixInAnalysisFrame() const {
  const Matrix &nativeAd = mExtractor.getDiscreteStateMatrix();

  if (mAnalysisFrame == StateSpaceAnalysisFrame::Native)
    return nativeAd;

  if (mAnalysisFrame == StateSpaceAnalysisFrame::GlobalDQ0) {
    if (!mExtractor.hasExtractionTime()) {
      throw std::logic_error(
          "GlobalDQ0 modal analysis requires a valid extraction timestamp.");
    }

    if (mGlobalOmega <= 0.0) {
      throw std::logic_error(
          "GlobalDQ0 modal analysis requires a positive frame angular speed.");
    }

    const Real time = mExtractor.getLastExtractionTime();
    const Real timeStep = mExtractor.getTimeStep();

    const Real thetaNow = mGlobalTheta0 + mGlobalOmega * time;
    const Real thetaNext = thetaNow + mGlobalOmega * timeStep;

    const Matrix transformNow = buildGlobalDq0Transformation(thetaNow);
    const Matrix transformNext = buildGlobalDq0Transformation(thetaNext);

    // For a time-dependent discrete coordinate transformation
    // xGlobalDq0[k] = T[k] xNative[k], the transformed transition matrix is
    // AdGlobalDq0[k] = T[k+1] AdNative[k] T[k]^{-1}.
    //
    // The Park transform is power-invariant, so T^{-1} = T^T.
    return transformNext * nativeAd * transformNow.transpose();
  }

  throw std::logic_error("Unsupported state-space analysis frame.");
}

Matrix StateSpaceModalAnalysis::buildGlobalDq0Transformation(Real theta) const {
  const UInt stateCount = mExtractor.getStateCount();

  Matrix transform = Matrix::Identity(stateCount, stateCount);

  const Matrix park = parkTransformDQ0(theta);

  for (const auto &abcBlock : mExtractor.getMetadata().abcStateBlocks) {
    for (UInt row = 0; row < 3; ++row) {
      for (UInt col = 0; col < 3; ++col) {
        transform(abcBlock.indices[row], abcBlock.indices[col]) =
            park(row, col);
      }
    }
  }

  return transform;
}

std::vector<String>
StateSpaceModalAnalysis::buildStateNamesInAnalysisFrame() const {
  const UInt stateCount = mExtractor.getStateCount();
  const auto &metadata = mExtractor.getMetadata();

  std::vector<String> stateNames(stateCount);

  for (UInt idx = 0; idx < stateCount; ++idx) {
    if (idx < metadata.stateNames.size() && !metadata.stateNames[idx].empty())
      stateNames[idx] = metadata.stateNames[idx];
    else
      stateNames[idx] = fallbackStateName(idx);
  }

  if (mAnalysisFrame == StateSpaceAnalysisFrame::Native)
    return stateNames;

  if (mAnalysisFrame == StateSpaceAnalysisFrame::GlobalDQ0) {
    for (const auto &abcBlock : metadata.abcStateBlocks) {
      if (abcBlock.name.empty())
        throw std::logic_error(
            "GlobalDQ0 modal analysis requires named abc state blocks.");

      stateNames[abcBlock.indices[0]] = abcBlock.name + "_d";
      stateNames[abcBlock.indices[1]] = abcBlock.name + "_q";
      stateNames[abcBlock.indices[2]] = abcBlock.name + "_0";
    }

    return stateNames;
  }

  throw std::logic_error("Unsupported state-space analysis frame.");
}

CPS::Complex
StateSpaceModalAnalysis::mapDiscreteToContinuous(const CPS::Complex &z) const {
  if (mPoleMapping == StateSpacePoleMapping::Logarithmic)
    return std::log(z) / mExtractor.getTimeStep();

  const CPS::Complex one(1.0, 0.0);
  const CPS::Complex denominator = z + one;

  if (std::abs(denominator) <= DOUBLE_EPSILON)
    return CPS::Complex(std::numeric_limits<Real>::infinity(), 0.0);

  return (2.0 / mExtractor.getTimeStep()) * (z - one) / denominator;
}

} // namespace DPsim
