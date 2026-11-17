/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim/MNASolverDirect.h>

#include <dpsim-models/MathUtils.h>
#include <dpsim/SequentialScheduler.h>

#include <algorithm>
#include <fstream>
#include <numeric>

using namespace DPsim;
using namespace CPS;

namespace DPsim {

template <typename VarType>
MnaSolverDirect<VarType>::MnaSolverDirect(String name, CPS::Domain domain,
                                          CPS::Logger::Level logLevel)
    : MnaSolver<VarType>(name, domain, logLevel) {
#ifdef WITH_KLU
  mImplementationInUse = DirectLinearSolverImpl::KLU;
#elif defined(WITH_SPARSE)
  mImplementationInUse = DirectLinearSolverImpl::SparseLU;
#else
  mImplementationInUse = DirectLinearSolverImpl::DenseLU;
#endif
}

template <typename VarType>
void MnaSolverDirect<VarType>::switchedMatrixEmpty(std::size_t index) {
  mSwitchedMatrices[std::bitset<SWITCH_NUM>(index)][0].setZero();
}

template <typename VarType>
void MnaSolverDirect<VarType>::switchedMatrixEmpty(std::size_t swIdx,
                                                   Int freqIdx) {
  mSwitchedMatrices[std::bitset<SWITCH_NUM>(swIdx)][freqIdx].setZero();
}

template <typename VarType>
void MnaSolverDirect<VarType>::switchedMatrixStamp(
    std::size_t index, std::vector<std::shared_ptr<CPS::MNAInterface>> &comp) {
  auto bit = std::bitset<SWITCH_NUM>(index);
  auto &sys = mSwitchedMatrices[bit][0];
  for (auto component : comp) {
    component->mnaApplySystemMatrixStamp(sys);
  }
  for (UInt i = 0; i < mSwitches.size(); ++i)
    mSwitches[i]->mnaApplySwitchSystemMatrixStamp(bit[i], sys, 0);

  // Compute LU-factorization for system matrix
  sys.prune([](Eigen::Index, Eigen::Index, Real value) { return value != 0.; });
  sys.makeCompressed();
  mDirectLinearSolvers[bit][0]->setConfiguration(mConfigurationInUse);
  mDirectLinearSolvers[bit][0]->preprocessing(sys,
                                              mListVariableSystemMatrixEntries);
  auto start = std::chrono::steady_clock::now();
  mDirectLinearSolvers[bit][0]->factorize(sys);
  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<Real> diff = end - start;
  mFactorizeTimes.push_back(diff.count());
}

template <typename VarType>
void MnaSolverDirect<VarType>::switchedMatrixStamp(
    std::size_t swIdx, Int freqIdx, CPS::MNAInterface::List &components,
    CPS::MNASwitchInterface::List &switches) {
  const auto bit = std::bitset<SWITCH_NUM>(swIdx);
  auto &systemMatrix = mSwitchedMatrices[bit][freqIdx];
  for (auto component : components)
    component->mnaApplySystemMatrixStampHarm(systemMatrix, freqIdx);
  systemMatrix.makeCompressed();
  for (UInt i = 0; i < switches.size(); ++i)
    switches[i]->mnaApplySwitchSystemMatrixStamp(bit[i], systemMatrix, freqIdx);

  auto &directSolver = mDirectLinearSolvers[bit][freqIdx];
  systemMatrix.prune(
      [](Eigen::Index, Eigen::Index, Real value) { return value != 0.; });
  systemMatrix.makeCompressed();
  directSolver->setConfiguration(mConfigurationInUse);
  directSolver->preprocessing(systemMatrix, mListVariableSystemMatrixEntries);
  const auto start = std::chrono::steady_clock::now();
  directSolver->factorize(systemMatrix);
  const auto end = std::chrono::steady_clock::now();
  mFactorizeTimes.push_back(std::chrono::duration<Real>(end - start).count());
}

template <typename VarType>
void MnaSolverDirect<VarType>::stampVariableSystemMatrix() {

  this->mDirectLinearSolverVariableSystemMatrix =
      createDirectSolverImplementation(mSLog);
  this->mDirectLinearSolverVariableSystemMatrix->setConfiguration(
      mConfigurationInUse);
  if (auto dense = std::dynamic_pointer_cast<DenseLUAdapter>(
          mDirectLinearSolverVariableSystemMatrix))
    dense->setDetailedProfiling(Solver::mLogSolveTimes);

  SPDLOG_LOGGER_INFO(mSLog,
                     "Number of variable Elements: {}"
                     "\nNumber of MNA components: {}",
                     mVariableComps.size(), mMNAComponents.size());

  // Build base matrix with only static elements
  mBaseSystemMatrix.setZero();
  CPS::Math::clearSparseMatrixValueCache(mVariableSystemMatrix);
  for (auto statElem : mMNAComponents)
    statElem->mnaApplySystemMatrixStamp(mBaseSystemMatrix);
  // Common exact-zero policy: static zeros are discarded once, before any
  // value-index cache exists. No threshold or runtime pruning is used.
  mBaseSystemMatrix.prune(
      [](Eigen::Index, Eigen::Index, Real value) { return value != 0.; });
  mBaseSystemMatrix.makeCompressed();
  if (Solver::mLogLevel <= CPS::Logger::Level::info)
    SPDLOG_LOGGER_INFO(mSLog, "Base matrix with only static elements: {}",
                       Logger::matrixToString(mBaseSystemMatrix));
  mSLog->flush();

  // Continue from base matrix
  mVariableSystemMatrix = mBaseSystemMatrix;
  std::sort(mListVariableSystemMatrixEntries.begin(),
            mListVariableSystemMatrixEntries.end());
  mListVariableSystemMatrixEntries.erase(
      std::unique(mListVariableSystemMatrixEntries.begin(),
                  mListVariableSystemMatrixEntries.end()),
      mListVariableSystemMatrixEntries.end());
  for (const auto &[row, column] : mListVariableSystemMatrixEntries)
    mVariableSystemMatrix.coeffRef(row, column) += 0.;

  // Now stamp initial state of variable elements and switches into matrix
  SPDLOG_LOGGER_INFO(mSLog, "Stamping variable elements");
  for (auto varElem : mMNAIntfVariableComps)
    varElem->mnaApplySystemMatrixStamp(mVariableSystemMatrix);

  mVariableSystemMatrix.prune(
      [this](Eigen::Index row, Eigen::Index column, Real value) {
        return value != 0. || mBaseSystemMatrix.coeff(row, column) != 0. ||
               std::binary_search(mListVariableSystemMatrixEntries.begin(),
                                  mListVariableSystemMatrixEntries.end(),
                                  std::make_pair(static_cast<UInt>(row),
                                                 static_cast<UInt>(column)));
      });
  mVariableSystemMatrix.makeCompressed();

  // Keep the reset matrix on exactly the same sparse pattern and in the same
  // coefficient order as the variable matrix. Variable-only entries must be
  // retained as explicit zeros here: recomputeSystemMatrix() copies values by
  // array position, and the static-only pattern can have a different ordering
  // once variable stamps have inserted additional entries.
  SparseMatrix resetMatrix = mVariableSystemMatrix;
  for (Eigen::Index row = 0; row < resetMatrix.outerSize(); ++row)
    for (SparseMatrix::InnerIterator it(resetMatrix, row); it; ++it)
      it.valueRef() = mBaseSystemMatrix.coeff(it.row(), it.col());
  mBaseSystemMatrix = std::move(resetMatrix);

  mVariableValueIndices.clear();
  mFactorizedVariableValues.clear();
  for (const auto &[row, column] : mListVariableSystemMatrixEntries) {
    const auto index = &mVariableSystemMatrix.coeffRef(row, column) -
                       mVariableSystemMatrix.valuePtr();
    mVariableValueIndices.push_back(index);
    mFactorizedVariableValues.push_back(
        mVariableSystemMatrix.valuePtr()[index]);
  }
  // If any variable component omits its entries, compare the full pattern.
  // Other components' declarations do not cover that component's changes.
  const bool incompleteDeclarations =
      std::any_of(mVariableComps.begin(), mVariableComps.end(),
                  [](const auto &component) {
                    return component->mVariableSystemMatrixEntries.empty();
                  });
  if (incompleteDeclarations) {
    mVariableValueIndices.clear();
    mFactorizedVariableValues.clear();
    for (Eigen::Index index = 0; index < mVariableSystemMatrix.nonZeros();
         ++index) {
      mVariableValueIndices.push_back(index);
      mFactorizedVariableValues.push_back(
          mVariableSystemMatrix.valuePtr()[index]);
    }
  }

  if (Solver::mLogLevel <= CPS::Logger::Level::info)
    SPDLOG_LOGGER_INFO(mSLog, "Initial system matrix with variable elements {}",
                       Logger::matrixToString(mVariableSystemMatrix));
  /* TODO: find replacement for flush() */
  mSLog->flush();

  // Calculate factorization of current matrix
  mDirectLinearSolverVariableSystemMatrix->preprocessing(
      mVariableSystemMatrix, mListVariableSystemMatrixEntries);
  CPS::Math::prepareSparseMatrixValueCache(mVariableSystemMatrix);

  auto start = std::chrono::steady_clock::now();
  mDirectLinearSolverVariableSystemMatrix->factorize(mVariableSystemMatrix);
  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<Real> diff = end - start;
  mFactorizeTimes.push_back(diff.count());
}

template <typename VarType>
void MnaSolverDirect<VarType>::solveWithSystemMatrixRecomputation(
    Real time, Int timeStepCount) {
  using Clock = std::chrono::steady_clock;
  const bool profile = Solver::mLogSolveTimes || Solver::mSolverPhaseProfiling;
  auto profileStart = profile ? Clock::now() : Clock::time_point{};
  // Reset source vector
  mRightSideVector.setZero();
  if (profile)
    recordProfileMeasurement("rhs_zero", profileStart, Clock::now());

  // Add together the right side vector (computed by the components'
  // pre-step tasks)
  profileStart = profile ? Clock::now() : Clock::time_point{};
  addRightVectorStamps(mRightSideVector);
  if (profile)
    recordProfileMeasurement("rhs_assembly", profileStart, Clock::now());

  // Get switch and variable comp status and update system matrix and lu factorization accordingly
  profileStart = profile ? Clock::now() : Clock::time_point{};
  mVariableComponentChanged = hasVariableComponentChanged();
  if (profile)
    recordProfileMeasurement("variable_change_detection", profileStart,
                             Clock::now());
  if (mVariableComponentChanged)
    recomputeSystemMatrix(time);
  else
    ++mUnchangedFactorReuseCount;

  // Calculate new solution vector
  const auto start = profile ? Clock::now() : Clock::time_point{};
  mDirectLinearSolverVariableSystemMatrix->solveInto(mRightSideVector,
                                                     *mLeftSideVectorValue);
  if (profile) {
    const auto end = Clock::now();
    mSolveTimes.push_back(std::chrono::duration<Real>(end - start).count());
    recordProfileMeasurement("numeric_solve", start, end);
  }

  // TODO split into separate task? (dependent on x, updating all v attributes)
  profileStart = profile ? Clock::now() : Clock::time_point{};
  for (UInt nodeIdx = 0; nodeIdx < mNumNetNodes; ++nodeIdx)
    mNodes[nodeIdx]->mnaUpdateVoltage(*mLeftSideVectorValue);
  if (profile)
    recordProfileMeasurement("node_voltage_propagation", profileStart,
                             Clock::now());

  // Components' states will be updated by the post-step tasks
}

template <typename VarType>
void MnaSolverDirect<VarType>::recomputeSystemMatrix(Real time) {
  using Clock = std::chrono::steady_clock;
  const bool profile = Solver::mLogSolveTimes || Solver::mSolverPhaseProfiling;
  const auto profileStart = profile ? Clock::now() : Clock::time_point{};
  // The sparse structure is fixed after initialization. Reset only the
  // coefficient values instead of copying/reallocating the full sparse
  // matrix on every timestep.
  std::copy_n(mBaseSystemMatrix.valuePtr(), mBaseSystemMatrix.nonZeros(),
              mVariableSystemMatrix.valuePtr());

  // Now stamp variable elements and switches into matrix
  for (auto comp : mMNAIntfVariableComps)
    comp->mnaApplySystemMatrixStamp(mVariableSystemMatrix);
  ++mMatrixStampCount;
  if (mVariableSystemMatrix.nonZeros() != mBaseSystemMatrix.nonZeros())
    throw std::runtime_error(
        "Variable component changed the finalized sparse pattern");
  if (profile)
    recordProfileMeasurement("matrix_reset_and_variable_stamping", profileStart,
                             Clock::now());

  // Components may advertise themselves as always dirty. Reuse the factors
  // when all declared variable coefficients are exactly unchanged.
  const auto compareStart = profile ? Clock::now() : Clock::time_point{};
  mVariableComponentChanged = false;
  for (std::size_t i = 0; i < mVariableValueIndices.size(); ++i) {
    const Real value =
        mVariableSystemMatrix.valuePtr()[mVariableValueIndices[i]];
    mVariableComponentChanged |= value != mFactorizedVariableValues[i];
    mFactorizedVariableValues[i] = value;
  }
  if (profile)
    recordProfileMeasurement("matrix_value_change_detection", compareStart,
                             Clock::now());
  if (!mVariableComponentChanged) {
    ++mUnchangedFactorReuseCount;
    return;
  }

  // Refactorization of matrix assuming that structure remained
  // constant by omitting analyzePattern
  const auto start = profile ? Clock::now() : Clock::time_point{};
  const bool partial =
      mConfigurationInUse.getPartialRefactorizationMethod() !=
      PARTIAL_REFACTORIZATION_METHOD::NO_PARTIAL_REFACTORIZATION;
  if (partial)
    mDirectLinearSolverVariableSystemMatrix->partialRefactorize(
        mVariableSystemMatrix, mListVariableSystemMatrixEntries);
  else
    mDirectLinearSolverVariableSystemMatrix->refactorize(mVariableSystemMatrix);
  if (profile) {
    const auto end = Clock::now();
    mRecomputationTimes.push_back(
        std::chrono::duration<Real>(end - start).count());
    recordProfileMeasurement(partial ? "partial_refactorization"
                                     : "numeric_refactorization",
                             start, end);
  }
  ++mNumRecomputations;
}

template <typename VarType>
void MnaSolverDirect<VarType>::extractStateSpace(Real time) {
  if (!mStateSpaceExtractor)
    return;

  if (mSystemMatrixRecomputationEnabled) {
    mStateSpaceExtractor->extract(*mDirectLinearSolverVariableSystemMatrix,
                                  mVariableComponentChanged,
                                  mVariableComponentChanged, time);
  } else {
    mStateSpaceExtractor->extract(
        *mDirectLinearSolvers[mCurrentSwitchStatus][0], false,
        mSystemMatrixChanged, time);
  }
}

template <> void MnaSolverDirect<Real>::createEmptySystemMatrix() {
  if (mSwitches.size() > SWITCH_NUM)
    throw SystemError("Too many Switches.");

  if (mSystemMatrixRecomputationEnabled) {
    mBaseSystemMatrix =
        SparseMatrix(mNumMatrixNodeIndices, mNumMatrixNodeIndices);
    mVariableSystemMatrix =
        SparseMatrix(mNumMatrixNodeIndices, mNumMatrixNodeIndices);
  } else {
    for (std::size_t i = 0; i < (1ULL << mSwitches.size()); i++) {
      auto bit = std::bitset<SWITCH_NUM>(i);
      mSwitchedMatrices[bit].push_back(
          SparseMatrix(mNumMatrixNodeIndices, mNumMatrixNodeIndices));
      mDirectLinearSolvers[bit].push_back(
          createDirectSolverImplementation(mSLog));
    }
  }
}

template <> void MnaSolverDirect<Complex>::createEmptySystemMatrix() {
  if (mSwitches.size() > SWITCH_NUM)
    throw SystemError("Too many Switches.");

  if (mFrequencyParallel) {
    for (UInt i = 0; i < std::pow(2, mSwitches.size()); ++i) {
      for (Int freq = 0; freq < mSystem.mFrequencies.size(); ++freq) {
        auto bit = std::bitset<SWITCH_NUM>(i);
        mSwitchedMatrices[bit].push_back(SparseMatrix(
            2 * (mNumMatrixNodeIndices), 2 * (mNumMatrixNodeIndices)));
        mDirectLinearSolvers[bit].push_back(
            createDirectSolverImplementation(mSLog));
      }
    }
  } else if (mSystemMatrixRecomputationEnabled) {
    mBaseSystemMatrix =
        SparseMatrix(2 * (mNumMatrixNodeIndices), 2 * (mNumMatrixNodeIndices));
    mVariableSystemMatrix =
        SparseMatrix(2 * (mNumMatrixNodeIndices), 2 * (mNumMatrixNodeIndices));
  } else {
    for (std::size_t i = 0; i < (1ULL << mSwitches.size()); i++) {
      auto bit = std::bitset<SWITCH_NUM>(i);
      mSwitchedMatrices[bit].push_back(SparseMatrix(
          2 * (mNumTotalMatrixNodeIndices), 2 * (mNumTotalMatrixNodeIndices)));
      mDirectLinearSolvers[bit].push_back(
          createDirectSolverImplementation(mSLog));
    }
  }
}

template <typename VarType>
std::shared_ptr<CPS::Task> MnaSolverDirect<VarType>::createSolveTask() {
  prepareCompactRightVectorStamps();
  return std::make_shared<MnaSolverDirect<VarType>::SolveTask>(*this);
}

template <typename VarType>
std::shared_ptr<CPS::Task> MnaSolverDirect<VarType>::createSolveTaskRecomp() {
  prepareCompactRightVectorStamps();
  return std::make_shared<MnaSolverDirect<VarType>::SolveTaskRecomp>(*this);
}

template <typename VarType>
void MnaSolverDirect<VarType>::prepareCompactRightVectorStamps() {
  std::unordered_map<const Matrix *, CPS::MNAInterface *> compactStamps;
  auto collect = [this, &compactStamps](const auto &components) {
    for (const auto &component : components) {
      component->mnaSetDetailedProfiling(Solver::mLogSolveTimes);
      if (!component->mnaHasCompactRightSideVectorStamp())
        continue;
      auto stamp = component->getRightVector();
      if (stamp->get().size() != 0)
        compactStamps.emplace(&stamp->get(), component.get());
    }
  };
  collect(mMNAComponents);
  collect(mMNAIntfVariableComps);

  // Resolve compact-vs-dense accumulation once. The solve loop then follows
  // this ordered list directly instead of hashing every RHS stamp each step.
  mRightVectorStampOperations.clear();
  mRightVectorStampOperations.reserve(mRightVectorStamps.size());
  for (const Matrix *stamp : mRightVectorStamps) {
    const auto compact = compactStamps.find(stamp);
    if (compact != compactStamps.end())
      mRightVectorStampOperations.push_back({nullptr, compact->second});
    else
      mRightVectorStampOperations.push_back({stamp, nullptr});
  }
}

template <typename VarType>
void MnaSolverDirect<VarType>::addRightVectorStamps(Matrix &rightVector) {
  using Clock = std::chrono::steady_clock;
  // The coarse solver mode times the complete RHS assembly at its caller.
  // Keep individual compact callbacks uninstrumented in that mode.
  const bool profile = Solver::mLogSolveTimes;
  const auto start = profile ? Clock::now() : Clock::time_point{};
  std::chrono::nanoseconds compactTime{};
  for (const auto &operation : mRightVectorStampOperations) {
    if (operation.compact) {
      const auto compactStart = profile ? Clock::now() : Clock::time_point{};
      operation.compact->mnaAddCompactRightSideVectorStamp(rightVector);
      if (profile)
        compactTime += std::chrono::duration_cast<std::chrono::nanoseconds>(
            Clock::now() - compactStart);
      ++mCompactRightVectorStampCalls;
    } else {
      rightVector += *operation.dense;
      ++mDenseRightVectorStampCalls;
    }
  }
  if (profile) {
    const auto end = Clock::now();
    recordProfileMeasurement("rhs_compact_stamp_callbacks", start,
                             start + compactTime);
    recordProfileMeasurement("rhs_dense_stamp_accumulation",
                             start + compactTime, end);
  }
}

template <typename VarType>
std::shared_ptr<CPS::Task>
MnaSolverDirect<VarType>::createSolveTaskHarm(UInt freqIdx) {
  return std::make_shared<MnaSolverDirect<VarType>::SolveTaskHarm>(*this,
                                                                   freqIdx);
}

template <typename VarType>
std::shared_ptr<CPS::Task>
MnaSolverDirect<VarType>::createStateSpaceExtractionTask() {
  return std::make_shared<
      typename MnaSolverDirect<VarType>::StateSpaceExtractionTask>(*this);
}

template <typename VarType>
std::shared_ptr<CPS::Task> MnaSolverDirect<VarType>::createLogTask() {
  return std::make_shared<MnaSolverDirect<VarType>::LogTask>(*this);
}

template <typename VarType>
void MnaSolverDirect<VarType>::solve(Real time, Int timeStepCount) {
  const bool profile = Solver::mLogSolveTimes || Solver::mSolverPhaseProfiling;
  const bool aggregateRhs =
      Solver::mSolverPhaseProfiling && !Solver::mLogSolveTimes;
  const auto previousSwitchStatus = mCurrentSwitchStatus;
  // Reset source vector
  auto profileStart = profile ? std::chrono::steady_clock::now()
                              : std::chrono::steady_clock::time_point{};
  mRightSideVector.setZero();
  if (profile)
    recordProfileMeasurement("rhs_zero", profileStart,
                             std::chrono::steady_clock::now());

  // Add together the right side vector (computed by the components' pre-step tasks)
  profileStart = aggregateRhs ? std::chrono::steady_clock::now()
                             : std::chrono::steady_clock::time_point{};
  addRightVectorStamps(mRightSideVector);
  if (aggregateRhs)
    recordProfileMeasurement("rhs_assembly", profileStart,
                             std::chrono::steady_clock::now());

  if (!mIsInInitialization) {
    profileStart = profile ? std::chrono::steady_clock::now()
                           : std::chrono::steady_clock::time_point{};
    MnaSolver<VarType>::updateSwitchStatus();
    if (profile)
      recordProfileMeasurement("switch_status_update", profileStart,
                               std::chrono::steady_clock::now());
  }

  if (mSwitchedMatrices.size() > 0) {
    std::chrono::steady_clock::time_point start;
    if (profile)
      start = std::chrono::steady_clock::now();

    mDirectLinearSolvers[mCurrentSwitchStatus][0]->solveInto(
        mRightSideVector, *mLeftSideVectorValue);

    if (profile) {
      auto end = std::chrono::steady_clock::now();
      std::chrono::duration<Real> diff = end - start;
      mSolveTimes.push_back(diff.count());
      recordProfileMeasurement("numeric_solve", start, end);
    }
  }

  // CHECK: Is this really required? Or can operations actually become part of
  // correctorStep and mnaPostStep?
  profileStart = profile ? std::chrono::steady_clock::now()
                         : std::chrono::steady_clock::time_point{};
  for (auto syncGen : mSyncGen)
    syncGen->updateVoltage(*mLeftSideVectorValue);
  if (profile)
    recordProfileMeasurement("iterative_generator_voltage_update", profileStart,
                             std::chrono::steady_clock::now());

  // Reset number of iterations
  mIter = 0;

  // Additional solve steps for iterative models
  if (mSyncGen.size() > 0) {
    UInt numCompsRequireIter;
    do {
      // count synchronous generators that require iteration
      numCompsRequireIter = 0;
      profileStart = profile ? std::chrono::steady_clock::now()
                             : std::chrono::steady_clock::time_point{};
      for (auto syncGen : mSyncGen)
        if (syncGen->requiresIteration())
          numCompsRequireIter++;
      if (profile)
        recordProfileMeasurement("iterative_generator_convergence_check",
                                 profileStart,
                                 std::chrono::steady_clock::now());

      // recompute solve step if at least one component demands iteration
      if (numCompsRequireIter > 0) {
        mIter++;

        // Reset source vector
        mRightSideVector.setZero();

        if (!mIsInInitialization)
          MnaSolver<VarType>::updateSwitchStatus();

        profileStart = profile ? std::chrono::steady_clock::now()
                               : std::chrono::steady_clock::time_point{};
        for (auto syncGen : mSyncGen)
          syncGen->correctorStep();
        if (profile)
          recordProfileMeasurement("iterative_generator_corrector",
                                   profileStart,
                                   std::chrono::steady_clock::now());

        // Add together the right side vector (computed by the components' pre-step tasks)
        profileStart = aggregateRhs ? std::chrono::steady_clock::now()
                                   : std::chrono::steady_clock::time_point{};
        addRightVectorStamps(mRightSideVector);
        if (aggregateRhs)
          recordProfileMeasurement("rhs_assembly", profileStart,
                                   std::chrono::steady_clock::now());

        if (mSwitchedMatrices.size() > 0) {
          const auto start = profile ? std::chrono::steady_clock::now()
                                     : std::chrono::steady_clock::time_point{};
          mDirectLinearSolvers[mCurrentSwitchStatus][0]->solveInto(
              mRightSideVector, *mLeftSideVectorValue);
          if (profile) {
            const auto end = std::chrono::steady_clock::now();
            mSolveTimes.push_back(
                std::chrono::duration<Real>(end - start).count());
            recordProfileMeasurement("numeric_solve", start, end);
          }
        }

        // CHECK: Is this really required? Or can operations actually become part of
        // correctorStep and mnaPostStep?
        profileStart = profile ? std::chrono::steady_clock::now()
                               : std::chrono::steady_clock::time_point{};
        for (auto syncGen : mSyncGen)
          syncGen->updateVoltage(*mLeftSideVectorValue);
        if (profile)
          recordProfileMeasurement("iterative_generator_voltage_update",
                                   profileStart,
                                   std::chrono::steady_clock::now());
      }
    } while (numCompsRequireIter > 0);
  }

  mSystemMatrixChanged = (mCurrentSwitchStatus != previousSwitchStatus);

  // TODO split into separate task? (dependent on x, updating all v attributes)
  profileStart = profile ? std::chrono::steady_clock::now()
                         : std::chrono::steady_clock::time_point{};
  for (UInt nodeIdx = 0; nodeIdx < mNumNetNodes; ++nodeIdx)
    mNodes[nodeIdx]->mnaUpdateVoltage(*mLeftSideVectorValue);
  if (profile)
    recordProfileMeasurement("node_voltage_propagation", profileStart,
                             std::chrono::steady_clock::now());

  // Components' states will be updated by the post-step tasks
}

template <typename VarType>
void MnaSolverDirect<VarType>::solveWithHarmonics(Real time, Int timeStepCount,
                                                  Int freqIdx) {
  mRightSideVectorHarm[freqIdx].setZero();

  // Sum of right side vectors (computed by the components' pre-step tasks)
  for (auto stamp : mRightVectorStamps)
    mRightSideVectorHarm[freqIdx] += stamp->col(freqIdx);

  // Frequency solve tasks execute concurrently. Derive the switch-state key
  // locally instead of mutating the shared mCurrentSwitchStatus bitset from
  // each frequency worker.
  std::bitset<SWITCH_NUM> switchStatus;
  for (UInt i = 0; i < mSwitches.size(); ++i)
    switchStatus.set(i, mSwitches[i]->mnaIsClosed());

  mDirectLinearSolvers[switchStatus][freqIdx]->solveInto(
      mRightSideVectorHarm[freqIdx], **mLeftSideVectorHarm[freqIdx]);
}

template <typename VarType> void MnaSolverDirect<VarType>::logSystemMatrices() {
  // Matrix formatting materializes sparse matrices as dense matrices. Avoid
  // allocating and formatting suppressed diagnostics during initialization.
  if (Solver::mLogLevel > CPS::Logger::Level::info)
    return;
  if (mFrequencyParallel) {
    for (UInt i = 0; i < mSwitchedMatrices[std::bitset<SWITCH_NUM>(0)].size();
         ++i) {
      SPDLOG_LOGGER_INFO(mSLog, "System matrix for frequency: {:d} \n{:s}", i,
                         Logger::matrixToString(
                             mSwitchedMatrices[std::bitset<SWITCH_NUM>(0)][i]));
    }

    for (UInt i = 0; i < mRightSideVectorHarm.size(); ++i) {
      SPDLOG_LOGGER_INFO(mSLog, "Right side vector for frequency: {:d} \n{:s}",
                         i, Logger::matrixToString(mRightSideVectorHarm[i]));
    }

  } else if (mSystemMatrixRecomputationEnabled) {
    SPDLOG_LOGGER_INFO(mSLog, "Summarizing matrices: ");
    SPDLOG_LOGGER_INFO(mSLog, "Base matrix with only static elements: {}",
                       Logger::matrixToString(mBaseSystemMatrix));
    SPDLOG_LOGGER_INFO(mSLog, "Initial system matrix with variable elements {}",
                       Logger::matrixToString(mVariableSystemMatrix));
    SPDLOG_LOGGER_INFO(mSLog, "Right side vector: {}",
                       Logger::matrixToString(mRightSideVector));
  } else {
    if (mSwitches.size() < 1) {
      SPDLOG_LOGGER_INFO(mSLog, "System matrix: \n{}",
                         mSwitchedMatrices[std::bitset<SWITCH_NUM>(0)][0]);
    } else {
      SPDLOG_LOGGER_INFO(mSLog, "Initial switch status: {:s}",
                         mCurrentSwitchStatus.to_string());

      for (auto sys : mSwitchedMatrices) {
        SPDLOG_LOGGER_INFO(mSLog, "Switching System matrix {:s} \n{:s}",
                           sys.first.to_string(),
                           Logger::matrixToString(sys.second[0]));
      }
    }
    SPDLOG_LOGGER_INFO(mSLog, "Right side vector: \n{}", mRightSideVector);
  }
}

template <typename VarType> void MnaSolverDirect<VarType>::logLUTimes() {
  logFactorizationTime();
  logRecomputationTime();
  logSolveTime();
  for (const auto &entry : mProfileMeasurements) {
    const auto &measurement = entry.second;
    const double averageNs =
        measurement.count == 0
            ? 0.
            : static_cast<double>(measurement.total.count()) /
                  measurement.count;
    SPDLOG_LOGGER_INFO(mSLog,
                       "DEEP_TIMING_PHASE,{},{:.9f},{:.3f},{:.3f},{:.3f},{}",
                       entry.first, measurement.total.count() * 1e-9,
                       averageNs * 1e-3, measurement.minimum.count() * 1e-3,
                       measurement.maximum.count() * 1e-3, measurement.count);
  }
  std::ofstream phaseFile(CPS::Logger::logDir() + "/" + this->mName +
                          "_SolverPhaseTimes.csv");
  phaseFile << "phase,total_s,avg_us,min_us,max_us,calls\n";
  for (const auto &entry : mProfileMeasurements) {
    const auto &measurement = entry.second;
    const double count = static_cast<double>(measurement.count);
    phaseFile << entry.first << "," << measurement.total.count() * 1e-9 << ","
              << (count ? measurement.total.count() / count * 1e-3 : 0.) << ","
              << measurement.minimum.count() * 1e-3 << ","
              << measurement.maximum.count() * 1e-3 << "," << measurement.count
              << "\n";
  }
  phaseFile << "rhs_dense_stamp_component_calls,0,0,0,0,"
            << mDenseRightVectorStampCalls << "\n";
  phaseFile << "rhs_compact_stamp_component_calls,0,0,0,0,"
            << mCompactRightVectorStampCalls << "\n";
  // Top-level composites recurse into their own primitive children. Do not
  // also traverse those children here, or their operation times double count.
  for (const auto &component : mMNAComponents)
    component->mnaWriteProfile(phaseFile);
  for (const auto &component : mMNAIntfVariableComps)
    component->mnaWriteProfile(phaseFile);
  std::shared_ptr<DirectLinearSolver> activeSolver =
      mDirectLinearSolverVariableSystemMatrix;
  if (!mSystemMatrixRecomputationEnabled) {
    const auto found = mDirectLinearSolvers.find(mCurrentSwitchStatus);
    if (found != mDirectLinearSolvers.end() && !found->second.empty())
      activeSolver = found->second.front();
  }
  if (auto dense = std::dynamic_pointer_cast<DenseLUAdapter>(activeSolver))
    dense->writeDetailedProfile(phaseFile);
  std::ofstream diagnostics(CPS::Logger::logDir() + "/" + this->mName +
                            "_SolverDiagnostics.csv");
  if (auto dense = std::dynamic_pointer_cast<DenseLUAdapter>(activeSolver))
    dense->writeDetailedDiagnostics(diagnostics);
#ifdef WITH_KLU
  else if (auto klu = std::dynamic_pointer_cast<KLUAdapter>(activeSolver))
    klu->writeDetailedDiagnostics(diagnostics);
#endif
  else
    diagnostics << "metric,value,unit\n";
  diagnostics << "unchanged_factor_reuses," << mUnchangedFactorReuseCount
              << ",steps\n"
              << "matrix_reset_and_stamp_calls," << mMatrixStampCount
              << ",calls\n"
              << "detailed_profiling_enabled," << Solver::mLogSolveTimes
              << ",bool\n"
              << "solver_phase_profiling_enabled,"
              << (Solver::mLogSolveTimes || Solver::mSolverPhaseProfiling)
              << ",bool\n"
              << "rhs_dense_component_calls," << mDenseRightVectorStampCalls
              << ",calls\n"
              << "rhs_compact_component_calls," << mCompactRightVectorStampCalls
              << ",calls\n";
  std::ofstream matrixFile(CPS::Logger::logDir() + "/" + this->mName +
                           "_SolverMatrixStats.csv");
  const SparseMatrix *profiledMatrix = &mVariableSystemMatrix;
  const char *matrixSource = "variable_system_matrix";
  if (!mSystemMatrixRecomputationEnabled) {
    const auto activeMatrix = mSwitchedMatrices.find(mCurrentSwitchStatus);
    if (activeMatrix != mSwitchedMatrices.end() &&
        !activeMatrix->second.empty()) {
      profiledMatrix = &activeMatrix->second.front();
      matrixSource = "active_switched_matrix";
    }
  }
  matrixFile << "source,rows,cols,nonzeros,structural_density,variable_entries,"
                "numerical_nonzeros\n";
  const auto numericalNonzeros =
      std::count_if(profiledMatrix->valuePtr(),
                    profiledMatrix->valuePtr() + profiledMatrix->nonZeros(),
                    [](Real value) { return value != 0.; });
  const double matrixSize = static_cast<double>(profiledMatrix->rows());
  const double matrixEntries = matrixSize * matrixSize;
  matrixFile << matrixSource << "," << profiledMatrix->rows() << ","
             << profiledMatrix->cols() << "," << profiledMatrix->nonZeros()
             << ","
             << (matrixEntries ? profiledMatrix->nonZeros() / matrixEntries
                               : 0.)
             << "," << mListVariableSystemMatrixEntries.size() << ","
             << numericalNonzeros << "\n";
  if (!mFactorizeTimes.empty()) {
    const Real total = std::accumulate(mFactorizeTimes.begin(),
                                       mFactorizeTimes.end(), Real(0.));
    const auto bounds =
        std::minmax_element(mFactorizeTimes.begin(), mFactorizeTimes.end());
    phaseFile << "initial_factorization," << total << ","
              << total * 1e6 / mFactorizeTimes.size() << ","
              << *bounds.first * 1e6 << "," << *bounds.second * 1e6 << ","
              << mFactorizeTimes.size() << "\n";
  }
}

template <typename VarType>
void MnaSolverDirect<VarType>::recordProfileMeasurement(
    const String &phase, std::chrono::steady_clock::time_point start,
    std::chrono::steady_clock::time_point end) {
  const auto elapsed =
      std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
  auto &measurement = mProfileMeasurements[phase];
  if (measurement.count == 0) {
    measurement.minimum = elapsed;
    measurement.maximum = elapsed;
  } else {
    measurement.minimum = std::min(measurement.minimum, elapsed);
    measurement.maximum = std::max(measurement.maximum, elapsed);
  }
  measurement.total += elapsed;
  ++measurement.count;
}

template <typename VarType> void MnaSolverDirect<VarType>::logSolveTime() {
  Real solveSum = 0.0;
  Real solveMax = 0.0;
  for (auto meas : mSolveTimes) {
    solveSum += meas;
    if (meas > solveMax)
      solveMax = meas;
  }
  SPDLOG_LOGGER_INFO(mSLog, "Cumulative solve times: {:.12f}", solveSum);
  SPDLOG_LOGGER_INFO(mSLog, "Average solve time: {:.12f}",
                     mSolveTimes.empty()
                         ? 0.
                         : solveSum / static_cast<double>(mSolveTimes.size()));
  SPDLOG_LOGGER_INFO(mSLog, "Maximum solve time: {:.12f}", solveMax);
  SPDLOG_LOGGER_INFO(mSLog, "Number of solves: {:d}", mSolveTimes.size());
}

template <typename VarType>
void MnaSolverDirect<VarType>::logFactorizationTime() {
  for (auto meas : mFactorizeTimes) {
    SPDLOG_LOGGER_INFO(mSLog, "LU factorization time: {:.12f}", meas);
  }
}

template <typename VarType>
void MnaSolverDirect<VarType>::logRecomputationTime() {
  Real recompSum = 0.0;
  Real recompMax = 0.0;
  for (auto meas : mRecomputationTimes) {
    recompSum += meas;
    if (meas > recompMax)
      recompMax = meas;
  }
  // Sometimes, refactorization is not used
  if (mRecomputationTimes.size() != 0) {
    SPDLOG_LOGGER_INFO(mSLog, "Cumulative refactorization times: {:.12f}",
                       recompSum);
    SPDLOG_LOGGER_INFO(mSLog, "Average refactorization time: {:.12f}",
                       recompSum / ((double)mRecomputationTimes.size()));
    SPDLOG_LOGGER_INFO(mSLog, "Maximum refactorization time: {:.12f}",
                       recompMax);
    SPDLOG_LOGGER_INFO(mSLog, "Number of refactorizations: {:d}",
                       mRecomputationTimes.size());
  }
}

template <typename VarType>
std::shared_ptr<DirectLinearSolver>
MnaSolverDirect<VarType>::createDirectSolverImplementation(
    CPS::Logger::Log mSLog) {
  switch (this->mImplementationInUse) {
  case DirectLinearSolverImpl::DenseLU:
    return std::make_shared<DenseLUAdapter>(mSLog);
#ifdef WITH_SPARSE
  case DirectLinearSolverImpl::SparseLU:
    return std::make_shared<SparseLUAdapter>(mSLog);
#endif // WITH_SPARSE
#ifdef WITH_KLU
  case DirectLinearSolverImpl::KLU:
    return std::make_shared<KLUAdapter>(mSLog);
#endif
#ifdef WITH_CUDA
  case DirectLinearSolverImpl::CUDADense:
    return std::make_shared<GpuDenseAdapter>(mSLog);
#ifdef WITH_CUDA_SPARSE
  case DirectLinearSolverImpl::CUDASparse:
    return std::make_shared<GpuSparseAdapter>(mSLog);
#endif
#ifdef WITH_MAGMA
  case DirectLinearSolverImpl::CUDAMagma:
    return std::make_shared<GpuMagmaAdapter>(mSLog);
#endif
#endif
#ifdef WITH_METAL
  case DirectLinearSolverImpl::MetalDense:
    return std::make_shared<MetalDenseAdapter>(mSLog);
  case DirectLinearSolverImpl::MetalSparseBiCGSTAB:
    return std::make_shared<MetalSparseBiCGSTABAdapter>(mSLog);
#endif
  default:
    throw CPS::SystemError("unsupported linear solver implementation.");
  }
}

template <typename VarType>
void MnaSolverDirect<VarType>::setDirectLinearSolverImplementation(
    DirectLinearSolverImpl implementation) {
  this->mImplementationInUse = implementation;
}

template <typename VarType>
void MnaSolverDirect<VarType>::setDirectLinearSolverConfiguration(
    DirectLinearSolverConfiguration &configuration) {
  this->mConfigurationInUse = configuration;
}

} // namespace DPsim

template class DPsim::MnaSolverDirect<Real>;
template class DPsim::MnaSolverDirect<Complex>;
