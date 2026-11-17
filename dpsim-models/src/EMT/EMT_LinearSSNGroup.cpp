// SPDX-License-Identifier: MPL-2.0
#include <dpsim-models/EMT/EMT_LinearSSNGroup.h>

#include <chrono>
#include <cmath>

#include <dpsim-models/EMT/EMT_Ph3_SSN_ImpedanceLoad.h>
#include <dpsim-models/EMT/EMT_Ph3_SSN_PiLine.h>
#include <dpsim-models/EMT/EMT_Ph3_SSN_Transformer.h>
#include <dpsim-models/MNAStampUtils.h>

#include <stdexcept>

using namespace CPS;

namespace {
using PhaseVector = Eigen::Matrix<Real, 3, 1>;

Matrix modelHistory(const EMT::LinearSSNModel &model,
                    Eigen::Ref<const Matrix> state) {
  if (model.kind == EMT::LinearSSNModel::Kind::PiLine) {
    const PhaseVector current = state.block<3, 1>(0, 0);
    const PhaseVector inductorVoltage = state.block<3, 1>(3, 0);
    const PhaseVector capVoltage0 = state.block<3, 1>(6, 0);
    const PhaseVector capCurrent0 = state.block<3, 1>(9, 0);
    const PhaseVector capVoltage1 = state.block<3, 1>(12, 0);
    const PhaseVector capCurrent1 = state.block<3, 1>(15, 0);
    const PhaseVector series =
        model.seriesSystemInverse *
        (model.inductorAdmittance * inductorVoltage + current);
    Matrix result(6, 1);
    result.block<3, 1>(0, 0) =
        series - (-capCurrent0 - model.capAdmittance * capVoltage0);
    result.block<3, 1>(3, 0) =
        -series - (-capCurrent1 - model.capAdmittance * capVoltage1);
    return result;
  }
  if (model.kind == EMT::LinearSSNModel::Kind::ImpedanceLoad) {
    Matrix result = Matrix::Zero(3, 1);
    for (UInt phase = 0; phase < 3; ++phase) {
      Real loadHistory =
          model.hasLoadInductance
              ? model.inductorAdmittanceScalar * state.coeff(phase, 0) +
                    state.coeff(3 + phase, 0)
              : 0.;
      Real shuntHistory =
          model.hasLoadCapacitance
              ? -state.coeff(9 + phase, 0) -
                    model.capacitorAdmittanceScalar * state.coeff(6 + phase, 0)
              : 0.;
      result(phase) = -(loadHistory + shuntHistory);
    }
    return result;
  }

  const PhaseVector h0 = model.seriesSystemInverse *
                         (model.inductorAdmittance * state.block<3, 1>(0, 0) +
                          state.block<3, 1>(3, 0));
  const PhaseVector h1 = model.seriesSystemInverse *
                         (model.inductorAdmittance * state.block<3, 1>(6, 0) +
                          state.block<3, 1>(9, 0));
  const PhaseVector hMag =
      model.magnetizingInductorAdmittance * state.block<3, 1>(12, 0) +
      state.block<3, 1>(15, 0);
  const PhaseVector midpointSource = model.midpointInverse * (-h0 + h1 + hMag);
  Matrix result(6, 1);
  result.block<3, 1>(0, 0) = h0 + model.seriesAdmittance * midpointSource;
  result.block<3, 1>(3, 0) =
      model.ratioFromReference * (model.seriesAdmittance * midpointSource - h1);
  if (model.referenceTerminal == 1)
    result.block<3, 1>(0, 0).swap(result.block<3, 1>(3, 0));
  return result;
}

Matrix modelNext(const EMT::LinearSSNModel &model,
                 Eigen::Ref<const Matrix> state,
                 Eigen::Ref<const Matrix> input) {
  if (model.kind == EMT::LinearSSNModel::Kind::PiLine) {
    const PhaseVector v0 = input.block<3, 1>(0, 0);
    const PhaseVector v1 = input.block<3, 1>(3, 0);
    const PhaseVector interfaceVoltage = v1 - v0;
    const PhaseVector series =
        model.seriesSystemInverse *
        (model.inductorAdmittance * state.block<3, 1>(3, 0) +
         state.block<3, 1>(0, 0));
    const PhaseVector current =
        model.seriesAdmittance * interfaceVoltage + series;
    const PhaseVector cap0 = model.capAdmittance * v0 +
                             (-state.block<3, 1>(9, 0) -
                              model.capAdmittance * state.block<3, 1>(6, 0));
    const PhaseVector cap1 = model.capAdmittance * v1 +
                             (-state.block<3, 1>(15, 0) -
                              model.capAdmittance * state.block<3, 1>(12, 0));
    Matrix result(18, 1);
    result.block<3, 1>(0, 0) = current;
    result.block<3, 1>(3, 0) =
        interfaceVoltage - model.seriesResistance * current;
    result.block<3, 1>(6, 0) = v0;
    result.block<3, 1>(9, 0) = cap0;
    result.block<3, 1>(12, 0) = v1;
    result.block<3, 1>(15, 0) = cap1;
    return result;
  }
  if (model.kind == EMT::LinearSSNModel::Kind::ImpedanceLoad) {
    Matrix result = state;
    for (UInt phase = 0; phase < 3; ++phase) {
      const Real voltage = input.coeff(phase, 0);
      if (model.hasLoadInductance) {
        result(phase) = voltage;
        result(3 + phase) =
            model.inductorAdmittanceScalar * voltage +
            (model.inductorAdmittanceScalar * state.coeff(phase, 0) +
             state.coeff(3 + phase, 0));
      }
      if (model.hasLoadCapacitance) {
        result(6 + phase) = voltage;
        result(9 + phase) =
            model.capacitorAdmittanceScalar * voltage -
            state.coeff(9 + phase, 0) -
            model.capacitorAdmittanceScalar * state.coeff(6 + phase, 0);
      }
    }
    return result;
  }

  PhaseVector v0 = input.block<3, 1>(0, 0);
  PhaseVector v1 = input.block<3, 1>(3, 0);
  const UInt ref = model.referenceTerminal;
  const PhaseVector referenceVoltage = ref == 0 ? v0 : v1;
  const PhaseVector otherVoltage = ref == 0 ? v1 : v0;
  const PhaseVector referredOther = model.ratioFromReference * otherVoltage;
  const PhaseVector h0 = model.seriesSystemInverse *
                         (model.inductorAdmittance * state.block<3, 1>(0, 0) +
                          state.block<3, 1>(3, 0));
  const PhaseVector h1 = model.seriesSystemInverse *
                         (model.inductorAdmittance * state.block<3, 1>(6, 0) +
                          state.block<3, 1>(9, 0));
  const PhaseVector hMag =
      model.magnetizingInductorAdmittance * state.block<3, 1>(12, 0) +
      state.block<3, 1>(15, 0);
  const PhaseVector midpoint =
      model.midpointInverse *
      (model.seriesAdmittance * referenceVoltage +
       model.seriesAdmittance * referredOther - h0 + h1 + hMag);
  const PhaseVector total0 = midpoint - referenceVoltage;
  const PhaseVector current0 = model.seriesAdmittance * total0 + h0;
  const PhaseVector total1 = referredOther - midpoint;
  const PhaseVector current1 = model.seriesAdmittance * total1 + h1;
  Matrix result(18, 1);
  result.block<3, 1>(0, 0) = total0 - model.seriesResistanceHalf * current0;
  result.block<3, 1>(3, 0) = current0;
  result.block<3, 1>(6, 0) = total1 - model.seriesResistanceHalf * current1;
  result.block<3, 1>(9, 0) = current1;
  result.block<3, 1>(12, 0) = -midpoint;
  result.block<3, 1>(15, 0) =
      model.magnetizingInductorAdmittance * (-midpoint) + hMag;
  return result;
}

void addBlock(Matrix &target, const std::vector<UInt> &rows,
              const std::vector<UInt> &cols, const Matrix &block) {
  for (UInt i = 0; i < rows.size(); ++i)
    for (UInt j = 0; j < cols.size(); ++j)
      target(rows[i], cols[j]) += block(i, j);
}
} // namespace

EMT::LinearSSNGroup::LinearSSNGroup(String name, ComponentList components,
                                    SimNode::List allNodes,
                                    SimNode::List retainedNodes,
                                    Logger::Level level)
    : MNASimPowerComp<Real>(name, true, true, level),
      mComponents(std::move(components)), mAllNodes(std::move(allNodes)) {
  mPhaseType = PhaseType::ABC;
  const UInt retainedCount = static_cast<UInt>(retainedNodes.size());
  setTerminalNumber(retainedCount);
  connect(std::move(retainedNodes));
  for (UInt i = 0; i < retainedCount; ++i)
    mRetainedNames.push_back(node(i)->name());
  for (UInt i = 0; i < mRetainedNames.size(); ++i)
    mRetainedNodeIndex[mRetainedNames[i]] = i;
  for (UInt i = 0; i < mAllNodes.size(); ++i)
    mLocalNodeIndex[mAllNodes[i]->name()] = i;
  mState = mAttributes->create<Matrix>("group_state");
  mState->get().resize(0, 1);
}

Bool EMT::LinearSSNGroup::publishesNodeVoltage(const SimNode::Ptr &candidate) {
  const auto retained = mRetainedNodeIndex.find(candidate->name());
  // The surrounding topology may use separate terminal objects sharing the
  // original bus names. Those original bus attributes also need publication.
  return retained == mRetainedNodeIndex.end() ||
         node(retained->second) != candidate;
}

void EMT::LinearSSNGroup::initializeFromNodesAndTerminals(Real frequency) {
  // Eliminated buses are absent from the surrounding topology's node list.
  // Initialize their published voltage just as the retained buses are initialized.
  for (const auto &node : mAllNodes)
    if (publishesNodeVoltage(node))
      node->initialize();
  for (const auto &component : mComponents)
    component->initializeFromNodesAndTerminals(frequency);
}

EMT::LinearSSNGroup::PrecomputedStateMaps
EMT::LinearSSNGroup::precomputedStateMaps() const {
  if (!mReady || !mUseFullyReducedInternalMaps)
    throw std::logic_error(
        "LinearSSNGroup accelerator snapshot requires initialized, fully "
        "reduced state maps.");

  PrecomputedStateMaps snapshot;
  snapshot.historyMap = mReducedHistoryMap;
  snapshot.stateMap = mReducedStateMap;
  snapshot.inputMap = mReducedInputMap;
  snapshot.terminalAdmittance = mDd;
  snapshot.initialState = mState->get();
  snapshot.terminalMatrixIndices = mRetainedMatrixIndices;
  return snapshot;
}

void EMT::LinearSSNGroup::buildRealization(Real omega, Real timeStep,
                                           Attribute<Matrix>::Ptr leftVector) {
  const auto preparationStart = std::chrono::steady_clock::now();
  mModels.clear();
  mModelOutputs.clear();
  mNodeVoltages.clear();
  for (const auto &node : mAllNodes) {
    // Actual retained terminal attributes are written by node propagation;
    // hidden buses and distinct original-bus aliases belong to this post-step.
    if (publishesNodeVoltage(node)) {
      mNodeVoltages.push_back(&node->mVoltage->get());
    } else
      mNodeVoltages.push_back(nullptr);
  }
  for (const auto &component : mComponents) {
    component->mnaInitialize(omega, timeStep, leftVector);
    LinearSSNModel model;
    if (auto line = std::dynamic_pointer_cast<Ph3::SSN::PiLine>(component))
      line->exportFusedModel(model);
    else if (auto load =
                 std::dynamic_pointer_cast<Ph3::SSN::ImpedanceLoad>(component))
      load->exportFusedModel(model);
    else if (auto transformer =
                 std::dynamic_pointer_cast<Ph3::SSN::Transformer>(component))
      transformer->exportFusedModel(model);
    else
      throw std::invalid_argument("Unsupported component in LinearSSNGroup.");
    ModelOutputs outputs;
    const auto cacheOutput = [&](const String &name) {
      const auto attr = component->attributeTyped<Matrix>(name);
      return &attr->get();
    };
    outputs.interfaceVoltage = cacheOutput("v_intf");
    outputs.interfaceCurrent = cacheOutput("i_intf");
    const auto updateCount = component->attributeTyped<UInt>("output_updates");
    outputs.updateCount = &updateCount->get();
    for (UInt terminal = 0; terminal < model.terminalCount; ++terminal)
      outputs.terminalCurrent[terminal] =
          cacheOutput("i_terminal_" + std::to_string(terminal));
    mModelOutputs.push_back(outputs);
    mModels.push_back(std::move(model));
  }

  const UInt localNodeCount = static_cast<UInt>(mAllNodes.size()) * 3;
  const UInt stateCount = [&] {
    UInt n = 0;
    for (const auto &model : mModels)
      n += model.state.rows();
    return n;
  }();
  Matrix G = Matrix::Zero(localNodeCount, localNodeCount);
  Matrix H = Matrix::Zero(localNodeCount, stateCount);
  Matrix F = Matrix::Zero(stateCount, stateCount);
  Matrix B = Matrix::Zero(stateCount, localNodeCount);
  mStateOffsets.clear();
  mStateOffsets.reserve(mModels.size());
  mModelInputIndices.clear();
  mModelInputIndices.reserve(mModels.size());
  Matrix initialState = Matrix::Zero(stateCount, 1);
  UInt stateOffset = 0;
  for (const auto &model : mModels) {
    mStateOffsets.push_back(stateOffset);
    std::vector<UInt> inputIndices(model.terminalCount * 3,
                                   std::numeric_limits<UInt>::max());
    for (UInt p = 0; p < inputIndices.size(); ++p) {
      const UInt terminal = p / 3;
      if (!model.nodeNames[p].empty() && !model.grounded[terminal])
        inputIndices[p] = mLocalNodeIndex.at(model.nodeNames[p]) * 3 + p % 3;
    }
    mModelInputIndices.push_back(std::move(inputIndices));
    initialState.block(stateOffset, 0, model.state.rows(), 1) = model.state;
    const UInt inputCount = model.terminalCount * 3;
    Matrix zeroState = Matrix::Zero(model.state.rows(), 1);
    Matrix zeroInput = Matrix::Zero(inputCount, 1);
    Matrix localF = Matrix::Zero(model.state.rows(), model.state.rows());
    Matrix localB = Matrix::Zero(model.state.rows(), inputCount);
    Matrix localH = Matrix::Zero(inputCount, model.state.rows());
    for (UInt i = 0; i < model.state.rows(); ++i) {
      Matrix basis = Matrix::Zero(model.state.rows(), 1);
      basis(i) = 1.;
      localF.col(i) = modelNext(model, basis, zeroInput);
      localH.col(i) = modelHistory(model, basis);
    }
    for (UInt i = 0; i < inputCount; ++i) {
      Matrix basis = Matrix::Zero(inputCount, 1);
      basis(i) = 1.;
      localB.col(i) = modelNext(model, zeroState, basis);
    }
    F.block(stateOffset, stateOffset, model.state.rows(), model.state.rows()) =
        localF;
    for (UInt p = 0; p < inputCount; ++p) {
      const String &name = model.nodeNames[p];
      const UInt terminal = p / 3;
      if (name.empty() || model.grounded[terminal])
        continue;
      const UInt base = mLocalNodeIndex.at(name) * 3 + (p % 3);
      B.block(stateOffset, base, model.state.rows(), 1) += localB.col(p);
    }
    for (UInt p = 0; p < inputCount; ++p) {
      const String &rowName = model.nodeNames[p];
      const UInt terminal = p / 3;
      if (rowName.empty() || model.grounded[terminal])
        continue;
      const UInt row = mLocalNodeIndex.at(rowName) * 3 + (p % 3);
      for (UInt q = 0; q < model.state.rows(); ++q)
        H(row, stateOffset + q) += localH(p, q);
    }
    for (UInt p = 0; p < inputCount; ++p) {
      const UInt terminal = p / 3;
      if (model.nodeNames[p].empty() || model.grounded[terminal])
        continue;
      const UInt row = mLocalNodeIndex.at(model.nodeNames[p]) * 3 + (p % 3);
      for (UInt q = 0; q < inputCount; ++q) {
        const UInt qTerminal = q / 3;
        if (model.nodeNames[q].empty() || model.grounded[qTerminal])
          continue;
        const UInt col = mLocalNodeIndex.at(model.nodeNames[q]) * 3 + (q % 3);
        G(row, col) += model.portAdmittance(p, q);
      }
    }
    stateOffset += model.state.rows();
  }
  mState->get() = initialState;

  mTerminalIndices.clear();
  for (const auto &name : mRetainedNames) {
    const UInt localBase = mLocalNodeIndex.at(name) * 3;
    for (UInt phase = 0; phase < 3; ++phase)
      mTerminalIndices.push_back(localBase + phase);
  }
  std::vector<Bool> isTerminal(localNodeCount, false);
  for (UInt index : mTerminalIndices)
    isTerminal[index] = true;
  mInternalIndices.clear();
  for (UInt i = 0; i < localNodeCount; ++i)
    if (!isTerminal[i])
      mInternalIndices.push_back(i);

  Matrix Gii = Matrix::Zero(mInternalIndices.size(), mInternalIndices.size());
  Matrix Gin = Matrix::Zero(mInternalIndices.size(), mTerminalIndices.size());
  Matrix Gni = Matrix::Zero(mTerminalIndices.size(), mInternalIndices.size());
  Matrix Gnn = Matrix::Zero(mTerminalIndices.size(), mTerminalIndices.size());
  Matrix Hi = Matrix::Zero(mInternalIndices.size(), stateCount);
  Matrix Hn = Matrix::Zero(mTerminalIndices.size(), stateCount);
  Matrix Bi = Matrix::Zero(stateCount, mInternalIndices.size());
  Matrix Bn = Matrix::Zero(stateCount, mTerminalIndices.size());
  for (UInt i = 0; i < mInternalIndices.size(); ++i) {
    Hi.row(i) = H.row(mInternalIndices[i]);
    for (UInt j = 0; j < mInternalIndices.size(); ++j)
      Gii(i, j) = G(mInternalIndices[i], mInternalIndices[j]);
    for (UInt j = 0; j < mTerminalIndices.size(); ++j)
      Gin(i, j) = G(mInternalIndices[i], mTerminalIndices[j]);
    Bi.col(i) = B.col(mInternalIndices[i]);
  }
  for (UInt i = 0; i < mTerminalIndices.size(); ++i) {
    Hn.row(i) = H.row(mTerminalIndices[i]);
    for (UInt j = 0; j < mInternalIndices.size(); ++j)
      Gni(i, j) = G(mTerminalIndices[i], mInternalIndices[j]);
    for (UInt j = 0; j < mTerminalIndices.size(); ++j)
      Gnn(i, j) = G(mTerminalIndices[i], mTerminalIndices[j]);
    Bn.col(i) = B.col(mTerminalIndices[i]);
  }
  if (!mInternalIndices.empty()) {
    Eigen::PartialPivLU<Matrix> factor(Gii);
    mInternalFactor = std::move(factor);
    // Q is needed for the terminal Schur complement in either runtime mode.
    mInternalInputMap = mInternalFactor.solve(Gin);
    mInternalInputMapSparse = mInternalInputMap.sparseView(0., 0.);
    if (mUsePrecomputedInternalMaps) {
      // Form both reduced response maps with the already-computed LU factors.
      // Do not form an explicit inverse: these are multiple-right-hand-side
      // solves, and P*x is reused by both the history and state-update paths.
      const Matrix denseHi = Hi;
      const Matrix denseStateMap = mInternalFactor.solve(denseHi);
      mInternalStateMap = denseStateMap.sparseView(0., 0.);
      mInternalStateMap.makeCompressed();
    } else {
      const Matrix identity =
          Matrix::Identity(mInternalIndices.size(), mInternalIndices.size());
      mInternalHistoryInverse = mInternalFactor.solve(identity);
      const Real inverseResidual =
          (Gii * mInternalHistoryInverse - identity).cwiseAbs().maxCoeff();
      // Use the cached response only when its factorization residual satisfies
      // the same absolute scale as the solver's refinement acceptance check.
      // Ill-conditioned realizations safely retain the original LU solve path.
      mUseCachedInternalHistoryInverse =
          std::isfinite(inverseResidual) && inverseResidual <= 1e-12;
    }
  }
  // Constant maps retain every exact nonzero, including roundoff-sized fill.
  // No magnitude threshold or timestep pruning is used in this realization.
  mF = F.sparseView(0., 0.);
  mBi = Bi.sparseView(0., 0.);
  mBn = Bn.sparseView(0., 0.);
  mHi = Hi.sparseView(0., 0.);
  mHn = Hn.sparseView(0., 0.);
  mGin = Gin.sparseView(0., 0.);
  mGni = Gni.sparseView(0., 0.);
  if (mInternalIndices.empty())
    mDd = Gnn;
  else
    mDd = Gnn - Gni * mInternalInputMap;
  if (mUseFullyReducedInternalMaps) {
    if (mInternalIndices.empty()) {
      mReducedHistoryMap = mHn;
      mReducedStateMap = mF;
      mReducedInputMap = mBn;
    } else {
      // Complete the Kron/state elimination during initialization. All
      // products remain sparse; no reduced dense state-space matrix is formed.
      const SparseMatrixRow eliminatedHistory = mGni * mInternalStateMap;
      const SparseMatrixRow eliminatedState = mBi * mInternalStateMap;
      const SparseMatrixRow eliminatedInput = mBi * mInternalInputMapSparse;
      mReducedHistoryMap = mHn - eliminatedHistory;
      mReducedStateMap = mF + eliminatedState;
      mReducedInputMap = mBn - eliminatedInput;
    }
    mReducedHistoryMap.makeCompressed();
    mReducedStateMap.makeCompressed();
    mReducedInputMap.makeCompressed();
  }
  mFullHistoryScratch = Matrix::Zero(localNodeCount, 1);
  mInternalHistoryScratch = Matrix::Zero(mInternalIndices.size(), 1);
  mInternalVoltageScratch = Matrix::Zero(mInternalIndices.size(), 1);
  mLocalVoltageScratch = Matrix::Zero(localNodeCount, 1);
  mTerminalInputScratch = Matrix::Zero(mTerminalIndices.size(), 1);
  mTerminalHistoryScratch = Matrix::Zero(mTerminalIndices.size(), 1);
  mStateUpdateScratch = Matrix::Zero(stateCount, 1);
  mReady = true;
  mRealizationPreparationTime =
      std::chrono::duration_cast<std::chrono::nanoseconds>(
          std::chrono::steady_clock::now() - preparationStart);
}

void EMT::LinearSSNGroup::mnaCompInitialize(Real omega, Real timeStep,
                                            Attribute<Matrix>::Ptr leftVector) {
  updateMatrixNodeIndices();
  buildRealization(omega, timeStep, leftVector);
}

void EMT::LinearSSNGroup::assembleFullHistory(Matrix &fullHistory) const {
  fullHistory.setZero();
  const Matrix &state = mState->get();
  for (UInt modelIndex = 0; modelIndex < mModels.size(); ++modelIndex) {
    const auto &model = mModels[modelIndex];
    const Matrix localHistory =
        modelHistory(model, state.block(mStateOffsets[modelIndex], 0,
                                        model.state.rows(), 1));
    const auto &indices = mModelInputIndices[modelIndex];
    for (UInt p = 0; p < indices.size(); ++p)
      if (indices[p] != std::numeric_limits<UInt>::max())
        fullHistory(indices[p]) += localHistory(p);
  }
}

void EMT::LinearSSNGroup::currentHistory(Matrix &result) const {
  using Clock = std::chrono::steady_clock;
  ++mHistoryCalls;
  auto start = mDetailedProfiling ? Clock::now() : Clock::time_point{};
  const auto record = [&](std::chrono::nanoseconds &duration) {
    if (mDetailedProfiling) {
      const auto now = Clock::now();
      duration += std::chrono::duration_cast<std::chrono::nanoseconds>(now - start);
      start = now;
    }
  };
  if (mInternalIndices.empty()) {
    result.noalias() = mHn * mState->get();
    record(mHistoryMultiplyTime);
    return;
  }
  if (mUseFullyReducedInternalMaps) {
    result.noalias() = mReducedHistoryMap * mState->get();
    record(mHistoryMultiplyTime);
    return;
  }

  if (mUsePrecomputedInternalMaps) {
    // P*x is shared by terminal history evaluation and subsequent recovery.
    mInternalVoltageScratch.noalias() = mInternalStateMap * mState->get();
    mInternalStateResponseReady = true;
    record(mHistoryMultiplyTime);
  } else {
    mInternalHistoryScratch.noalias() = mHi * mState->get();
    record(mHistoryMultiplyTime);
    if (mUseCachedInternalHistoryInverse)
      mInternalVoltageScratch.noalias() =
          mInternalHistoryInverse * mInternalHistoryScratch;
    else {
      mInternalVoltageScratch = mInternalFactor.solve(mInternalHistoryScratch);
      ++mInternalSolveCalls;
      ++mHistoryInternalSolveCalls;
    }
    record(mInternalSolveTime);
  }
  result.noalias() = mHn * mState->get();
  result.noalias() -= mGni * mInternalVoltageScratch;
  record(mReducedHistoryTime);
}

void EMT::LinearSSNGroup::terminalVoltage(const Matrix &leftVector,
                                          Matrix &voltage) {
  voltage.setZero();
  for (UInt terminal = 0; terminal < mRetainedNames.size(); ++terminal)
    for (UInt phase = 0; phase < 3; ++phase)
      voltage(terminal * 3 + phase) = Math::realFromVectorElement(
          leftVector, mRetainedMatrixIndices[terminal * 3 + phase]);
}

void EMT::LinearSSNGroup::mnaCompApplySystemMatrixStamp(
    SparseMatrixRow &systemMatrix) {
  updateMatrixNodeIndices();
  if (!mReady)
    throw std::logic_error("LinearSSNGroup used before initialization.");
  mRetainedMatrixIndices.resize(mRetainedNames.size() * 3);
  for (UInt i = 0; i < mRetainedMatrixIndices.size(); ++i)
    mRetainedMatrixIndices[i] = matrixNodeIndex(i / 3, i % 3);
  for (UInt i = 0; i < mRetainedNames.size() * 3; ++i)
    for (UInt j = 0; j < mRetainedNames.size() * 3; ++j)
      if (mDd(i, j) != 0.)
        Math::addToMatrixElement(systemMatrix, mRetainedMatrixIndices[i],
                                 mRetainedMatrixIndices[j], mDd(i, j));
}

void EMT::LinearSSNGroup::mnaCompApplyRightSideVectorStamp(
    Matrix &rightVector) {
  (void)rightVector;
}

void EMT::LinearSSNGroup::mnaAddCompactRightSideVectorStamp(
    Matrix &rightVector) {
  using Clock = std::chrono::steady_clock;
  // History and scatter clocks are exclusive. Both execute inside the
  // enclosing solver RHS-assembly phase and must not be added to it again.
  currentHistory(mTerminalHistoryScratch);
  const auto start = mDetailedProfiling ? Clock::now() : Clock::time_point{};
  for (UInt i = 0; i < mTerminalHistoryScratch.rows(); ++i)
    Math::addToVectorElement(rightVector, mRetainedMatrixIndices[i],
                             mTerminalHistoryScratch(i));
  if (mDetailedProfiling)
    mRhsScatterTime += std::chrono::duration_cast<std::chrono::nanoseconds>(
        Clock::now() - start);
  ++mRhsScatterCalls;
}

void EMT::LinearSSNGroup::mnaWriteProfile(std::ostream &profile) const {
  const auto write = [&profile](const char *name, std::chrono::nanoseconds time,
                                std::size_t calls) {
    const double total = time.count() * 1e-9;
    const double averageUs = calls ? time.count() * 1e-3 / calls : 0.;
    profile << name << "," << total << "," << averageUs << ",0,0," << calls
            << "\n";
  };
  write("group_history_state_multiply", mHistoryMultiplyTime, mHistoryCalls);
  write("group_initialization_constant_maps", mRealizationPreparationTime, 1);
  write("group_history_internal_solve", mInternalSolveTime,
        mHistoryInternalSolveCalls);
  write("group_history_reduced_output", mReducedHistoryTime, mHistoryCalls);
  write("group_rhs_scatter", mRhsScatterTime, mRhsScatterCalls);
  write("group_state_update", mStateUpdateTime, mStateUpdateCalls);
  write("group_voltage_recovery", mVoltageRecoveryTime, mVoltageRecoveryCalls);
  write("group_output_publication", mOutputPublicationTime,
        mOutputPublicationCalls);
  profile << "group_internal_timestep_lu_solve_count,0,0,0,0,"
          << mInternalSolveCalls << "\n";
  profile << "group_output_component_count,0,0,0,0," << mModels.size() << "\n";
  profile << "group_output_bus_count,0,0,0,0," << mAllNodes.size() << "\n";
  profile << "group_hidden_bus_count,0,0,0,0," << mInternalIndices.size() / 3
          << "\n";
  profile << "group_original_coordinates_retained,0,0,0,0,1\n";
  profile << "group_initialization_internal_dense_lu_factorizations,0,0,0,0,"
          << (mInternalIndices.empty() ? 0 : 1) << "\n";
  profile << "group_base_state_map_nonzeros,0,0,0,0," << mF.nonZeros() << "\n";
  profile << "group_internal_state_input_nonzeros,0,0,0,0," << mBi.nonZeros()
          << "\n";
  profile << "group_terminal_state_input_nonzeros,0,0,0,0," << mBn.nonZeros()
          << "\n";
  profile << "group_terminal_history_map_nonzeros,0,0,0,0," << mHn.nonZeros()
          << "\n";
  profile << "group_precomputed_internal_maps_selected,0,0,0,0,"
          << (mUsePrecomputedInternalMaps ? 1 : 0) << "\n";
  profile << "group_fully_reduced_maps_selected,0,0,0,0,"
          << (mUseFullyReducedInternalMaps ? 1 : 0) << "\n";
  profile << "group_internal_state_map_nonzeros,0,0,0,0,"
          << mInternalStateMap.nonZeros() << "\n";
  profile << "group_internal_input_map_nonzeros,0,0,0,0,"
          << mInternalInputMapSparse.nonZeros() << "\n";
  profile << "group_reduced_history_map_nonzeros,0,0,0,0,"
          << mReducedHistoryMap.nonZeros() << "\n";
  profile << "group_reduced_state_map_nonzeros,0,0,0,0,"
          << mReducedStateMap.nonZeros() << "\n";
  profile << "group_reduced_input_map_nonzeros,0,0,0,0,"
          << mReducedInputMap.nonZeros() << "\n";
  profile << "group_state_count,0,0,0,0," << mState->get().rows() << "\n";
  profile << "group_reduced_state_map_rows,0,0,0,0," << mReducedStateMap.rows()
          << "\n";
  profile << "group_reduced_state_map_columns,0,0,0,0,"
          << mReducedStateMap.cols() << "\n";
  if (mReducedStateMap.rows() > 0) {
    Eigen::Index minNonzerosPerRow = mReducedStateMap.cols();
    Eigen::Index maxNonzerosPerRow = 0;
    for (Eigen::Index row = 0; row < mReducedStateMap.rows(); ++row) {
      const Eigen::Index rowNonzeros =
          mReducedStateMap.outerIndexPtr()[row + 1] -
          mReducedStateMap.outerIndexPtr()[row];
      minNonzerosPerRow = std::min(minNonzerosPerRow, rowNonzeros);
      maxNonzerosPerRow = std::max(maxNonzerosPerRow, rowNonzeros);
    }
    profile << "group_reduced_state_map_min_row_nonzeros,0,0,0,0,"
            << minNonzerosPerRow << "\n";
    profile << "group_reduced_state_map_max_row_nonzeros,0,0,0,0,"
            << maxNonzerosPerRow << "\n";
  }
  profile << "group_history_cached_response_selected,0,0,0,0,"
          << (mUseCachedInternalHistoryInverse ? 1 : 0) << "\n";
}

void EMT::LinearSSNGroup::mnaCompPreStep(Real, Int) {}

void EMT::LinearSSNGroup::mnaCompPostStep(Real, Int,
                                          Attribute<Matrix>::Ptr &leftVector) {
  using Clock = std::chrono::steady_clock;
  auto start = mDetailedProfiling ? Clock::now() : Clock::time_point{};
  terminalVoltage(**leftVector, mTerminalInputScratch);
  const Matrix &input = mTerminalInputScratch;
  if (!mInternalIndices.empty()) {
    if (mUsePrecomputedInternalMaps) {
      if (mUseFullyReducedInternalMaps)
        mInternalVoltageScratch.noalias() = mInternalStateMap * mState->get();
      else if (!mInternalStateResponseReady)
        throw std::logic_error(
            "Precomputed internal state response missing before post-step.");
      mInternalVoltageScratch.noalias() -= mInternalInputMapSparse * input;
    } else {
      mInternalHistoryScratch.noalias() = mHi * mState->get();
      mInternalHistoryScratch.noalias() -= mGin * input;
      mInternalVoltageScratch = mInternalFactor.solve(mInternalHistoryScratch);
      ++mInternalSolveCalls;
    }
  }
  for (UInt i = 0; i < mTerminalIndices.size(); ++i)
    mLocalVoltageScratch(mTerminalIndices[i]) = input(i);
  for (UInt i = 0; i < mInternalIndices.size(); ++i)
    mLocalVoltageScratch(mInternalIndices[i]) = mInternalVoltageScratch(i);
  for (UInt node = 0; node < mNodeVoltages.size(); ++node)
    if (mNodeVoltages[node])
      *mNodeVoltages[node] = mLocalVoltageScratch.block<3, 1>(3 * node, 0);
  if (mDetailedProfiling) {
    const auto now = Clock::now();
    mVoltageRecoveryTime +=
        std::chrono::duration_cast<std::chrono::nanoseconds>(now - start);
    start = now;
  }
  ++mVoltageRecoveryCalls;

  // Every row remains an original component companion state/history
  // coordinate. No rank reduction or state-coordinate transformation occurs.
  if (mUseFullyReducedInternalMaps) {
    mStateUpdateScratch.noalias() = mReducedStateMap * mState->get();
    mStateUpdateScratch.noalias() += mReducedInputMap * input;
  } else {
    mStateUpdateScratch.noalias() = mF * mState->get();
    if (!mInternalIndices.empty())
      mStateUpdateScratch.noalias() += mBi * mInternalVoltageScratch;
    mStateUpdateScratch.noalias() += mBn * input;
  }
  mState->get().swap(mStateUpdateScratch);
  mInternalStateResponseReady = false;
  if (mDetailedProfiling) {
    const auto now = Clock::now();
    mStateUpdateTime +=
        std::chrono::duration_cast<std::chrono::nanoseconds>(now - start);
    start = now;
  }
  ++mStateUpdateCalls;

  publishPhysicalOutputs();
  if (mDetailedProfiling)
    mOutputPublicationTime += std::chrono::duration_cast<std::chrono::nanoseconds>(
        Clock::now() - start);
  ++mOutputPublicationCalls;
}

void EMT::LinearSSNGroup::publishPhysicalOutputs() {
  const Matrix &state = mState->get();
  for (UInt index = 0; index < mModels.size(); ++index) {
    const auto &model = mModels[index];
    const auto &outputs = mModelOutputs[index];
    const auto &indices = mModelInputIndices[index];
    const UInt offset = mStateOffsets[index];
    PhaseVector voltage0 = PhaseVector::Zero();
    PhaseVector voltage1 = PhaseVector::Zero();
    for (UInt phase = 0; phase < 3; ++phase) {
      if (!model.grounded[0])
        voltage0(phase) = mLocalVoltageScratch(indices[phase]);
      if (model.terminalCount == 2 && !model.grounded[1])
        voltage1(phase) = mLocalVoltageScratch(indices[3 + phase]);
    }
    if (model.kind == LinearSSNModel::Kind::PiLine) {
      const PhaseVector current = state.block<3, 1>(offset, 0);
      *outputs.interfaceVoltage = voltage1 - voltage0;
      *outputs.interfaceCurrent = current;
      // Both currents point into the physical PiLine and include its end shunts.
      *outputs.terminalCurrent[0] =
          -current + 0.5 * model.parallelConductance * voltage0 +
          state.block<3, 1>(offset + 9, 0);
      *outputs.terminalCurrent[1] =
          current + 0.5 * model.parallelConductance * voltage1 +
          state.block<3, 1>(offset + 15, 0);
    } else if (model.kind == LinearSSNModel::Kind::ImpedanceLoad) {
      PhaseVector current = model.totalConductance * voltage0;
      if (model.hasLoadInductance)
        current += state.block<3, 1>(offset + 3, 0);
      if (model.hasLoadCapacitance)
        current += state.block<3, 1>(offset + 9, 0);
      *outputs.interfaceVoltage = voltage0;
      *outputs.interfaceCurrent = -current;
      *outputs.terminalCurrent[0] = current;
    } else {
      const PhaseVector current0 = state.block<3, 1>(offset + 3, 0);
      const PhaseVector current1 = state.block<3, 1>(offset + 9, 0);
      *outputs.interfaceVoltage = voltage1 - voltage0;
      *outputs.interfaceCurrent = model.orientationSign * current0;
      *outputs.terminalCurrent[model.referenceTerminal] = -current0;
      *outputs.terminalCurrent[1 - model.referenceTerminal] =
          model.ratioFromReference * current1;
    }
    ++*outputs.updateCount;
  }
}

void EMT::LinearSSNGroup::mnaCompAddPreStepDependencies(
    AttributeBase::List &prevStepDependencies, AttributeBase::List &,
    AttributeBase::List &modifiedAttributes) {
  prevStepDependencies.push_back(mState);
  modifiedAttributes.push_back(mRightVector);
}

void EMT::LinearSSNGroup::mnaCompAddPostStepDependencies(
    AttributeBase::List &, AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes,
    Attribute<Matrix>::Ptr &leftVector) {
  attributeDependencies.push_back(leftVector);
  modifiedAttributes.push_back(mState);
  // MNA task dependencies are collected before mnaCompInitialize builds the
  // realization. Enumerate constructor-owned attributes here, not cached maps.
  for (const auto &node : mAllNodes)
    if (publishesNodeVoltage(node))
      modifiedAttributes.push_back(node->mVoltage);
  for (const auto &component : mComponents) {
    modifiedAttributes.push_back(component->attribute("v_intf"));
    modifiedAttributes.push_back(component->attribute("i_intf"));
    modifiedAttributes.push_back(component->attribute("output_updates"));
    for (UInt terminal = 0; terminal < component->terminalNumber(); ++terminal)
      modifiedAttributes.push_back(
          component->attribute("i_terminal_" + std::to_string(terminal)));
  }
}
