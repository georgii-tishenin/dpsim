// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0
#pragma once

#include <cmath>
#include <dpsim-models/MNASimPowerComp.h>
#include <stdexcept>

namespace IEEEBenchmark {

/// Static turns-ratio constraint for the flat transformer T circuit.
/// Reuses its referred-voltage node and auxiliary current unknown: no states,
/// history source, or per-step tasks. It can own the original circuit's virtual
/// nodes in their original order, preserving the MNA unknown numbering.
class IdealTransformer : public CPS::MNASimPowerComp<CPS::Real> {
  CPS::EMT::SimNode::Ptr mCurrentUnknown;
  CPS::Real mRatio;
  CPS::Bool mReferenceConnected;

public:
  IdealTransformer(CPS::String name, CPS::EMT::SimNode::Ptr referredNode,
                   CPS::EMT::SimNode::Ptr currentUnknown,
                   CPS::EMT::SimNode::Ptr otherWindingNode,
                   CPS::Real ratioFromReference,
                   CPS::Bool referenceConnected = true,
                   const CPS::EMT::SimNode::List &virtualNodes = {})
      : MNASimPowerComp(name, false, false), mCurrentUnknown(currentUnknown),
        mRatio(ratioFromReference), mReferenceConnected(referenceConnected) {
    if (!referredNode || !currentUnknown || !otherWindingNode ||
        currentUnknown->isGround() || !std::isfinite(mRatio) || mRatio <= 0.)
      throw std::invalid_argument("Invalid ideal-transformer nodes or ratio.");
    mPhaseType = CPS::PhaseType::ABC;
    mVirtualNodes = virtualNodes;
    mNumVirtualNodes = static_cast<CPS::UInt>(mVirtualNodes.size());
    setTerminalNumber(2);
    connect({referredNode, otherWindingNode});
  }

  void mnaCompInitialize(CPS::Real, CPS::Real,
                         CPS::Attribute<CPS::Matrix>::Ptr) override {
    // An empty stamp also excludes this source-free constraint from RHS sums.
    (**mRightVector).resize(0, 0);
  }

  void
  mnaCompApplySystemMatrixStamp(CPS::SparseMatrixRow &systemMatrix) override {
    updateMatrixNodeIndices();
    for (CPS::UInt phase = 0; phase < 3; ++phase) {
      const auto current =
          mCurrentUnknown->matrixNodeIndex(static_cast<CPS::PhaseType>(phase));
      // v_referred - ratio * v_other = 0; the winding currents are
      // -i_aux and ratio * i_aux, respectively (the native stamp convention).
      if (mReferenceConnected && terminalNotGrounded(0)) {
        const auto referred = matrixNodeIndex(0, phase);
        CPS::Math::setMatrixElement(systemMatrix, referred, current, -1.);
        CPS::Math::setMatrixElement(systemMatrix, current, referred, 1.);
      }
      if (terminalNotGrounded(1)) {
        const auto other = matrixNodeIndex(1, phase);
        CPS::Math::setMatrixElement(systemMatrix, other, current, mRatio);
        CPS::Math::setMatrixElement(systemMatrix, current, other, -mRatio);
      }
    }
  }
};

} // namespace IEEEBenchmark
