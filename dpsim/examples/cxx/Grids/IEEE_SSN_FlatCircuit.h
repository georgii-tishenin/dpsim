// SPDX-License-Identifier: MPL-2.0
#pragma once

#include "IEEE_SSN_IdealTransformer.h"
#include <algorithm>
#include <dpsim-models/CompositePowerComp.h>
#include <dpsim-models/EMT/EMT_Ph3_Transformer.h>
#include <dpsim-models/SystemTopology.h>

namespace IEEEBenchmark {

// Keep the established equipment output formulas, without executing the
// composite pre/post dispatch or its intermediate RHS aggregation. Children
// are ordinary top-level MNA components; this task only reconstructs outputs.
class FlatEquipmentOutput final : public CPS::MNASimPowerComp<CPS::Real> {
  std::shared_ptr<CPS::CompositePowerComp<CPS::Real>> mEquipment;

public:
  explicit FlatEquipmentOutput(
      std::shared_ptr<CPS::CompositePowerComp<CPS::Real>> equipment)
      : MNASimPowerComp(equipment->name() + ".Outputs", false, true),
        mEquipment(std::move(equipment)) {
    mPhaseType = CPS::PhaseType::ABC;
    setTerminalNumber(0);
    // Internal voltages and transformer current unknowns stay virtual and in
    // the same order. They must not become extra propagated network voltages.
    mVirtualNodes = mEquipment->virtualNodes();
    mNumVirtualNodes = static_cast<CPS::UInt>(mVirtualNodes.size());
  }

  void mnaCompInitialize(CPS::Real, CPS::Real,
                         CPS::Attribute<CPS::Matrix>::Ptr) override {
    mEquipment->updateMatrixNodeIndices();
    (**mRightVector).resize(0, 0);
  }

  void mnaCompPostStep(CPS::Real time, CPS::Int step,
                       CPS::Attribute<CPS::Matrix>::Ptr &left) override {
    mEquipment->mnaParentPostStep(time, step, left);
  }

  void mnaCompAddPostStepDependencies(
      CPS::AttributeBase::List &previous,
      CPS::AttributeBase::List &dependencies,
      CPS::AttributeBase::List &modified,
      CPS::Attribute<CPS::Matrix>::Ptr &left) override {
    mEquipment->mnaParentAddPostStepDependencies(previous, dependencies,
                                                 modified, left);
    // Preserve every primitive update, including stateless resistors whose
    // current may duplicate a series inductor's current.
    for (const auto &child : mEquipment->subComponents()) {
      dependencies.push_back(child->mIntfVoltage);
      dependencies.push_back(child->mIntfCurrent);
    }
  }
};

// Use the native circuit builders only for parameterization, internal-node
// phasors and output records. No CompositePowerComp remains in the topology.
inline void flattenPassiveCircuit(CPS::SystemTopology &system,
                                  const CPS::IdentifiedObject::List &passive) {
  CPS::IdentifiedObject::List components;
  for (const auto &component : system.mComponents) {
    if (std::find(passive.begin(), passive.end(), component) == passive.end()) {
      components.push_back(component);
      continue;
    }
    auto equipment =
        std::dynamic_pointer_cast<CPS::CompositePowerComp<CPS::Real>>(
            component);
    if (!equipment)
      throw std::invalid_argument("Flat circuit requires primitive equipment.");
    equipment->initializeFromNodesAndTerminals(system.mSystemFrequency);
    components.push_back(std::make_shared<FlatEquipmentOutput>(equipment));

    if (auto transformer =
            std::dynamic_pointer_cast<CPS::EMT::Ph3::Transformer>(equipment)) {
      // Only the reference external winding touches an R/L child; the other
      // winding connects through the ideal constraint. This also respects
      // explicitly selected reference windings, not only voltage-base order.
      CPS::UInt reference = 2;
      for (CPS::UInt terminal = 0; terminal < 2; ++terminal)
        for (const auto &child : transformer->subComponents())
          for (CPS::UInt port = 0; port < child->terminalNumberConnected();
               ++port)
            if (!child->node(port)->isGround() &&
                child->node(port) == transformer->node(terminal))
              reference = terminal;
      if (reference > 1)
        throw std::invalid_argument("Cannot identify transformer winding.");
      const CPS::Complex ratio = reference == 0
                                     ? transformer->effectiveRatio()
                                     : 1. / transformer->effectiveRatio();
      components.push_back(std::make_shared<IdealTransformer>(
          transformer->name() + ".Ideal", transformer->virtualNode(0),
          transformer->virtualNode(1), transformer->node(1 - reference),
          ratio.real(), transformer->terminalNotGrounded(reference)));
    }
    for (const auto &child : equipment->subComponents()) {
      if (child->hasSubComponents() || child->hasVirtualNodes())
        throw std::invalid_argument(
            "Expected elementary R/L/C circuit children.");
      components.push_back(child);
    }
  }
  system.mComponents.clear();
  system.addComponents(components);
}

} // namespace IEEEBenchmark
