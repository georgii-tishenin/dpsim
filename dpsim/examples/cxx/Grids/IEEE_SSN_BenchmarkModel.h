// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0
#pragma once

#include "IEEE_SSN_FlatCircuit.h"

// Shared electrical construction for all benchmark grids and representations.
// Each .cpp selects its grid before including this shared model.
#ifndef IEEE_SSN_BUS_COUNT
#define IEEE_SSN_BUS_COUNT 118
#endif
#include "IEEE118_Data.h"
#if IEEE_SSN_BUS_COUNT == 9
#include "IEEE9_BenchmarkData.h"
#elif IEEE_SSN_BUS_COUNT == 39
#include "IEEE39_BenchmarkData.h"
#elif IEEE_SSN_BUS_COUNT == 300
#include "IEEE300_Data.h"
#elif IEEE_SSN_BUS_COUNT != 118
#error "IEEE_SSN_BUS_COUNT must be 9, 39, 118, or 300"
#endif

#include <DPsim.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <dpsim-models/EMT/EMT_LinearSSNGroup.h>
#include <dpsim-models/EMT/EMT_Ph3_SSN_ImpedanceLoad.h>
#include <dpsim-models/EMT/EMT_Ph3_SSN_PiLine.h>
#include <dpsim-models/EMT/EMT_Ph3_SSN_Transformer.h>
#include <dpsim-models/EMT/EMT_Ph3_SynchronGenerator4OrderVBR.h>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace IEEEBenchmark {
using namespace CPS;
using namespace DPsim;
#if IEEE_SSN_BUS_COUNT == 9
namespace Data = IEEE9BenchmarkData;
inline constexpr char gridName[] = "IEEE9";
inline constexpr int eventBus = 2;
#elif IEEE_SSN_BUS_COUNT == 39
namespace Data = IEEE39BenchmarkData;
inline constexpr char gridName[] = "IEEE39";
inline constexpr int eventBus = 31;
#elif IEEE_SSN_BUS_COUNT == 300
namespace Data = IEEE300Data;
inline constexpr char gridName[] = "IEEE300";
// Largest positive generator dispatch among buses rated at least 100 kV.
inline constexpr int eventBus = 170;
#else
namespace Data = IEEE118Data;
inline constexpr char gridName[] = "IEEE118";
inline constexpr int eventBus = 59;
#endif
using namespace Data;
inline constexpr int gridBusCount = IEEE_SSN_BUS_COUNT;
// A negative demand in a power-flow case specifies injection, not a physical
// negative resistor. Keep that injection prescribed at its solved PF phasor.
// The signed series equivalent at IEEE300 bus245 also requires its exact zero
// shunt conductance; adding numerical shunts creates an artificial local mode.
#if IEEE_SSN_BUS_COUNT == 300
inline constexpr char modelAdaptation[] =
    "prescribed_negative_load_currents_exact_series_equivalent";
#else
inline constexpr char modelAdaptation[] = "none";
#endif
inline constexpr Real frequency = 60.;
inline constexpr Real omega = 2. * PI * frequency;
inline constexpr Real baseMVA = 100.;
// Generic machine parameters, deliberately identical in every grid. Only
// the source nameplate/dispatch base and bus voltage change. No controllers.
inline constexpr Real H = 4., Xd = 1.75, Xq = 1.75, Xl = .15;
inline constexpr Real XdTransient = .425, XqTransient = .425;
inline constexpr Real Td0Transient = 5., Tq0Transient = .5;

inline Real baseV(const Bus &bus) { return bus.baseKV * 1e3; }
inline Real zBase(const Bus &bus) {
  return baseV(bus) * baseV(bus) / (baseMVA * 1e6);
}
inline const Bus &busFor(int id) {
  return buses.at(static_cast<std::size_t>(id - 1));
}
inline Real dynamicMachineBase(const Generator &gen) {
  return std::max(gen.mbase, gen.pmax) * 1e6;
}
inline Matrix parameter3(Real value) {
  return Math::singlePhaseParameterToThreePhase(value);
}
// Signed sum of source Pd values: existing negative demands are retained.
inline constexpr Real totalLoadPower() {
  Real total = 0.;
  for (const auto &bus : buses)
    total += bus.pd * 1e6;
  return total;
}
inline constexpr Real loadStepPower() { return .20 * totalLoadPower(); }

inline PowerflowBusType pfBusType(int type) {
  if (type == 3)
    return PowerflowBusType::VD;
  if (type == 2)
    return PowerflowBusType::PV;
  return PowerflowBusType::PQ;
}

inline SystemTopology makePowerflowSystem() {
  SystemNodeList nodes;
  SystemComponentList components;
  std::vector<SP::SimNode::Ptr> busNodes;
  busNodes.reserve(buses.size());
  for (const auto &bus : buses) {
    auto node =
        SP::SimNode::make("N" + std::to_string(bus.id), PhaseType::Single);
    busNodes.push_back(node);
    nodes.push_back(node);
  }

  for (const auto &gen : generators) {
    const auto &bus = busFor(gen.bus);
    auto machine =
        SP::Ph1::SynchronGenerator::make("Gen" + std::to_string(gen.bus));
    machine->setParameters(gen.mbase * 1e6, baseV(bus), gen.pg * 1e6,
                           gen.vg * baseV(bus), pfBusType(bus.type),
                           gen.qg * 1e6, gen.qmax * 1e6, gen.qmin * 1e6);
    machine->setBaseVoltage(baseV(bus));
    machine->connect({busNodes.at(gen.bus - 1)});
    components.push_back(machine);
  }

  int loadIndex = 0;
  for (const auto &bus : buses) {
    if (bus.pd != 0. || bus.qd != 0.) {
      auto load = SP::Ph1::Load::make("PFLoad" + std::to_string(++loadIndex));
      load->setParameters(bus.pd * 1e6, bus.qd * 1e6, baseV(bus));
      load->connect({busNodes.at(bus.id - 1)});
      components.push_back(load);
    }
    if (bus.gs != 0. || bus.bs != 0.) {
      auto shunt = SP::Ph1::Shunt::make("PFShunt" + std::to_string(bus.id));
      shunt->setParameters(bus.gs * 1e6 / (baseV(bus) * baseV(bus)),
                           bus.bs * 1e6 / (baseV(bus) * baseV(bus)));
      shunt->setBaseVoltage(baseV(bus));
      shunt->connect({busNodes.at(bus.id - 1)});
      components.push_back(shunt);
    }
  }

  int lineIndex = 0;
  int transformerIndex = 0;
  for (const auto &branch : branches) {
    const auto &from = busFor(branch.from);
    const auto &to = busFor(branch.to);
    const Real vf = baseV(from), vt = baseV(to);
    const Real tap = branch.ratio == 0. ? 1. : branch.ratio;
    const Real ratio = tap * vf / vt;
    const Real zTo = zBase(to);
    if (branch.ratio == 0. && vf == vt) {
      auto line = SP::Ph1::PiLine::make("PFLine" + std::to_string(++lineIndex));
      line->setParameters(branch.r * zBase(to), branch.x * zBase(to) / omega,
                          branch.b / zBase(to) / omega, 0.);
      line->setBaseVoltage(vt);
      line->connect({busNodes.at(branch.from - 1), busNodes.at(branch.to - 1)});
      components.push_back(line);
    } else {
      auto transformer = SP::Ph1::Transformer::make(
          "PFTrafo" + std::to_string(++transformerIndex));
      const Real highV = std::max(vf, vt);
      const Real transformerBaseZ = zTo;
      const Real transformerScale = (vf > vt) ? ratio : tap;
      const Real transformerR =
          branch.r * transformerBaseZ * transformerScale * transformerScale;
      const Real transformerL = branch.x * transformerBaseZ * transformerScale *
                                transformerScale / omega;
      const Real rate = branch.rateA > 0. ? branch.rateA * 1e6 : baseMVA * 1e6;
      transformer->setParameters(vf, vt, rate, std::abs(ratio),
                                 branch.angle * PI / 180., transformerR,
                                 transformerL);
      transformer->setBaseVoltage(highV);
      transformer->connect(
          {busNodes.at(branch.from - 1), busNodes.at(branch.to - 1)});
      components.push_back(transformer);
    }
  }
  return SystemTopology(frequency, nodes, components);
}

// Preserve the benchmark's electrical conversion exactly: series impedance
// is referred to the appropriate transformer winding, and no resistance or
// reactance is clipped. Native branches retain their primitive subcomponents.
inline void addBranches(const std::vector<EMT::SimNode::Ptr> &nodes,
                        SystemComponentList &components, bool ssn,
                        bool elementary) {
  int lines = 0, transformers = 0;
  for (const auto &branch : branches) {
    const auto &from = busFor(branch.from), &to = busFor(branch.to);
    const Real vf = baseV(from), vt = baseV(to);
    const Real tap = branch.ratio == 0. ? 1. : branch.ratio;
    const Real ratio = tap * vf / vt, zTo = zBase(to);
    if (branch.ratio == 0. && vf == vt) {
      const String name = "Line" + std::to_string(++lines);
      const Real r = branch.r * zTo, l = branch.x * zTo / omega;
      const Real c = branch.b / zTo / omega;
      const bool exactSeriesEquivalent =
          gridBusCount == 300 && ((branch.from == 97 && branch.to == 245) ||
                                  (branch.from == 245 && branch.to == 99));
      if (ssn) {
        auto line = EMT::Ph3::SSN::PiLine::make(name);
        line->setParameters(parameter3(r), parameter3(l), parameter3(c),
                            parameter3(0.));
        if (exactSeriesEquivalent)
          line->setUseDefaultShuntConductance(false);
        line->connect({nodes.at(branch.from - 1), nodes.at(branch.to - 1)});
        components.push_back(line);
      } else {
        auto line = EMT::Ph3::PiLine::make(name);
        line->setUseElementarySSN(elementary);
        line->setParameters(parameter3(r), parameter3(l), parameter3(c),
                            parameter3(0.));
        if (exactSeriesEquivalent)
          line->setUseDefaultShuntConductance(false);
        line->connect({nodes.at(branch.from - 1), nodes.at(branch.to - 1)});
        components.push_back(line);
      }
    } else {
      const String name = "Trafo" + std::to_string(++transformers);
      const Real scale = vf > vt ? ratio : tap;
      const Real r = branch.r * zTo * scale * scale;
      const Real l = branch.x * zTo * scale * scale / omega;
      const Real rate = branch.rateA > 0. ? branch.rateA * 1e6 : baseMVA * 1e6;
      if (ssn) {
        auto transformer = EMT::Ph3::SSN::Transformer::make(name);
        transformer->setParameters(vf, vt, rate, std::abs(ratio),
                                   branch.angle * PI / 180., parameter3(r),
                                   parameter3(l));
        transformer->connect(
            {nodes.at(branch.from - 1), nodes.at(branch.to - 1)});
        components.push_back(transformer);
      } else {
        auto transformer = EMT::Ph3::Transformer::make(name, true);
        transformer->setUseElementarySSN(elementary);
        transformer->setParameters(vf, vt, rate, std::abs(ratio),
                                   branch.angle * PI / 180., parameter3(r),
                                   parameter3(l));
        transformer->connect(
            {nodes.at(branch.from - 1), nodes.at(branch.to - 1)});
        components.push_back(transformer);
      }
    }
  }
}

struct Network {
  SystemTopology system{60.};
  SystemTopology powerflow{60.};
  std::vector<EMT::SimNode::Ptr> nodes;
  std::unordered_map<int, EMT::SimNode::Ptr> retained;
  SystemComponentList passive, machines, injections;
  EMT::LinearSSNGroup::Ptr group;
  String formulation, monitoredBranch;
  UInt monitoredTerminal = 0;
  int monitoredFrom = 0, monitoredTo = 0;
  std::size_t loadCount = 0, shuntCount = 0;
  std::size_t lineCount = 0, transformerCount = 0;
  std::size_t activeCoordinates = 0;
  double powerflowSeconds = 0.;

  EMT::SimNode::Ptr terminal(int bus) const {
    return retained.empty() ? nodes.at(bus - 1) : retained.at(bus);
  }
};

// The same circuit, machines and physical outputs at three SSN granularities.
// Native companions are the accuracy reference. A network group retains the
// machine/source buses and recovers every internal voltage and branch current.
inline Network makeNetwork(String formulation) {
  const bool native = formulation == "native";
  const bool elementary =
      formulation == "elementary" || formulation == "elementary_wrapped";
  const bool grouped = formulation == "network";
  if (!native && !elementary && !grouped && formulation != "equipment")
    throw std::invalid_argument("Expected native, elementary, "
                                "elementary_wrapped, equipment, or network");

  Network network;
  network.formulation = formulation;
  const auto powerflowStart = std::chrono::steady_clock::now();
  network.powerflow = makePowerflowSystem();
  Simulation pf(String(gridName) + "_Benchmark_PF", Logger::Level::off);
  pf.setSystem(network.powerflow);
  pf.setTimeStep(.01);
  pf.setFinalTime(.01);
  pf.setDomain(Domain::SP);
  pf.setSolverType(Solver::Type::NRP);
  pf.setSolverAndComponentBehaviour(Solver::Behaviour::Initialization);
  pf.doInitFromNodesAndTerminals(false);
  pf.run();
  network.powerflowSeconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                    powerflowStart)
          .count();

  SystemNodeList systemNodes;
  for (const auto &bus : buses) {
    auto node =
        EMT::SimNode::make("N" + std::to_string(bus.id), PhaseType::ABC);
    node->setInitialVoltage(
        network.powerflow.node<SP::SimNode>(node->name())->singleVoltage());
    network.nodes.push_back(node);
    if (!grouped)
      systemNodes.push_back(node);
  }
  std::vector<int> retainedBusIds;
  if (grouped) {
    const auto retain = [&](int bus) {
      if (network.retained.count(bus))
        return;
      auto node = EMT::SimNode::make("N" + std::to_string(bus), PhaseType::ABC);
      network.retained.emplace(bus, node);
      systemNodes.push_back(node);
      retainedBusIds.push_back(bus);
    };
    for (const auto &gen : generators)
      retain(gen.bus);
    if (gridBusCount == 300)
      for (const auto &bus : buses)
        if (bus.pd < 0.)
          retain(bus.id);
  }

  SystemComponentList components;
  for (const auto &gen : generators) {
    auto machine = EMT::Ph3::SynchronGenerator4OrderVBR::make(
        "Gen" + std::to_string(gen.bus));
    machine->setOperationalParametersPerUnit(
        dynamicMachineBase(gen), baseV(busFor(gen.bus)), frequency, H, Xd, Xq,
        Xl, XdTransient, XqTransient, Td0Transient, Tq0Transient);
    machine->connect({network.terminal(gen.bus)});
    network.machines.push_back(machine);
    components.push_back(machine);
  }

  // Keep loads and fixed shunts separate, as in the native network.
  for (const auto &bus : buses) {
    const auto &node = network.nodes.at(bus.id - 1);
    const bool prescribedInjection = gridBusCount == 300 && bus.pd < 0.;
    const Real loadActivePower = prescribedInjection ? 0. : bus.pd;
    if (prescribedInjection) {
      const Complex voltage =
          network.powerflow.node<SP::SimNode>("N" + std::to_string(bus.id))
              ->singleVoltage();
      // CurrentSource's positive current flows from its node toward ground.
      // Negative Pd therefore gives the correct injecting current phasor.
      const Complex current =
          std::conj(Complex(bus.pd * 1e6, 0.) / (std::sqrt(3.) * voltage));
      auto injection =
          EMT::Ph3::CurrentSource::make("Injection" + std::to_string(bus.id));
      injection->setParameters(Math::singlePhaseVariableToThreePhase(current),
                               frequency);
      injection->connect({EMT::SimNode::GND, network.terminal(bus.id)});
      network.injections.push_back(injection);
      components.push_back(injection);
    }
    if (loadActivePower != 0. || bus.qd != 0.) {
      const String name = "Load" + std::to_string(bus.id);
      if (native || elementary) {
        auto load = EMT::Ph3::RXLoad::make(name, true);
        load->setUseElementarySSN(elementary);
        load->setParameters(parameter3(loadActivePower * 1e6 / 3.),
                            parameter3(bus.qd * 1e6 / 3.), bus.vm * baseV(bus));
        load->connect({node});
        network.passive.push_back(load);
      } else {
        auto load = EMT::Ph3::SSN::ImpedanceLoad::make(name);
        load->setParameters(loadActivePower * 1e6, bus.qd * 1e6,
                            bus.vm * baseV(bus), 0., 0.);
        load->connect({node});
        network.passive.push_back(load);
      }
      ++network.loadCount;
      if (bus.qd != 0.)
        network.activeCoordinates += 6;
    }

    const Real g = bus.gs * 1e6 / std::pow(baseV(bus), 2);
    const Real b = bus.bs * 1e6 / std::pow(baseV(bus), 2);
    if (g != 0. || b != 0.) {
      const String name = "Shunt" + std::to_string(bus.id);
      if (native || elementary) {
        auto shunt = EMT::Ph3::Shunt::make(name);
        shunt->setUseElementarySSN(elementary);
        shunt->setParameters(g, b);
        shunt->connect({node});
        network.passive.push_back(shunt);
      } else {
        auto shunt = EMT::Ph3::SSN::ImpedanceLoad::make(name);
        shunt->setParameters(0., 0., bus.vm * baseV(bus), g, b);
        shunt->connect({node});
        network.passive.push_back(shunt);
      }
      ++network.shuntCount;
      if (b != 0.)
        network.activeCoordinates += 6;
    }
  }

  // Native uses the ordinary primitive composite branches; SSN uses the
  // corresponding companion-coordinate line and transformer realizations.
  addBranches(network.nodes, network.passive, !native && !elementary,
              elementary);
  for (const auto &branch : branches) {
    if (branch.ratio == 0. &&
        baseV(busFor(branch.from)) == baseV(busFor(branch.to))) {
      ++network.lineCount;
      if ((branch.from == eventBus || branch.to == eventBus) &&
          (network.monitoredBranch.empty() ||
           network.monitoredBranch.rfind("Trafo", 0) == 0)) {
        network.monitoredBranch = "Line" + std::to_string(network.lineCount);
        network.monitoredTerminal = branch.from == eventBus ? 0 : 1;
        network.monitoredFrom = branch.from;
        network.monitoredTo = branch.to;
      }
      network.activeCoordinates += 6 + (branch.b != 0. ? 12 : 0);
    } else {
      ++network.transformerCount;
      if ((branch.from == eventBus || branch.to == eventBus) &&
          network.monitoredBranch.empty()) {
        network.monitoredBranch =
            "Trafo" + std::to_string(network.transformerCount);
        network.monitoredTerminal = branch.from == eventBus ? 0 : 1;
        network.monitoredFrom = branch.from;
        network.monitoredTo = branch.to;
      }
      network.activeCoordinates += 18;
    }
  }

  if (grouped) {
    EMT::LinearSSNGroup::ComponentList children;
    for (const auto &component : network.passive)
      children.push_back(
          std::dynamic_pointer_cast<EMT::LinearSSNGroup::Component>(component));
    EMT::SimNode::List terminals;
    for (const int bus : retainedBusIds)
      terminals.push_back(network.retained.at(bus));
    network.group = EMT::LinearSSNGroup::make("PassiveGroup", children,
                                              network.nodes, terminals);
    network.group->setPrecomputedInternalMaps(true);
    network.group->setFullyReducedInternalMaps(false);
    components.push_back(network.group);
  } else
    components.insert(components.end(), network.passive.begin(),
                      network.passive.end());

  network.system = SystemTopology(frequency, systemNodes, components);
  network.system.initWithPowerflow(network.powerflow, Domain::EMT);
  if (formulation == "elementary")
    flattenPassiveCircuit(network.system, network.passive);
  return network;
}

// Call after Simulation::initialize(): the exported companion realization is
// then available. Grouping must preserve all allocated SSN coordinates.
inline std::size_t allocatedCoordinates(Network &network) {
  if (network.formulation == "native" || network.formulation == "elementary" ||
      network.formulation == "elementary_wrapped")
    return network.activeCoordinates;
  std::size_t count = 0;
  for (const auto &component : network.passive) {
    EMT::LinearSSNModel model;
    if (auto line = std::dynamic_pointer_cast<EMT::Ph3::SSN::PiLine>(component))
      line->exportFusedModel(model);
    else if (auto load =
                 std::dynamic_pointer_cast<EMT::Ph3::SSN::ImpedanceLoad>(
                     component))
      load->exportFusedModel(model);
    else if (auto transformer =
                 std::dynamic_pointer_cast<EMT::Ph3::SSN::Transformer>(
                     component))
      transformer->exportFusedModel(model);
    else
      throw std::runtime_error("Unsupported passive SSN component");
    count += model.state.rows();
  }
  if (network.group && count != network.group->stateCount())
    throw std::runtime_error("Group has changed companion coordinate count");
  return count;
}

} // namespace IEEEBenchmark
