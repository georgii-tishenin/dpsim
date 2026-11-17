// SPDX-License-Identifier: MPL-2.0
#pragma once

#include "IEEE_BenchmarkData.h"
#include <array>
#include <limits>

// Adapted from CPS::CIM::Examples::Grids::IEEE9::ScenarioConfig in
// dpsim/examples/cxx/Examples.h and EMT_Ph3_IEEE9_4Order.cpp. This retains that
// checkout's six line parameter sets, three step-up transformers, three loads,
// generator dispatch/voltage setpoints, and 100 MVA generator nameplate bases.
// It is not claimed to be a byte-for-byte copy of MATPOWER case9. No external
// data provenance is specified by ScenarioConfig. All branch impedances use
// the 100 MVA / 230 kV reference base, with transformer taps relative to the
// nominal winding-voltage ratio. Dynamic machine parameters/controllers from
// the original example are replaced by the common benchmark machine model.
// Pmax is not supplied by the source: zero means unspecified, so the common
// max(mbase,pmax) dynamic-base rule uses each generator's source rated base.
namespace IEEE9BenchmarkData {
using IEEEBenchmarkData::Branch;
using IEEEBenchmarkData::Bus;
using IEEEBenchmarkData::Generator;
inline constexpr double unlimited = std::numeric_limits<double>::infinity();

inline constexpr std::array<Bus, 9> buses = {{
    {1, 3, 0, 0, 0, 0, 1.04, 0, 16.5},
    {2, 2, 0, 0, 0, 0, 1.025, 0, 18},
    {3, 2, 0, 0, 0, 0, 1.025, 0, 13.8},
    {4, 1, 0, 0, 0, 0, 1, 0, 230},
    {5, 1, 125, 50, 0, 0, 1, 0, 230},
    {6, 1, 90, 30, 0, 0, 1, 0, 230},
    {7, 1, 0, 0, 0, 0, 1, 0, 230},
    {8, 1, 100, 35, 0, 0, 1, 0, 230},
    {9, 1, 0, 0, 0, 0, 1, 0, 230},
}};
inline constexpr std::array<Generator, 3> generators = {{
    {1, 71.6, 27, unlimited, -unlimited, 1.04, 100, 0},
    {2, 163, 6.7, unlimited, -unlimited, 1.025, 100, 0},
    {3, 85, -10.9, unlimited, -unlimited, 1.025, 100, 0},
}};
inline constexpr std::array<Branch, 9> branches = {{
    {5, 4, 0.01005, 0.08521, 1.0 / 5.688921, 0, 0, 0},
    {6, 4, 0.017083, 0.092216, 1.0 / 6.336801, 0, 0, 0},
    {7, 5, 0.032533, 0.162281, 1.0 / 3.28151, 0, 0, 0},
    {9, 6, 0.039806, 0.171651, 1.0 / 2.807618, 0, 0, 0},
    {7, 8, 0.00853, 0.072127, 1.0 / 6.717421, 0, 0, 0},
    {8, 9, 0.011984, 0.10115, 1.0 / 4.793121, 0, 0, 0},
    {1, 4, 0, 0.0576, 0, 100, 1, 0},
    {2, 7, 0, 0.0625, 0, 100, 1, 0},
    {3, 9, 0, 0.0586, 0, 100, 1, 0},
}};
} // namespace IEEE9BenchmarkData
