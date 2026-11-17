// SPDX-License-Identifier: MPL-2.0
#pragma once

// The same power-flow fields/units as IEEE118_Data.h and IEEE300_Data.h.
// This small schema is used to adapt the existing IEEE9 and IEEE39 examples.
namespace IEEEBenchmarkData {
struct Bus {
  int id, type;
  double pd, qd, gs, bs, vm, va, baseKV;
};
struct Generator {
  int bus;
  double pg, qg, qmax, qmin, vg, mbase, pmax;
};
struct Branch {
  int from, to;
  double r, x, b, rateA, ratio, angle;
};
} // namespace IEEEBenchmarkData
