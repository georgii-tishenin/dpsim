// SPDX-License-Identifier: MPL-2.0
#pragma once

#include <array>
#include <chrono>
#include <cstdint>
#include <ostream>
#include <string>

namespace CPS::EMT {

/// Optional exclusive clocks for passive-component kernels. Disabled scopes
/// never read the clock. Rows are nested in scheduler/solver timings.
class PassiveOperationProfile {
public:
  enum Operation {
    History,
    RhsLocal,
    VoltageRecovery,
    StateAdvance,
    VoltageState,
    Output,
    OperationCount
  };
  bool enabled = false;

  class Scope {
  public:
    Scope(PassiveOperationProfile &profile, Operation operation)
        : mProfile(profile), mOperation(operation) {
      if (mProfile.enabled)
        mStart = Clock::now();
    }
    ~Scope() { finish(); }
    void next(Operation operation) {
      finish();
      mOperation = operation;
      if (mProfile.enabled)
        mStart = Clock::now();
    }

  private:
    using Clock = std::chrono::steady_clock;
    PassiveOperationProfile &mProfile;
    Operation mOperation;
    Clock::time_point mStart{};
    void finish() {
      if (mProfile.enabled) {
        mProfile.mTime[mOperation] +=
            std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() -
                                                                 mStart)
                .count();
        ++mProfile.mCalls[mOperation];
      }
    }
  };

  void write(std::ostream &stream, const std::string &name) const {
    static constexpr const char *labels[] = {
        "history", "rhs_local", "voltage", "state", "voltage_state", "output"};
    for (unsigned i = 0; i < OperationCount; ++i)
      if (mCalls[i])
        stream << "passive_" << labels[i] << "_" << name << ','
               << mTime[i] * 1e-9 << ',' << mTime[i] * 1e-3 / mCalls[i]
               << ",0,0," << mCalls[i] << '\n';
  }

private:
  std::array<std::int64_t, OperationCount> mTime{};
  std::array<std::uint64_t, OperationCount> mCalls{};
};
} // namespace CPS::EMT
