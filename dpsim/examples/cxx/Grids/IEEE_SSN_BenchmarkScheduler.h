// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include <dpsim/Scheduler.h>

#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <vector>

namespace IEEEBenchmark {

// Keep the ordinary sequential schedule. Optional clocks surround contiguous
// blocks of the same mathematical work, rather than every component task.
class BenchmarkScheduler final : public DPsim::Scheduler {
public:
  explicit BenchmarkScheduler(bool profiling = false)
      : Scheduler(CPS::Logger::Level::off), mProfiling(profiling) {}

  void createSchedule(const CPS::Task::List &tasks, const Edges &inEdges,
                      const Edges &outEdges) override {
    topologicalSort(tasks, inEdges, outEdges, mSchedule);
    mBlocks.clear();
    mMeasurements = {};
    for (std::size_t i = 0; i < mSchedule.size(); ++i) {
      const auto category = classify(mSchedule[i]->toString());
      ++mMeasurements[category].taskCount;
      if (mBlocks.empty() || mBlocks.back().category != category)
        mBlocks.push_back({category, i, i + 1});
      else
        mBlocks.back().end = i + 1;
    }
  }

  void step(CPS::Real time, CPS::Int timeStepCount) override {
    if (!mProfiling) {
      for (const auto &task : mSchedule)
        task->executeFast(time, timeStepCount);
      return;
    }

    for (const auto &block : mBlocks) {
      const auto start = Clock::now();
      for (std::size_t i = block.begin; i < block.end; ++i)
        mSchedule[i]->executeFast(time, timeStepCount);
      const auto end = Clock::now();
      auto &measurement = mMeasurements[block.category];
      measurement.total +=
          std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
      ++measurement.blocks;
      measurement.calls += block.end - block.begin;
    }
  }

  // These are scheduled tasks and timed blocks per timestep, not accumulated
  // call counts. The CSV contains accumulated calls for each category.
  std::size_t executedTaskCount() const { return mSchedule.size(); }
  std::size_t blockCount() const { return mBlocks.size(); }

  void writeCSV(const std::string &filename) const {
    std::ofstream output(filename);
    if (!output)
      throw std::runtime_error("Cannot write scheduler measurements: " +
                               filename);
    output << "category,total_s,blocks,calls,task_count\n"
           << std::setprecision(17);
    static constexpr std::array<const char *, CategoryCount> labels = {
        "passive_history", "passive_update", "generators", "mna", "other"};
    for (std::size_t i = 0; i < CategoryCount; ++i) {
      const auto &measurement = mMeasurements[i];
      output << labels[i] << ','
             << std::chrono::duration<double>(measurement.total).count() << ','
             << measurement.blocks << ',' << measurement.calls << ','
             << measurement.taskCount << '\n';
    }
    if (!output)
      throw std::runtime_error("Failed writing scheduler measurements: " +
                               filename);
  }

private:
  using Clock = std::chrono::steady_clock;
  enum Category {
    PassiveHistory,
    PassiveUpdate,
    Generators,
    Mna,
    Other,
    CategoryCount
  };
  struct Block {
    Category category;
    std::size_t begin, end;
  };
  struct Measurement {
    std::chrono::nanoseconds total{};
    std::uint64_t blocks = 0, calls = 0;
    std::size_t taskCount = 0;
  };

  static bool endsWith(const std::string &name, const std::string &suffix) {
    return name.size() >= suffix.size() &&
           name.compare(name.size() - suffix.size(), suffix.size(), suffix) ==
               0;
  }

  static Category classify(const std::string &name) {
    if (endsWith(name, ".Solve"))
      return Mna;
    if (name.rfind("Gen", 0) == 0)
      return Generators;
    const bool passive =
        name.rfind("Load", 0) == 0 || name.rfind("Shunt", 0) == 0 ||
        name.rfind("Line", 0) == 0 || name.rfind("Trafo", 0) == 0 ||
        name.rfind("PassiveGroup", 0) == 0;
    if (passive && endsWith(name, ".MnaPreStep"))
      return PassiveHistory;
    if (passive && endsWith(name, ".MnaPostStep"))
      return PassiveUpdate;
    return Other;
  }

  const bool mProfiling;
  CPS::Task::List mSchedule;
  std::vector<Block> mBlocks;
  std::array<Measurement, CategoryCount> mMeasurements{};
};

} // namespace IEEEBenchmark
