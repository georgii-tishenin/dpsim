// SPDX-License-Identifier: MPL-2.0
#pragma once

#include "IEEE_SSN_BenchmarkModel.h"
#include "IEEE_SSN_BenchmarkScheduler.h"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <nlohmann/json.hpp>

// Shared output recording and optional operation measurements.
// Network construction is in BenchmarkModel; each .cpp configures a simulation.
namespace IEEEBenchmark {
using BenchmarkClock = std::chrono::steady_clock;
inline double elapsed(BenchmarkClock::time_point begin,
                      BenchmarkClock::time_point end = BenchmarkClock::now()) {
  return std::chrono::duration<double>(end - begin).count();
}
constexpr Real timeStep = 50e-6;

struct Options {
  String formulation, mode, event;
  size_t steps;
  std::filesystem::path directory;
  bool profile() const { return mode == "profile"; }
  bool validate() const { return mode == "validate"; }

  Options(int argc, char **argv, const String &caseName, size_t defaultSteps,
          const String &defaultMode = "timing",
          const String &defaultFormulation = "network",
          const String &defaultEvent = "")
      : formulation(argc > 1 ? argv[1] : defaultFormulation),
        mode(argc > 2 ? argv[2] : defaultMode),
        event(argc > 5 ? argv[5] : defaultEvent),
        steps(argc > 3 ? std::stoul(argv[3]) : defaultSteps),
        directory(argc > 4
                      ? argv[4]
                      : "logs/" + caseName + "/" + formulation + "_" + mode) {
    if (argc > (defaultEvent.empty() ? 5 : 6) ||
        (!event.empty() && event != "load" && event != "fault") ||
        (formulation != "native" && formulation != "elementary" &&
         formulation != "elementary_wrapped" && formulation != "equipment" &&
         formulation != "network") ||
        (mode != "timing" && mode != "profile" && mode != "validate") ||
        steps == 0 || steps > 10000000)
      throw std::invalid_argument(
          "Arguments: [native|elementary|elementary_wrapped|equipment|network] "
          "[timing|profile|validate] "
          "[steps] [output-directory] (accuracy also accepts "
          "[load|fault])");
    std::filesystem::create_directories(directory);
    Logger::setLogDir((directory / "logs").string());
    Eigen::setNbThreads(1);
  }
};

// An output consumer is required even without CSV logging: otherwise the
// dependency scheduler can discard output-only tasks (e.g. a resistive load).
class OutputDemand final : public DataLoggerInterface {
public:
  AttributeBase::List attributes;
  class DemandTask final : public Task {
  public:
    explicit DemandTask(const AttributeBase::List &attributes)
        : Task("PhysicalOutputDemand") {
      mAttributeDependencies = attributes;
      mModifiedAttributes.push_back(Scheduler::external);
    }
    void execute(Real, Int) override {}
  };
  void start() override {}
  void stop() override {}
  void log(Real, Int) override {}
  Task::Ptr getTask() override {
    return std::make_shared<DemandTask>(attributes);
  }
};

class Measurements {
  struct Block {
    String name;
    Attribute<Matrix>::Ptr matrix;
    Attribute<Real>::Ptr scalar;
    size_t offset, size;
  };
  Network &mNetwork;
  Options mOptions;
  String mCase, mSolver;
  std::shared_ptr<BenchmarkScheduler> mScheduler;
  std::shared_ptr<EMT::Ph3::Switch> mEventSwitch;
  std::vector<Block> mBlocks;
  std::vector<Real> mValues;
  std::ofstream mTrajectory, mWaveforms;
  size_t mPublications = 0;
  double mOutputSeconds = 0., mFirstClosedTime = -1., mFirstOpenTime = -1.;
  Attribute<Matrix>::Ptr mBranchCurrent;
  BenchmarkClock::time_point mPreparationStart, mInitializeStart, mInitializeEnd;
  BenchmarkClock::time_point mLoopStart, mLoopEnd;

  void add(const String &name, Attribute<Matrix>::Ptr value) {
    const size_t size = value->get().size();
    if (!size)
      throw std::runtime_error("Uninitialized output: " + name);
    mBlocks.push_back({name, value, nullptr, mValues.size(), size});
    mValues.resize(mValues.size() + size);
  }
  void add(const String &name, Attribute<Real>::Ptr value) {
    mBlocks.push_back({name, nullptr, value, mValues.size(), 1});
    mValues.resize(mValues.size() + 1);
  }
  void publish() {
    for (const auto &block : mBlocks) {
      if (block.matrix.getPtr())
        std::copy_n(block.matrix->get().data(), block.size,
                    mValues.data() + block.offset);
      else
        mValues[block.offset] = block.scalar->get();
    }
    ++mPublications;
  }
  void checkFinite() const {
    for (const Real value : mValues)
      if (!std::isfinite(value))
        throw std::runtime_error("Non-finite physical output");
  }
  void trajectoryRow(size_t step, Real time) {
    mTrajectory << step << ',' << time;
    for (Real value : mValues)
      mTrajectory << ',' << value;
    mTrajectory << '\n';
  }

public:
  Measurements(Simulation &sim, Network &network, const Options &options,
               const String &caseName, const String &solver,
               BenchmarkClock::time_point preparationStart,
               std::shared_ptr<EMT::Ph3::Switch> eventSwitch = nullptr)
      : mNetwork(network), mOptions(options), mCase(caseName), mSolver(solver),
        mScheduler(std::make_shared<BenchmarkScheduler>(options.profile())),
        mEventSwitch(eventSwitch), mPreparationStart(preparationStart) {
    // These benchmarks are connected. Subnet discovery assumes all terminal
    // nodes are physical buses, which does not hold for flat internal circuits.
    sim.doSplitSubnets(false);
    sim.setScheduler(mScheduler);
    sim.setSolverPhaseProfiling(options.profile());
    auto demand = std::make_shared<OutputDemand>();
    for (const auto &node : network.nodes)
      demand->attributes.push_back(node->mVoltage);
    for (const auto &machine : network.machines)
      for (const String attr : {"i_intf", "v_intf", "w_r", "Tm", "Te"})
        demand->attributes.push_back(machine->attribute(attr));
    for (const auto &source : network.injections)
      demand->attributes.push_back(source->attribute("i_intf"));
    for (const auto &component : network.passive) {
      demand->attributes.push_back(component->attribute("i_terminal_0"));
      if (component->name().find("Line") == 0 ||
          component->name().find("Trafo") == 0)
        demand->attributes.push_back(component->attribute("i_terminal_1"));
    }
    if (eventSwitch)
      for (const String attr : {"v_intf", "i_intf"})
        demand->attributes.push_back(eventSwitch->attribute(attr));
    sim.addLogger(demand);
    mInitializeStart = BenchmarkClock::now();
  }

  // Called immediately after Simulation.initialize(), before sim.start().
  void bindOutputs() {
    mInitializeEnd = BenchmarkClock::now();
    for (const auto &node : mNetwork.nodes)
      add(node->name() + ".v", node->mVoltage);
    for (const auto &machine : mNetwork.machines) {
      for (const String attr : {"i_intf", "v_intf"})
        add(machine->name() + "." + attr,
            machine->attributeTyped<Matrix>(attr));
      for (const String attr : {"w_r", "Tm", "Te"})
        add(machine->name() + "." + attr, machine->attributeTyped<Real>(attr));
    }
    for (const auto &component : mNetwork.passive) {
      add(component->name() + ".i_terminal_0",
          component->attributeTyped<Matrix>("i_terminal_0"));
      if (component->name().find("Line") == 0 ||
          component->name().find("Trafo") == 0)
        add(component->name() + ".i_terminal_1",
            component->attributeTyped<Matrix>("i_terminal_1"));
    }
    for (const auto &source : mNetwork.injections)
      add(source->name() + ".i_intf", source->attributeTyped<Matrix>("i_intf"));
    if (mEventSwitch) {
      for (const String attr : {"v_intf", "i_intf"})
        add("Event." + attr, mEventSwitch->attributeTyped<Matrix>(attr));
      for (const auto &component : mNetwork.passive)
        if (component->name() == mNetwork.monitoredBranch)
          mBranchCurrent = component->attributeTyped<Matrix>(
              "i_terminal_" + std::to_string(mNetwork.monitoredTerminal));
      if (!mBranchCurrent.getPtr())
        throw std::runtime_error("Missing monitored branch current");
    }
    const size_t expectedValues =
        3 * mNetwork.nodes.size() + 9 * mNetwork.machines.size() +
        3 * mNetwork.injections.size() +
        6 * (mNetwork.lineCount + mNetwork.transformerCount) +
        3 * (mNetwork.loadCount + mNetwork.shuntCount) + (mEventSwitch ? 6 : 0);
    if (mValues.size() != expectedValues)
      throw std::runtime_error("Physical-output count changed");
    publish();
    checkFinite();
    if (mOptions.validate()) {
      mTrajectory.open(mOptions.directory / "trajectory.csv");
      mTrajectory << std::setprecision(17) << "step,time";
      for (const auto &block : mBlocks)
        for (size_t i = 0; i < block.size; ++i)
          mTrajectory << ',' << block.name << '[' << i << ']';
      mTrajectory << '\n';
      trajectoryRow(0, 0.); // Initial conditions, distinct from solved step 1.
      if (mEventSwitch) {
        mWaveforms.open(mOptions.directory / "waveforms.csv");
        mWaveforms
            << std::setprecision(17)
            << "time,voltage_a_V,voltage_b_V,voltage_c_V,"
               "branch_current_a_A,branch_current_b_A,branch_current_c_A,"
               "event_current_a_A,event_power_W,generator_speed_pu\n";
      }
    }
  }

  void beginLoop() { mLoopStart = BenchmarkClock::now(); }
  void endLoop() { mLoopEnd = BenchmarkClock::now(); }

  // Same cached physical-output publication for timing and validation runs.
  void observe(size_t step, Real time) {
    const auto start =
        mOptions.profile() ? BenchmarkClock::now() : BenchmarkClock::time_point{};
    publish();
    if (mOptions.profile())
      mOutputSeconds += elapsed(start);
    if (mOptions.validate()) {
      checkFinite();
      trajectoryRow(step, time);
      if (mEventSwitch) {
        const Matrix &v = mEventSwitch->mIntfVoltage->get();
        const Matrix &i = mEventSwitch->mIntfCurrent->get();
        auto generator =
            mNetwork.system.component<EMT::Ph3::SynchronGenerator4OrderVBR>(
                "Gen" + std::to_string(eventBus));
        const Matrix &voltage =
            mNetwork.nodes.at(eventBus - 1)->mVoltage->get();
        const Matrix &branch = mBranchCurrent->get();
        mWaveforms << time;
        for (UInt phase = 0; phase < 3; ++phase)
          mWaveforms << ',' << voltage(phase);
        for (UInt phase = 0; phase < 3; ++phase)
          mWaveforms << ',' << branch(phase);
        mWaveforms << ',' << i(0) << ',' << v.cwiseProduct(i).sum() << ','
                   << generator->attributeTyped<Real>("w_r")->get() << '\n';
        const bool closed = mEventSwitch->mPoleClosedA->get();
        if (mFirstClosedTime < 0. && closed)
          mFirstClosedTime = time;
        if (mFirstClosedTime >= 0. && mFirstOpenTime < 0. && !closed)
          mFirstOpenTime = time;
      }
    }
  }

  void write(Simulation &sim, nlohmann::json extra = {}) {
    checkFinite();
    if (mPublications != mOptions.steps + 1 ||
        sim.timeStepCount() != static_cast<Int>(mOptions.steps))
      throw std::runtime_error("Step/output count mismatch");
    size_t componentUpdates = 0;
    for (const auto &component : mNetwork.passive) {
      const auto count =
          component->attributeTyped<UInt>("output_updates")->get();
      if (count != mOptions.steps)
        throw std::runtime_error("Missing output update: " + component->name());
      componentUpdates += count;
    }
    if (mNetwork.group &&
        (mNetwork.group->voltageUpdateCount() != mOptions.steps ||
         mNetwork.group->outputUpdateCount() != mOptions.steps))
      throw std::runtime_error("Missing grouped voltage/current publication");
    const size_t coordinates = allocatedCoordinates(mNetwork);
    size_t primitives = 0;
    for (const auto &component : mNetwork.passive) {
      if (auto line = std::dynamic_pointer_cast<EMT::Ph3::PiLine>(component))
        primitives += line->primitiveCount();
      else if (auto trafo =
                   std::dynamic_pointer_cast<EMT::Ph3::Transformer>(component))
        primitives += trafo->primitiveCount();
      else if (auto load =
                   std::dynamic_pointer_cast<EMT::Ph3::RXLoad>(component))
        primitives += load->primitiveCount();
      else if (auto shunt =
                   std::dynamic_pointer_cast<EMT::Ph3::Shunt>(component))
        primitives += shunt->primitiveCount();
    }
    const size_t ssnCount = mOptions.formulation == "native" ? 0
                            : (mOptions.formulation == "elementary" ||
                               mOptions.formulation == "elementary_wrapped")
                                ? primitives
                            : mNetwork.group ? 1
                                             : mNetwork.passive.size();
    sim.logLUTimes();
    if (mOptions.profile())
      mScheduler->writeCSV(
          (mOptions.directory / "scheduler_phases.csv").string());
    nlohmann::json result = {
        {"case", mCase},
        {"grid", gridBusCount},
        {"event_bus", eventBus},
        {"voltage_bus", eventBus},
        {"monitored_branch", mNetwork.monitoredBranch},
        {"monitored_branch_terminal", mNetwork.monitoredTerminal},
        {"monitored_branch_from_bus", mNetwork.monitoredFrom},
        {"monitored_branch_to_bus", mNetwork.monitoredTo},
        {"voltage_base_V", baseV(busFor(eventBus))},
        {"current_base_A",
         baseMVA * 1e6 / (std::sqrt(3.) * baseV(busFor(eventBus)))},
        {"frequency_Hz", frequency},
        {"total_load_power_W", totalLoadPower()},
        {"load_step_fraction", loadStepPower() / totalLoadPower()},
        {"generator_model", "SynchronGenerator4OrderVBR"},
        {"generator_controls", "none"},
        {"generator_dynamic_parameters",
         {{"H_s", H},
          {"Xd_pu", Xd},
          {"Xq_pu", Xq},
          {"Xl_pu", Xl},
          {"Xd_transient_pu", XdTransient},
          {"Xq_transient_pu", XqTransient},
          {"Td0_transient_s", Td0Transient},
          {"Tq0_transient_s", Tq0Transient}}},
        {"formulation", mOptions.formulation},
        {"elementary_layout", mOptions.formulation == "elementary" ? "flat"
                              : mOptions.formulation == "elementary_wrapped"
                                  ? "wrapped"
                                  : "not_applicable"},
        {"runtime_composite_count",
         std::count_if(
             mNetwork.system.mComponents.begin(),
             mNetwork.system.mComponents.end(),
             [](const auto &component) {
               return static_cast<bool>(
                   std::dynamic_pointer_cast<CPS::CompositePowerComp<Real>>(
                       component));
             })},
        {"solver", mSolver},
        {"mode", mOptions.mode},
        {"steps", mOptions.steps},
        {"timestep_s", timeStep},
        {"horizon_s", mOptions.steps * timeStep},
        {"bus_count", mNetwork.nodes.size()},
        {"generator_count", mNetwork.machines.size()},
        {"injection_count", mNetwork.injections.size()},
        {"model_adaptation", modelAdaptation},
        {"passive_count", mNetwork.passive.size()},
        {"rlc_primitive_count", primitives},
        {"ssn_count", ssnCount},
        {"load_count", mNetwork.loadCount},
        {"shunt_count", mNetwork.shuntCount},
        {"line_count", mNetwork.lineCount},
        {"transformer_count", mNetwork.transformerCount},
        {"active_companion_coordinates", mNetwork.activeCoordinates},
        {"allocated_companion_coordinates", coordinates},
        {"group_states", mNetwork.group ? mNetwork.group->stateCount() : 0},
        {"retained_bus_count",
         mNetwork.group ? mNetwork.retained.size() : mNetwork.nodes.size()},
        {"internal_bus_count",
         mNetwork.group ? mNetwork.nodes.size() - mNetwork.retained.size() : 0},
        {"group_map_mode", mNetwork.group ? "precomputed_maps" : "none"},
        {"output_values", mValues.size()},
        {"output_publications", mPublications},
        {"component_output_updates", componentUpdates},
        {"task_count", sim.taskCount()},
        {"executed_task_count", mScheduler->executedTaskCount()},
        {"scheduler_blocks", mScheduler->blockCount()},
        {"powerflow_preparation_s", mNetwork.powerflowSeconds},
        {"simulation_initialize_s", elapsed(mInitializeStart, mInitializeEnd)},
        {"preparation_total_s", elapsed(mPreparationStart, mLoopStart)},
        {"loop_s", elapsed(mLoopStart, mLoopEnd)},
        {"output_vector_publication_s", mOutputSeconds},
        {"group_constant_map_preparation_s",
         mNetwork.group ? mNetwork.group->realizationPreparationSeconds() : 0.},
        {"partial_refactorization", false},
        {"dense_static_preconditioner", false},
        {"threads", 1},
        {"first_closed_sample_time_s", mFirstClosedTime},
        {"first_open_sample_time_s", mFirstOpenTime}};
    result.update(extra);
    std::ofstream(mOptions.directory / "result.json") << result.dump(2) << '\n';
    SPDLOG_LOGGER_INFO(Logger::get("ssn-benchmark"),
                       "{} {}: loop {:.6f} s; results in {}", mCase,
                       mOptions.formulation, result["loop_s"].get<double>(),
                       mOptions.directory.string());
  }
};
} // namespace IEEEBenchmark
