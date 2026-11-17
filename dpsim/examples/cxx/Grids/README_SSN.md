# IEEE network SSN examples

These C++ examples simulate IEEE9, IEEE39, IEEE118 and IEEE300 in the EMT
domain with a 50 µs timestep. Network data are embedded in the headers;
no external grid files or Python interface are required.

## Build and run IEEE300

Use a C++17 compiler, CMake 3.14 or newer, Ninja, Git, and the Eigen3,
spdlog (1.10 or newer) and fmt (8.1.1 or newer) development packages.
CMake fetches SuiteSparse/KLU, nlohmann/json and readerwriterqueue on the
first configuration, so that step normally requires internet access.

From the repository root:

```bash
scripts/build_ssn_examples.sh IEEE300_SSN_Sparse
build_ssn_examples/dpsim/examples/cxx/IEEE300_SSN_Sparse \
  network timing 400 /tmp/ieee300-network
```

This runs IEEE300 with one SSN for the linear network and KLU for the global
MNA solution. It simulates 20 ms without a switching event. Results are saved
under `/tmp/ieee300-network`.

The helper configures a Release build with KLU enabled, OpenMP disabled,
native CPU optimization enabled and LTO disabled. With no target arguments,
it builds all twelve examples. Set `DPSIM_SSN_BUILD` to choose another build
directory or `BUILD_JOBS` to change the default four compilation jobs.
An existing CMake dependency cache can be reused by setting
`DPSIM_DEPENDENCY_CACHE` to its `_deps` directory.

## Targets and arguments

Replace `N` with `9`, `39`, `118` or `300`:

| Target | Solver and operation | Defaults |
| --- | --- | --- |
| `IEEEN_SSN_Sparse` | KLU, no event | `network timing` |
| `IEEEN_SSN_Dense` | DenseLU, no event | `network timing` |
| `IEEEN_SSN_Accuracy` | KLU, load step or cleared fault | `network validate 4000`, load event |

All arguments are positional and optional:

```text
executable [representation] [mode] [steps] [output-directory]
```

| Argument | Values |
| --- | --- |
| Representation | `network`: one linear-network SSN; `equipment`: one SSN per line/transformer/load; `elementary`: individual R/L/C SSNs; `native`: companion-circuit components |
| Mode | `timing`: overall runtime; `profile`: operation timings; `validate`: save electrical trajectories |
| Steps | Positive integer; sparse/dense defaults are 20,000 / 4,000 / 2,000 / 400 for IEEE9 / 39 / 118 / 300 |
| Output directory | Defaults to `logs/<target>/<representation>_<mode>`, relative to the process working directory |

`elementary_wrapped` additionally retains composite equipment wrappers for
comparison with the flat `elementary` construction.

Accuracy targets accept a final `load` or `fault` argument and require
`validate` with more than 2,000 steps. The event starts at 50 ms; a fault
clears at 100 ms. For example:

```bash
scripts/build_ssn_examples.sh IEEE300_SSN_Accuracy
build_ssn_examples/dpsim/examples/cxx/IEEE300_SSN_Accuracy \
  network validate 4000 /tmp/ieee300-fault fault
```

## Outputs and profiling

Every run writes `result.json` with the configuration and overall timings,
plus solver diagnostics under `logs/`. Profile mode adds
`scheduler_phases.csv`. Validation mode writes `trajectory.csv`; accuracy
targets also write the smaller `waveforms.csv` containing representative
bus voltages and branch currents. The executables do not open plots.

For an external profiler, use `timing` mode: detailed solver/scheduler
profiling and trajectory writing are disabled. The sequential scheduler
still calculates the physical outputs and copies them to memory each step.
Keep the output dependencies in the benchmark helper; otherwise the
scheduler can omit tasks whose outputs are not requested.

Profile initialization separately from timestep execution. Network SSN maps
are prepared during initialization and reused in the loop. Increase the step
count if the profiler needs more samples. For source attribution, add debug
symbols to the Release build, for example by setting
`CMAKE_CXX_FLAGS_RELEASE="-O3 -DNDEBUG -g"` when configuring. The project also
applies its optimized-build flags. The helper's compilation job count does
not change the single-threaded simulation schedule.

## CMake Tools in VS Code

Select **Release**, enable `WITH_SPARSE` and `WITH_KLU`, and disable
`WITH_OPENMP`. Configure the project, select `IEEE300_SSN_Sparse` as both
build and launch target, then use **CMake: Run Without Debugging**.
Use the optional-dependency settings in the build helper when configuring
a minimal C++ environment.

To pass arguments, merge this into the workspace's `.vscode/settings.json`:

```json
{
  "cmake.debugConfig": {
    "cwd": "${workspaceFolder}",
    "args": ["network", "timing", "400", "/tmp/ieee300-network"]
  }
}
```

Without explicit arguments, the sparse target runs `network timing 400`.
Its default output directory follows the working directory selected by
CMake Tools, which may be the executable directory.
