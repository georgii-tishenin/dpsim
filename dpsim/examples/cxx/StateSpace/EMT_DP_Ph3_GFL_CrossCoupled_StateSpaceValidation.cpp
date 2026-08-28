// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

// Reuse the common benchmark and result writers while providing a distinct
// executable whose Studies 1--3 all enable nominal-frequency current
// cross-coupling.  The original validation executable remains unchanged.
#define DP_SIM_GFL_CROSS_COUPLED_VALIDATION
#include "EMT_DP_Ph3_GFL_StateSpaceValidation.cpp"

int main(int argc, char **argv) {
  Bool runTimeDomain = false;
  Bool runBenchmarkScan = false;
  Bool runWeakGridStudy = false;
  Bool runCandidateCheck = false;
  Bool runEmtVariableDiagnostics = false;
  Bool runAll = false;
  for (Int idx = 1; idx < argc; ++idx) {
    if (String(argv[idx]) == "--time-domain" || String(argv[idx]) == "--all")
      runTimeDomain = true;
    if (String(argv[idx]) == "--all")
      runAll = true;
    if (String(argv[idx]) == "--benchmark-scan")
      runBenchmarkScan = true;
    if (String(argv[idx]) == "--weak-grid-study")
      runWeakGridStudy = true;
    if (String(argv[idx]) == "--candidate-check")
      runCandidateCheck = true;
    if (String(argv[idx]) == "--emt-variable-diagnostics")
      runEmtVariableDiagnostics = true;
  }
  EMTDPPh3GFLStateSpaceValidation example(
      runTimeDomain, true,
      "EMT_DP_Ph3_GFL_CrossCoupled_StateSpaceValidation");
  if (runBenchmarkScan)
    example.runBenchmarkSelectionScan();
  else if (runWeakGridStudy)
    example.runWeakGridStudy();
  else if (runCandidateCheck)
    example.runCandidateCheck();
  else if (runEmtVariableDiagnostics)
    example.runEmtVariableDiagnostics();
  else {
    example.run();
    if (runAll)
      example.runWeakGridStudy();
  }
  return 0;
}
