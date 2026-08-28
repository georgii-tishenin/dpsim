// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include "../Examples.h"
#include "../GeneratorFactory.h"

#include <DPsim.h>

using namespace DPsim;
using namespace CPS;

namespace {
// Controller parameters and their equation-based derivation live in Examples.h.
using GfmParams = CPS::CIM::Examples::Components::GFM::Ieee9SsnGridForming;
using GflParams = CPS::CIM::Examples::Components::GFL::Ieee9AvVsi;
} // namespace

SystemTopology buildTopology(CommandLineArgs &args,
                             std::shared_ptr<DataLoggerInterface> logger) {

  String simName = args.name;

  CPS::CIM::Examples::Grids::IEEE9::ScenarioConfig ieee9(args.sysFreq);

  const auto optionBool = [&](const String &key, Bool defaultValue) {
    return args.options.find(key) != args.options.end()
               ? args.getOptionBool(key)
               : defaultValue;
  };
  const auto optionReal = [&](const String &key, Real defaultValue) {
    return args.options.find(key) != args.options.end()
               ? args.getOptionReal(key)
               : defaultValue;
  };
  const Bool useNetworkEquivalent =
      optionBool("network_equivalent", false);
  const Real networkEquivalentImpedanceScale =
      optionReal("network_equivalent_impedance_scale", 1.0);
  if (!(networkEquivalentImpedanceScale > 0.0))
    throw std::invalid_argument(
        "network_equivalent_impedance_scale must be positive.");
  if (!useNetworkEquivalent && networkEquivalentImpedanceScale != 1.0)
    throw std::invalid_argument(
        "network_equivalent_impedance_scale requires network_equivalent.");
  const Bool useMagnetizingBranches =
      optionBool("transformer_magnetizing", false);
  const String gflFormulation =
      args.options.find("gfl_formulation") != args.options.end()
          ? args.getOptionString("gfl_formulation")
          : "legacy_variable";
  const Bool enableGflCurrentCrossCoupling =
      optionBool("gfl_current_cross_coupling", false);
  const Bool enableDefaultLineConductance =
      optionBool("line_default_conductance", true);
  const Bool enableLine96Breaker = optionBool("line96_breaker", false);
  const Bool compactLog = optionBool("compact_log", false);
  // Representative no-load values, configurable for sensitivity studies.
  const Real transformerCoreLossPu =
      optionReal("transformer_p0_pu", 1.0e-3);
  const Real transformerMagnetizingQPu =
      optionReal("transformer_q0_pu", 1.0e-2);
  const Real steadyLoadStepPowerPu =
      optionReal("steady_load_step_p_pu", 0.0);
  if (steadyLoadStepPowerPu < 0.0)
    throw std::invalid_argument(
        "steady_load_step_p_pu must be nonnegative.");
  const Real load8ActivePower =
      ieee9.load8.RealPower + steadyLoadStepPowerPu * 100.0e6;
  const String outagedLine =
      args.options.find("outaged_line") != args.options.end()
          ? args.getOptionString("outaged_line")
          : "none";
  const std::vector<String> lineNames{
      ieee9.line54.Name, ieee9.line64.Name, ieee9.line75.Name,
      ieee9.line96.Name, ieee9.line78.Name, ieee9.line89.Name};
  if (outagedLine != "none" &&
      std::find(lineNames.begin(), lineNames.end(), outagedLine) ==
          lineNames.end()) {
    throw std::invalid_argument("Unknown outaged_line: " + outagedLine);
  }
  const auto lineInService = [&](const auto &line) {
    return line->name() != outagedLine;
  };

  // POWER FLOW FOR INITIALIZATION
  CPS::Logger::get(args.name)->info("Creating power flow initialization.");

  String simNamePF = simName + "_PF";
  CPS::Logger::setLogDir("logs/" + simNamePF);

  // Nodes
  auto n1PF = SimNode<Complex>::make("BUS1", PhaseType::Single);
  auto n2PF = SimNode<Complex>::make("BUS2", PhaseType::Single);
  auto n3PF = SimNode<Complex>::make("BUS3", PhaseType::Single);
  auto n4PF = SimNode<Complex>::make("BUS4", PhaseType::Single);
  auto n5PF = SimNode<Complex>::make("BUS5", PhaseType::Single);
  auto n6PF = SimNode<Complex>::make("BUS6", PhaseType::Single);
  auto n7PF = SimNode<Complex>::make("BUS7", PhaseType::Single);
  auto n8PF = SimNode<Complex>::make("BUS8", PhaseType::Single);
  auto n9PF = SimNode<Complex>::make("BUS9", PhaseType::Single);

  auto gen1PF = SP::Ph1::SynchronGenerator::make(ieee9.gen1.Name,
                                                 CPS::Logger::Level::off);
  gen1PF->setParameters(ieee9.gen1.RatedPower, ieee9.gen1.RatedVoltage,
                        ieee9.gen1.InitialPower, ieee9.gen1.InitialVoltage,
                        ieee9.gen1.BusType);
  gen1PF->setBaseVoltage(ieee9.gen1.RatedVoltage);

  const GfmParams gfmSeed;
  const auto [gfm2PPcc, unusedGfm2QPcc] =
      Math::pccPowerFromFilterPowerReference(
          ieee9.gen2.InitialPower, ieee9.gen2.InitialPowerReactive,
          gfmSeed.Rc, ieee9.gen2.InitialVoltage);
  (void)unusedGfm2QPcc;
  auto gen2PF = SP::Ph1::SynchronGenerator::make(ieee9.gen2.Name,
                                                 CPS::Logger::Level::off);
  gen2PF->setParameters(ieee9.gen2.RatedPower, ieee9.gen2.RatedVoltage,
                        gfm2PPcc, ieee9.gen2.InitialVoltage,
                        ieee9.gen2.BusType);
  gen2PF->setBaseVoltage(ieee9.gen2.RatedVoltage);

  // gen3's PF image: a negative PQ load injecting the PCC-side power (the
  // rc-corrected filter reference).
  const GflParams gfl;
  const auto [gfl3PPcc, gfl3QPcc] = Math::pccPowerFromFilterPowerReference(
      ieee9.gen3.InitialPower, ieee9.gen3.InitialPowerReactive, gfl.Rc,
      ieee9.gen3.RatedVoltage);
  auto gfl3PF = SP::Ph1::Load::make(ieee9.gen3.Name, CPS::Logger::Level::off);
  gfl3PF->setParameters(-gfl3PPcc, -gfl3QPcc, ieee9.gen3.RatedVoltage);
  gfl3PF->modifyPowerFlowBusType(PowerflowBusType::PQ);

  // Loads
  auto load5PF = SP::Ph1::Load::make(ieee9.load5.Name, CPS::Logger::Level::off);
  load5PF->setParameters(ieee9.load5.RealPower, ieee9.load5.ReactivePower,
                         ieee9.load5.BaseVoltage);
  load5PF->modifyPowerFlowBusType(PowerflowBusType::PQ);

  auto load6PF = SP::Ph1::Load::make(ieee9.load6.Name, CPS::Logger::Level::off);
  load6PF->setParameters(ieee9.load6.RealPower, ieee9.load6.ReactivePower,
                         ieee9.load6.BaseVoltage);
  load6PF->modifyPowerFlowBusType(PowerflowBusType::PQ);

  auto load8PF = SP::Ph1::Load::make(ieee9.load8.Name, CPS::Logger::Level::off);
  load8PF->setParameters(load8ActivePower, ieee9.load8.ReactivePower,
                          ieee9.load8.BaseVoltage);
  load8PF->modifyPowerFlowBusType(PowerflowBusType::PQ);

  // When the EMT transformers use physical magnetizing branches, include the
  // same no-load powers in the power flow so the dynamic initialization uses
  // the same operating point.
  std::vector<std::shared_ptr<SP::Ph1::Load>> magnetizingLoadsPF;
  if (useMagnetizingBranches) {
    const auto makeMagnetizingLoad = [&](const String &name, Real ratedPower,
                                         Real baseVoltage,
                                         const SimNode<Complex>::Ptr &node) {
      auto load = SP::Ph1::Load::make(name, CPS::Logger::Level::off);
      load->setParameters(transformerCoreLossPu * ratedPower,
                          transformerMagnetizingQPu * ratedPower,
                          baseVoltage);
      load->modifyPowerFlowBusType(PowerflowBusType::PQ);
      load->connect({node});
      magnetizingLoadsPF.push_back(load);
    };
    makeMagnetizingLoad("TR14_MAG_PF", ieee9.transf14.RatedPower,
                        ieee9.transf14.VoltageHVSide, n4PF);
    makeMagnetizingLoad("TR27_MAG_PF", ieee9.transf27.RatedPower,
                        ieee9.transf27.VoltageHVSide, n7PF);
    makeMagnetizingLoad("TR39_MAG_PF", ieee9.transf39.RatedPower,
                        ieee9.transf39.VoltageHVSide, n9PF);
  }

  // Transmission Lines

  auto line54PF =
      SP::Ph1::PiLine::make(ieee9.line54.Name, CPS::Logger::Level::off);
  line54PF->setParameters(ieee9.line54.Resistance, ieee9.line54.Inductance,
                          ieee9.line54.Capacitance, ieee9.line54.Conductance);
  line54PF->setBaseVoltage(ieee9.line54.BaseVoltage);

  auto line64PF =
      SP::Ph1::PiLine::make(ieee9.line64.Name, CPS::Logger::Level::off);
  line64PF->setParameters(ieee9.line64.Resistance, ieee9.line64.Inductance,
                          ieee9.line64.Capacitance, ieee9.line64.Conductance);
  line64PF->setBaseVoltage(ieee9.line64.BaseVoltage);

  auto line75PF =
      SP::Ph1::PiLine::make(ieee9.line75.Name, CPS::Logger::Level::off);
  line75PF->setParameters(ieee9.line75.Resistance, ieee9.line75.Inductance,
                          ieee9.line75.Capacitance, ieee9.line75.Conductance);
  line75PF->setBaseVoltage(ieee9.line75.BaseVoltage);

  auto line96PF =
      SP::Ph1::PiLine::make(ieee9.line96.Name, CPS::Logger::Level::off);
  line96PF->setParameters(ieee9.line96.Resistance, ieee9.line96.Inductance,
                          ieee9.line96.Capacitance, ieee9.line96.Conductance);
  line96PF->setBaseVoltage(ieee9.line96.BaseVoltage);

  auto line78PF =
      SP::Ph1::PiLine::make(ieee9.line78.Name, CPS::Logger::Level::off);
  line78PF->setParameters(ieee9.line78.Resistance, ieee9.line78.Inductance,
                          ieee9.line78.Capacitance, ieee9.line78.Conductance);
  line78PF->setBaseVoltage(ieee9.line78.BaseVoltage);

  auto line89PF =
      SP::Ph1::PiLine::make(ieee9.line89.Name, CPS::Logger::Level::off);
  line89PF->setParameters(ieee9.line89.Resistance, ieee9.line89.Inductance,
                          ieee9.line89.Capacitance, ieee9.line89.Conductance);
  line89PF->setBaseVoltage(ieee9.line89.BaseVoltage);

  // Transformers

  auto transf14PF =
      SP::Ph1::Transformer::make(ieee9.transf14.Name, CPS::Logger::Level::off);
  transf14PF->setParameters(
      ieee9.transf14.VoltageLVSide, ieee9.transf14.VoltageHVSide,
      ieee9.transf14.Ratio, 0.0, // No phase shift (ratioPhase = 0.0)
      networkEquivalentImpedanceScale * ieee9.transf14.Resistance,
      networkEquivalentImpedanceScale * ieee9.transf14.Inductance);
  transf14PF->setBaseVoltage(ieee9.transf14.VoltageHVSide);

  auto transf27PF =
      SP::Ph1::Transformer::make(ieee9.transf27.Name, CPS::Logger::Level::off);
  transf27PF->setParameters(ieee9.transf27.VoltageLVSide,
                            ieee9.transf27.VoltageHVSide, ieee9.transf27.Ratio,
                            0.0, ieee9.transf27.Resistance,
                            ieee9.transf27.Inductance);
  transf27PF->setBaseVoltage(ieee9.transf27.VoltageHVSide);

  auto transf39PF =
      SP::Ph1::Transformer::make(ieee9.transf39.Name, CPS::Logger::Level::off);
  transf39PF->setParameters(ieee9.transf39.VoltageLVSide,
                            ieee9.transf39.VoltageHVSide, ieee9.transf39.Ratio,
                            0.0, ieee9.transf39.Resistance,
                            ieee9.transf39.Inductance);
  transf39PF->setBaseVoltage(ieee9.transf39.VoltageHVSide);

  // Connect components
  gen1PF->connect({n1PF});
  gen2PF->connect({n2PF});
  gfl3PF->connect({n3PF});

  load5PF->connect({n5PF});
  load6PF->connect({n6PF});
  load8PF->connect({n8PF});

  if (lineInService(line54PF))
    line54PF->connect({n5PF, n4PF});
  if (lineInService(line64PF))
    line64PF->connect({n6PF, n4PF});
  if (lineInService(line75PF))
    line75PF->connect({n7PF, n5PF});
  if (lineInService(line96PF))
    line96PF->connect({n9PF, n6PF});
  if (lineInService(line78PF))
    line78PF->connect({n7PF, n8PF});
  if (lineInService(line89PF))
    line89PF->connect({n8PF, n9PF});

  transf14PF->connect({n1PF, n4PF});
  transf27PF->connect({n2PF, n7PF});
  transf39PF->connect({n3PF, n9PF});

  // Create system topology
  SystemComponentList componentsPF{gen1PF,     gen2PF,     gfl3PF,
                                   load5PF,    load6PF,    load8PF,
                                   transf14PF, transf27PF, transf39PF};
  for (const auto &line :
       {line54PF, line64PF, line75PF, line96PF, line78PF, line89PF}) {
    if (lineInService(line))
      componentsPF.push_back(line);
  }
  componentsPF.insert(componentsPF.end(), magnetizingLoadsPF.begin(),
                      magnetizingLoadsPF.end());
  auto systemPF = SystemTopology(
      ieee9.nomFreq,
      SystemNodeList{n1PF, n2PF, n3PF, n4PF, n5PF, n6PF, n7PF, n8PF, n9PF},
      componentsPF);

  // Logger
  auto loggerPF = DataLogger::make(simNamePF, CPS::Logger::Level::off);
  // Log node voltages
  loggerPF->logAttribute("v_bus1", n1PF->attribute("v"));
  loggerPF->logAttribute("v_bus2", n2PF->attribute("v"));
  loggerPF->logAttribute("v_bus3", n3PF->attribute("v"));
  loggerPF->logAttribute("v_bus4", n4PF->attribute("v"));
  loggerPF->logAttribute("v_bus5", n5PF->attribute("v"));
  loggerPF->logAttribute("v_bus6", n6PF->attribute("v"));
  loggerPF->logAttribute("v_bus7", n7PF->attribute("v"));
  loggerPF->logAttribute("v_bus8", n8PF->attribute("v"));
  loggerPF->logAttribute("v_bus9", n9PF->attribute("v"));
  // Log node powers
  loggerPF->logAttribute("s_bus1", n1PF->attribute("s"));
  loggerPF->logAttribute("s_bus2", n2PF->attribute("s"));
  loggerPF->logAttribute("s_bus3", n3PF->attribute("s"));
  loggerPF->logAttribute("s_bus4", n4PF->attribute("s"));
  loggerPF->logAttribute("s_bus5", n5PF->attribute("s"));
  loggerPF->logAttribute("s_bus6", n6PF->attribute("s"));
  loggerPF->logAttribute("s_bus7", n7PF->attribute("s"));
  loggerPF->logAttribute("s_bus8", n8PF->attribute("s"));
  loggerPF->logAttribute("s_bus9", n9PF->attribute("s"));

  // Run power flow simulation
  Simulation simPF(simNamePF, CPS::Logger::Level::off);
  simPF.setSystem(systemPF);
  simPF.setTimeStep(args.timeStep);
  simPF.setFinalTime(1 * args.timeStep);
  simPF.setDomain(Domain::SP);
  simPF.setSolverType(Solver::Type::NRP);
  simPF.setSolverAndComponentBehaviour(Solver::Behaviour::Simulation);
  simPF.addLogger(loggerPF);
  simPF.run();
  CPS::Logger::get(args.name)->info("Power flow simulation finished.");

  // Seed voltages for the step-2 inverters: converged n2 (GFM PCC) and
  // n3 (GFL PCC) node voltages, magnitude [kV] and angle [deg]. Dedicated
  // console-enabled logger so the seeds survive the setLogDir churn.
  Complex v2PF = n2PF->singleVoltage();
  Complex v3PF = n3PF->singleVoltage();
  const Complex gen2PowerPF = gen2PF->getApparentPower();
  auto seedLog = CPS::Logger::get(simNamePF + "_seed", CPS::Logger::Level::info,
                                  CPS::Logger::Level::info);
  seedLog->info(
      "PF seed n2 (BUS2): |V| = {:.4f} kV ({:.4f} pu), angle = {:.4f} deg",
      std::abs(v2PF) / 1e3, std::abs(v2PF) / ieee9.gen2.RatedVoltage,
      std::arg(v2PF) * 180.0 / PI);
  seedLog->info(
      "PF seed n3 (BUS3): |V| = {:.4f} kV ({:.4f} pu), angle = {:.4f} deg",
      std::abs(v3PF) / 1e3, std::abs(v3PF) / ieee9.gen3.RatedVoltage,
      std::arg(v3PF) * 180.0 / PI);
  seedLog->info("PF seed GEN2: P = {:.4f} MW, Q = {:.4f} Mvar",
                gen2PowerPF.real() / 1e6, gen2PowerPF.imag() / 1e6);

  // DYNAMIC SIMULATION - EMT
  CPS::Logger::get(args.name)->info("Dynamic simulation initialization.");
  String simNameEMT = simName + "_EMT";
  CPS::Logger::setLogDir("logs/" + simNameEMT);

  // Nodes
  auto n1EMT = SimNode<Real>::make("BUS1", PhaseType::ABC);
  auto n2EMT = SimNode<Real>::make("BUS2", PhaseType::ABC);
  auto n3EMT = SimNode<Real>::make("BUS3", PhaseType::ABC);
  auto n4EMT = SimNode<Real>::make("BUS4", PhaseType::ABC);
  auto n5EMT = SimNode<Real>::make("BUS5", PhaseType::ABC);
  auto n6EMT = SimNode<Real>::make("BUS6", PhaseType::ABC);
  auto n7EMT = SimNode<Real>::make("BUS7", PhaseType::ABC);
  auto n8EMT = SimNode<Real>::make("BUS8", PhaseType::ABC);
  auto n9EMT = SimNode<Real>::make("BUS9", PhaseType::ABC);
  auto n96BreakerEMT =
      SimNode<Real>::make("LINE96_BREAKER_NODE", PhaseType::ABC);

  std::shared_ptr<SimPowerComp<Real>> gen1EMT;
  if (useNetworkEquivalent) {
    auto source = EMT::Ph3::NetworkInjection::make(
        ieee9.gen1.Name, CPS::Logger::Level::off);
    source->setParameters(
        Math::singlePhaseVariableToThreePhase(n1PF->singleVoltage()),
        ieee9.nomFreq);
    gen1EMT = source;
  } else {
    auto generator = EMT::Ph3::SynchronGenerator4OrderVBR::make(
        ieee9.gen1.Name, CPS::Logger::Level::off);

    generator->setOperationalParametersPerUnit(
        ieee9.gen1.RatedPower, ieee9.gen1.RatedVoltage, ieee9.nomFreq,
        ieee9.gen1.H, ieee9.gen1.Xd, ieee9.gen1.Xq, ieee9.gen1.Xa,
        ieee9.gen1.XdPrime, ieee9.gen1.XqPrime, ieee9.gen1.TdoPrime,
        ieee9.gen1.TqoPrime);

    auto exciter1Params = std::make_shared<Signal::ExciterDC1SimpParameters>();
    exciter1Params->Ta = ieee9.exc1.TA;
    exciter1Params->Ka = ieee9.exc1.KA;
    exciter1Params->Tef = ieee9.exc1.TE;
    exciter1Params->Kef = ieee9.exc1.KE;
    exciter1Params->Tf = ieee9.exc1.TF;
    exciter1Params->Kf = ieee9.exc1.KF;
    exciter1Params->Tr = 0.01;
    exciter1Params->MaxVa = ieee9.exc1.VRmax;
    exciter1Params->MinVa = ieee9.exc1.VRmin;
    exciter1Params->Bef = std::log(ieee9.exc1.S_EX2 / ieee9.exc1.S_EX1) /
                          (ieee9.exc1.EX2 - ieee9.exc1.EX1);
    exciter1Params->Aef =
        ieee9.exc1.S_EX1 / std::exp(exciter1Params->Bef * ieee9.exc1.EX1);
    auto exciter1 =
        Signal::ExciterDC1Simp::make("Gen1_Exciter", CPS::Logger::Level::off);
    exciter1->setParameters(exciter1Params);
    generator->addExciter(exciter1);

    auto turbineGovernor1 = Signal::TurbineGovernorType1::make(
        "Gen1_TurbineGovernor", CPS::Logger::Level::off);
    turbineGovernor1->setParameters(
        ieee9.gov1.T2, 1.0, 1.0, ieee9.gov1.T3, ieee9.gov1.T1, ieee9.gov1.R,
        ieee9.gov1.Vmin, ieee9.gov1.Vmax, 1.0);
    generator->addGovernor(turbineGovernor1);
    gen1EMT = generator;
  }

  const Real omegaN = 2.0 * PI * ieee9.nomFreq;

  // gen2 replaced by a grid-forming SSN inverter, keeping the GEN2 identity.
  // nominalVoltage is the peak phase target at the 1.025 pu PV setpoint.
  const Real gfmNominalVoltage = RMS3PH_TO_PEAK1PH * ieee9.gen2.InitialVoltage;
  GfmParams gfm = gfmSeed;
  // Optional overrides for studying the grid-forming tuning from a notebook.
  auto opt = [&](const String &key, Real def) {
    return args.options.find(key) != args.options.end()
               ? args.getOptionReal(key)
               : def;
  };
  gfm.dampingCoefficient = opt("gfm_d", gfm.dampingCoefficient);
  gfm.KpVoltage = opt("gfm_kpv", gfm.KpVoltage);
  gfm.KiVoltage = opt("gfm_kiv", gfm.KiVoltage);
  gfm.reactivePowerDroop = opt("gfm_dq", gfm.reactivePowerDroop);
  gfm.reactiveDroopCutoff = opt("gfm_dqc", gfm.reactiveDroopCutoff);
  const Real gfmFeedforward = opt("gfm_ff", gfm.gridCurrentFeedforward);
  const Real gfmVirtualResistance = opt("gfm_rv", 0.0);

  auto gen2EMT = EMT::Ph3::SSN_GFM::make(ieee9.gen2.Name, ieee9.gen2.Name,
                                         CPS::Logger::Level::off);
  gen2EMT->setNumericalLinearizationParameters(1e-6, 1e-8);
  gen2EMT->setParameters(gfm.Lf, gfm.Cf, gfm.Rf, gfm.Rc, gfmNominalVoltage,
                         omegaN, ieee9.gen2.InitialPower, gen2PowerPF.imag(),
                         gfm.virtualInertia,
                         gfm.dampingCoefficient, gfm.voltageDroopGain,
                         gfm.reactiveIntegralGain, gfm.KpVoltage, gfm.KiVoltage,
                         gfm.KpCurrent, gfm.KiCurrent, gfm.activeDampingGain,
                         gfm.powerFilterCutoff, gfm.delayBandwidth);
  // Grid-connected control: no grid-current feedforward, proportional Q-V droop.
  gen2EMT->setGridCurrentFeedforward(gfmFeedforward);
  gen2EMT->setVirtualImpedance(gfmVirtualResistance, 0.0);
  gen2EMT->setReactivePowerDroop(gfm.reactivePowerDroop,
                                 gfm.reactiveDroopCutoff);

  // gen3 replaced by a grid-following averaged VSI (SSN), keeping the GEN3
  // identity so the topology wiring is unchanged.
  std::shared_ptr<SimPowerComp<Real>> gen3EMT;
  const auto configureGfl = [&](const auto &inverter) {
    inverter->setParameters(
        gfl.Lf, gfl.Cf, gfl.Rf, gfl.Rc, omegaN,
        opt("gfl_kppll", gfl.KpPLL), opt("gfl_kipll", gfl.KiPLL), omegaN,
        ieee9.gen3.InitialPower, ieee9.gen3.InitialPowerReactive,
        opt("gfl_kpp", gfl.KpPowerCtrl), opt("gfl_kip", gfl.KiPowerCtrl),
        opt("gfl_kpi", gfl.KpCurrCtrl), opt("gfl_kii", gfl.KiCurrCtrl));
  };
  if (gflFormulation == "variable") {
    auto inverter =
        EMT::Ph3::SSN_GFL::make(ieee9.gen3.Name, CPS::Logger::Level::off);
    configureGfl(inverter);
    inverter->setEnableCurrentCrossCoupling(enableGflCurrentCrossCoupling);
    gen3EMT = inverter;
  } else if (gflFormulation == "split") {
    auto inverter = EMT::Ph3::SSN_GFL_Split::make(
        ieee9.gen3.Name, CPS::Logger::Level::off);
    configureGfl(inverter);
    inverter->setEnableCurrentCrossCoupling(enableGflCurrentCrossCoupling);
    gen3EMT = inverter;
  } else if (gflFormulation == "legacy_variable") {
    if (enableGflCurrentCrossCoupling)
      throw std::invalid_argument(
          "Current cross-coupling is unavailable for legacy_variable GFL.");
    auto inverter = EMT::Ph3::AvVoltSourceInverterStateSpace::make(
        ieee9.gen3.Name, CPS::Logger::Level::off);
    configureGfl(inverter);
    gen3EMT = inverter;
  } else {
    throw std::invalid_argument("Unknown gfl_formulation: " +
                                gflFormulation);
  }

  // Loads
  auto load5EMT =
      EMT::Ph3::RXLoad::make(ieee9.load5.Name, CPS::Logger::Level::off);
  load5EMT->setParameters(
      Math::singlePhasePowerToThreePhase(ieee9.load5.RealPower),
      Math::singlePhasePowerToThreePhase(ieee9.load5.ReactivePower),
      std::abs(n5PF->singleVoltage()));

  auto load6EMT =
      EMT::Ph3::RXLoad::make(ieee9.load6.Name, CPS::Logger::Level::off);
  load6EMT->setParameters(
      Math::singlePhasePowerToThreePhase(ieee9.load6.RealPower),
      Math::singlePhasePowerToThreePhase(ieee9.load6.ReactivePower),
      std::abs(n6PF->singleVoltage()));

  auto load8EMT =
      EMT::Ph3::RXLoad::make(ieee9.load8.Name, CPS::Logger::Level::off);
  load8EMT->setParameters(
      Math::singlePhasePowerToThreePhase(load8ActivePower),
      Math::singlePhasePowerToThreePhase(ieee9.load8.ReactivePower),
      std::abs(n8PF->singleVoltage()));

  // Lines
  auto line54EMT =
      EMT::Ph3::PiLine::make(ieee9.line54.Name, CPS::Logger::Level::off);
  line54EMT->setParameters(
      Math::singlePhaseParameterToThreePhase(ieee9.line54.Resistance),
      Math::singlePhaseParameterToThreePhase(ieee9.line54.Inductance),
      Math::singlePhaseParameterToThreePhase(ieee9.line54.Capacitance),
      Math::singlePhaseParameterToThreePhase(ieee9.line54.Conductance));

  auto line64EMT =
      EMT::Ph3::PiLine::make(ieee9.line64.Name, CPS::Logger::Level::off);
  line64EMT->setParameters(
      Math::singlePhaseParameterToThreePhase(ieee9.line64.Resistance),
      Math::singlePhaseParameterToThreePhase(ieee9.line64.Inductance),
      Math::singlePhaseParameterToThreePhase(ieee9.line64.Capacitance),
      Math::singlePhaseParameterToThreePhase(ieee9.line64.Conductance));

  auto line75EMT =
      EMT::Ph3::PiLine::make(ieee9.line75.Name, CPS::Logger::Level::off);
  line75EMT->setParameters(
      Math::singlePhaseParameterToThreePhase(ieee9.line75.Resistance),
      Math::singlePhaseParameterToThreePhase(ieee9.line75.Inductance),
      Math::singlePhaseParameterToThreePhase(ieee9.line75.Capacitance),
      Math::singlePhaseParameterToThreePhase(ieee9.line75.Conductance));

  auto line96EMT =
      EMT::Ph3::PiLine::make(ieee9.line96.Name, CPS::Logger::Level::off);
  line96EMT->setParameters(
      Math::singlePhaseParameterToThreePhase(ieee9.line96.Resistance),
      Math::singlePhaseParameterToThreePhase(ieee9.line96.Inductance),
      Math::singlePhaseParameterToThreePhase(ieee9.line96.Capacitance),
      Math::singlePhaseParameterToThreePhase(ieee9.line96.Conductance));

  auto line78EMT =
      EMT::Ph3::PiLine::make(ieee9.line78.Name, CPS::Logger::Level::off);
  line78EMT->setParameters(
      Math::singlePhaseParameterToThreePhase(ieee9.line78.Resistance),
      Math::singlePhaseParameterToThreePhase(ieee9.line78.Inductance),
      Math::singlePhaseParameterToThreePhase(ieee9.line78.Capacitance),
      Math::singlePhaseParameterToThreePhase(ieee9.line78.Conductance));

  auto line89EMT =
      EMT::Ph3::PiLine::make(ieee9.line89.Name, CPS::Logger::Level::off);
  line89EMT->setParameters(
      Math::singlePhaseParameterToThreePhase(ieee9.line89.Resistance),
      Math::singlePhaseParameterToThreePhase(ieee9.line89.Inductance),
      Math::singlePhaseParameterToThreePhase(ieee9.line89.Capacitance),
      Math::singlePhaseParameterToThreePhase(ieee9.line89.Conductance));

  std::shared_ptr<EMT::Ph3::Switch> line96BreakerEMT;
  if (enableLine96Breaker) {
    if (outagedLine == ieee9.line96.Name)
      throw std::invalid_argument(
          "line96_breaker cannot be combined with outaged_line=LINE96.");
    line96BreakerEMT = EMT::Ph3::Switch::make(
        "LINE96_BREAKER", CPS::Logger::Level::off);
    line96BreakerEMT->setParameters(Matrix::Identity(3, 3) * 1.0e9,
                                    Matrix::Identity(3, 3) * 1.0e-3, true);
  }

  for (const auto &line : {line54EMT, line64EMT, line75EMT, line96EMT,
                           line78EMT, line89EMT}) {
    line->setDefaultParallelConductanceEnabled(enableDefaultLineConductance);
  }

  // Transformers
  auto transf14EMT =
      EMT::Ph3::Transformer::make(ieee9.transf14.Name, CPS::Logger::Level::off);
  transf14EMT->setParameters(
      ieee9.transf14.VoltageLVSide, ieee9.transf14.VoltageHVSide,
      ieee9.transf14.RatedPower, ieee9.transf14.Ratio, 0.0,
      Math::singlePhaseParameterToThreePhase(
          networkEquivalentImpedanceScale * ieee9.transf14.Resistance),
      Math::singlePhaseParameterToThreePhase(
          networkEquivalentImpedanceScale * ieee9.transf14.Inductance));

  auto transf27EMT =
      EMT::Ph3::Transformer::make(ieee9.transf27.Name, CPS::Logger::Level::off);
  transf27EMT->setParameters(
      ieee9.transf27.VoltageLVSide, ieee9.transf27.VoltageHVSide,
      ieee9.transf27.RatedPower, ieee9.transf27.Ratio, 0.0,
      Math::singlePhaseParameterToThreePhase(ieee9.transf27.Resistance),
      Math::singlePhaseParameterToThreePhase(ieee9.transf27.Inductance));

  auto transf39EMT =
      EMT::Ph3::Transformer::make(ieee9.transf39.Name, CPS::Logger::Level::off);
  transf39EMT->setParameters(
      ieee9.transf39.VoltageLVSide, ieee9.transf39.VoltageHVSide,
      ieee9.transf39.RatedPower, ieee9.transf39.Ratio, 0.0,
      Math::singlePhaseParameterToThreePhase(ieee9.transf39.Resistance),
      Math::singlePhaseParameterToThreePhase(ieee9.transf39.Inductance));

  if (useMagnetizingBranches) {
    for (const auto &transformer :
         {transf14EMT, transf27EMT, transf39EMT}) {
      transformer->setMagnetizingBranch(transformerCoreLossPu,
                                         transformerMagnetizingQPu);
    }
  }

  // Connect components to nodes
  gen1EMT->connect({n1EMT});
  // Inverter terminals: 0 = GND, 1 = PCC.
  gen2EMT->connect({SimNode<Real>::GND, n2EMT});
  gen3EMT->connect({SimNode<Real>::GND, n3EMT});

  load5EMT->connect({n5EMT});
  load6EMT->connect({n6EMT});
  load8EMT->connect({n8EMT});

  if (lineInService(line54EMT))
    line54EMT->connect({n5EMT, n4EMT});
  if (lineInService(line64EMT))
    line64EMT->connect({n6EMT, n4EMT});
  if (lineInService(line75EMT))
    line75EMT->connect({n7EMT, n5EMT});
  if (lineInService(line96EMT)) {
    if (line96BreakerEMT) {
      line96BreakerEMT->connect({n9EMT, n96BreakerEMT});
      line96EMT->connect({n96BreakerEMT, n6EMT});
    } else {
      line96EMT->connect({n9EMT, n6EMT});
    }
  }
  if (lineInService(line78EMT))
    line78EMT->connect({n7EMT, n8EMT});
  if (lineInService(line89EMT))
    line89EMT->connect({n8EMT, n9EMT});

  transf14EMT->connect({n1EMT, n4EMT});
  transf27EMT->connect({n2EMT, n7EMT});
  transf39EMT->connect({n3EMT, n9EMT});

  // Create system topology
  SystemComponentList componentsEMT{gen1EMT,   gen2EMT,    gen3EMT,
                                    load5EMT,  load6EMT,   load8EMT,
                                    transf14EMT, transf27EMT, transf39EMT};
  for (const auto &line : {line54EMT, line64EMT, line75EMT, line96EMT,
                           line78EMT, line89EMT}) {
    if (lineInService(line))
      componentsEMT.push_back(line);
  }
  SystemNodeList nodesEMT{n1EMT, n2EMT, n3EMT, n4EMT, n5EMT,
                          n6EMT, n7EMT, n8EMT, n9EMT};
  if (line96BreakerEMT) {
    nodesEMT.push_back(n96BreakerEMT);
    componentsEMT.push_back(line96BreakerEMT);
  }
  auto systemEMT =
      SystemTopology(ieee9.nomFreq, nodesEMT, componentsEMT);

  systemEMT.initWithPowerflow(systemPF, Domain::EMT);
  if (line96BreakerEMT)
    n96BreakerEMT->setInitialVoltage(n9EMT->initialVoltage());

  // Logger
  if (logger) {
    logger->logAttribute("BUS2", n2EMT->attribute("v"));
    logger->logAttribute("BUS3", n3EMT->attribute("v"));
    if (!compactLog) {
      logger->logAttribute("BUS1", n1EMT->attribute("v"));
      logger->logAttribute("BUS4", n4EMT->attribute("v"));
      logger->logAttribute("BUS5", n5EMT->attribute("v"));
      logger->logAttribute("BUS6", n6EMT->attribute("v"));
      logger->logAttribute("BUS7", n7EMT->attribute("v"));
      logger->logAttribute("BUS8", n8EMT->attribute("v"));
      logger->logAttribute("BUS9", n9EMT->attribute("v"));
    }

    if (useNetworkEquivalent && !compactLog) {
      logger->logAttribute("GEN1.I", gen1EMT->attribute("i_intf"));
      logger->logAttribute("GEN1.V", gen1EMT->attribute("v_intf"));
    }

    // GFM inverter (gen2) signals
    if (!compactLog) {
      logger->logAttribute("GEN2.I", gen2EMT->attribute("i_intf"));
      logger->logAttribute("GEN2.V", gen2EMT->attribute("v_intf"));
      logger->logAttribute("GEN2.vc_d", gen2EMT->attribute("vc_d"));
      logger->logAttribute("GEN2.vc_q", gen2EMT->attribute("vc_q"));
    }
    logger->logAttribute("GEN2.p_inst", gen2EMT->attribute("p_inst"));
    logger->logAttribute("GEN2.q_inst", gen2EMT->attribute("q_inst"));
    logger->logAttribute("GEN2.omega", gen2EMT->attribute("omega_gfm"));

    // GFL inverter (gen3) signals
    if (!compactLog) {
      logger->logAttribute("GEN3.I", gen3EMT->attribute("i_intf"));
      logger->logAttribute("GEN3.V", gen3EMT->attribute("v_intf"));
      logger->logAttribute("GEN3.vc_d", gen3EMT->attribute("vc_d"));
      logger->logAttribute("GEN3.vc_q", gen3EMT->attribute("vc_q"));
    }
    logger->logAttribute("GEN3.p_inst", gen3EMT->attribute("p_inst"));
    logger->logAttribute("GEN3.q_inst", gen3EMT->attribute("q_inst"));
    logger->logAttribute("GEN3.omega_pll", gen3EMT->attribute("omega_pll"));

    if (!compactLog) {
      // log generator's current
      for (auto comp : systemEMT.mComponents) {
        if (std::dynamic_pointer_cast<
                CPS::EMT::Ph3::SynchronGenerator4OrderVBR>(comp)) {
          logger->logAttribute(comp->name() + ".I",
                               comp->attribute("i_intf"));
          logger->logAttribute(comp->name() + ".V",
                               comp->attribute("v_intf"));
          logger->logAttribute(comp->name() + ".omega",
                               comp->attribute("w_r"));
          logger->logAttribute(comp->name() + ".delta",
                               comp->attribute("delta"));
        }
      }

      // log transformer voltages and currents
      for (auto comp : systemEMT.mComponents) {
        if (std::dynamic_pointer_cast<CPS::EMT::Ph3::Transformer>(comp)) {
          logger->logAttribute(comp->name() + ".I",
                               comp->attribute("i_intf"));
          logger->logAttribute(comp->name() + ".V",
                               comp->attribute("v_intf"));
        }
      }

      // log line voltages and currents
      for (auto comp : systemEMT.mComponents) {
        if (std::dynamic_pointer_cast<CPS::EMT::Ph3::PiLine>(comp)) {
          logger->logAttribute(comp->name() + ".I",
                               comp->attribute("i_intf"));
          logger->logAttribute(comp->name() + ".V",
                               comp->attribute("v_intf"));
        }
      }
    }
  }

  return systemEMT;
}

#ifndef DPSIM_IEEE9_INVERTER_MIX_LIBRARY
int main(int argc, char *argv[]) {

  CommandLineArgs args(argc, argv, "EMT_Ph3_IEEE9_SSN_InverterMix", 0.00005,
                       0.01 * 60, 60, -1, CPS::Logger::Level::info,
                       CPS::Logger::Level::off, false, false, false,
                       CPS::Domain::EMT);

  CPS::Logger::setLogDir("./logs/" + args.name);
  bool log = args.options.find("log") != args.options.end() &&
             args.getOptionBool("log");

  std::filesystem::path logFilename =
      "./logs/" + args.name + "/" + args.name + ".csv";
  std::shared_ptr<DataLoggerInterface> logger = nullptr;

  if (log) {
    logger =
        RealTimeDataLogger::make(logFilename, args.duration, args.timeStep);
  }

  auto sys = buildTopology(args, logger);

  Simulation sim(args.name, args);
  sim.setSystem(sys);
  sim.setDomain(Domain::EMT);
  sim.doSystemMatrixRecomputation(true);
  if (log) {
    sim.addLogger(logger);
  }
  sim.run();

  CPS::Logger::get(args.name)->info("Simulation finished.");
}
#endif
