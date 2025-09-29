#include <DPsim.h>
#include <iostream>
#include <string>
#include <vector>

using namespace DPsim;
using namespace CPS::EMT;
using namespace CPS::EMT::Ph1;

void syncronousGeneratorOldTest(Real timeStep, Real finalTime,
                             bool doEigenvalueExtraction) {

  // Nodes
  auto n1 = SimNode::make("n1");
  auto n2 = SimNode::make("n2");
  auto n3 = SimNode::make("n3");
  auto n4 = SimNode::make("n4");
  auto n5 = SimNode::make("n5");
  auto n6 = SimNode::make("n6");
  auto n7 = SimNode::make("n7");
  auto n8 = SimNode::make("n8");
  auto n9 = SimNode::make("n9");
  auto n10 = SimNode::make("n10");
  auto n11 = SimNode::make("n11");
  auto n12 = SimNode::make("n12");
  auto n13 = SimNode::make("n13");
  auto n14 = SimNode::make("n14");
  auto n15 = SimNode::make("n15");
  auto n16 = SimNode::make("n16");
  auto n17 = SimNode::make("n17");
  auto n18 = SimNode::make("n18");
  auto n19 = SimNode::make("n19");

  //Components
  auto r1 = Resistor::make("r1", Logger::Level::debug);
  r1->setParameters(100);
  auto r2 = Resistor::make("r2", Logger::Level::debug);
  r2->setParameters(100);
  auto r3 = Resistor::make("r3", Logger::Level::debug);
  r3->setParameters(100);

  auto R0 =
      Resistor::make("R0", Logger::Level::debug); // zero sequence resistor
  R0->setParameters(0.00311);

  auto L0 =
      Inductor::make("L0", Logger::Level::debug); // zero sequence inductor
  L0->setParameters(4.9553e-4);

  // d-axis components
  auto Rsd = Resistor::make(
      "Rsd", Logger::Level::debug); //stator winding resistor d-axis
  Rsd->setParameters(0.00311);
  auto Ls_sigma_d = Inductor::make(
      "Ls_sigma", Logger::Level::debug); // leakage stator inductance d-axis
  Ls_sigma_d->setParameters(4.9553e-4);
  auto Lmd =
      Inductor::make("Lmd", Logger::Level::debug); // mutual inductance d-axis
  Lmd->setParameters(0.00548);
  auto L1d = Inductor::make(
      "L1_d", Logger::Level::debug); // d-axis damping winding inductance
  L1d->setParameters(5.659e-4);
  auto R1d = Resistor::make(
      "R1d", Logger::Level::debug); // d-axis damping winding resistor
  R1d->setParameters(0.02947);
  auto Lf_d =
      Inductor::make("Lf_d", Logger::Level::debug); // Field winding inductance
  Lf_d->setParameters(5.451e-4);
  auto Rf_d =
      Resistor::make("Rf_d", Logger::Level::debug); // Field winding resistor
  Rf_d->setParameters(6.227e-4);
  auto Vf_d = VoltageSource::make(
      "Vf_d", Logger::Level::debug); // Generator excitation voltage
  Vf_d->setParameters(Complex(7.0829, 0.0));

  // q-axis components
  auto Rsq = Resistor::make(
      "Rsq", Logger::Level::debug); //stator winding resistor q-axis
  Rsq->setParameters(0.00311);
  auto Ls_sigma_q = Inductor::make(
      "Lssigma_q", Logger::Level::debug); // leakage stator inductance q-axis
  Ls_sigma_q->setParameters(4.9553e-4);
  auto Lmq =
      Inductor::make("Lmq", Logger::Level::debug); // mutual inductance q-axis
  Lmq->setParameters(0.00532);
  auto L1q = Inductor::make(
      "L1q", Logger::Level::debug); // q-axis 1st damping winding inductance
  L1q->setParameters(0.00239);
  auto R1q = Resistor::make(
      "R1q", Logger::Level::debug); // q-axis 1st damping winding resistor
  R1q->setParameters(0.00642);
  auto L2q = Inductor::make(
      "L2q", Logger::Level::debug); // q-axis 2nd damping winding inductance
  L2q->setParameters(4.129e-4);
  auto R2q = Resistor::make(
      "R2q", Logger::Level::debug); // q-axis 2nd damping winding resistor
  R2q->setParameters(0.02458);

  auto ElectroMechanicalConverter_d = ElectroMechanicalConverter::make(
      "ElectroMechanicalConverter_d", Logger::Level::debug);
  ElectroMechanicalConverter_d->setVoltageReferenceNode(n13);

  auto ElectroMechanicalConverter_q = ElectroMechanicalConverter::make(
      "ElectroMechanicalConverter_q", Logger::Level::debug);
  ElectroMechanicalConverter_q->setVoltageReferenceNode(n6);

  auto ParkTrafo =
      ParkTransformer::make("ParkTransformer", Logger::Level::debug);
  ParkTrafo->setIsOmegaConstant(
      true); // Set the Park transformer to use a constant angular frequency
  ParkTrafo->setInitialValues(
      2 * M_PI * 50,
      M_PI / 2); // Set the angular frequency and initial angle

  //Mechanical components
  auto InertiaMoment =
      InertiaMoment::make("InertiaMoment", Logger::Level::debug);
  InertiaMoment->setParameters(0.0325); // Set the inertia moment value

  auto ConstantOmegaSource =
      VoltageSource::make("ConstantOmegaSource", Logger::Level::debug);
  ConstantOmegaSource->setParameters(
      Complex(2 * M_PI * 50, 0.0)); // Set the constant angular frequency source

  //Topology
  // abc side
  r1->connect({n1, SimNode::GND});
  r2->connect({n2, SimNode::GND});
  r3->connect({n3, SimNode::GND});
  ParkTrafo->connect({n1, n2, n3, n4, n11, n17});

  // d-axis circuit
  Rsd->connect({n4, n5});
  ElectroMechanicalConverter_d->connect({n6, n5, n18, SimNode::GND});
  Ls_sigma_d->connect({n6, n7});
  Lmd->connect({n7, SimNode::GND});
  L1d->connect({n7, n8});
  R1d->connect({n8, SimNode::GND});
  Lf_d->connect({n7, n9});
  Rf_d->connect({n9, n10});
  Vf_d->connect({SimNode::GND, n10});

  // q-axis circuit
  Rsq->connect({n11, n12});
  ElectroMechanicalConverter_q->connect({n12, n13, n18, SimNode::GND});
  Ls_sigma_q->connect({n13, n14});
  Lmq->connect({n14, SimNode::GND});
  L1q->connect({n14, n15});
  R1q->connect({n15, SimNode::GND});
  L2q->connect({n14, n16});
  R2q->connect({n16, SimNode::GND});

  // 0 axis circuit
  R0->connect({n17, n19});
  L0->connect({n19, SimNode::GND});

  //Mechanical side
  InertiaMoment->connect({n18, SimNode::GND});
  ConstantOmegaSource->connect({SimNode::GND, n18});

  // Define system topology
  SystemTopology system(50,
                        SystemNodeList{n1, n2, n3, n4, n5, n6, n7, n8, n9, n10,
                                       n11, n12, n13, n14, n15, n16, n17, n18,
                                       n19},
                        SystemComponentList{R1d,
                                            Rsd,
                                            Ls_sigma_d,
                                            Lmd,
                                            L1d,
                                            R1q,
                                            Rsq,
                                            Ls_sigma_q,
                                            Lmq,
                                            L1q,
                                            L2q,
                                            R2q,
                                            ElectroMechanicalConverter_d,
                                            ElectroMechanicalConverter_q,
                                            ParkTrafo,
                                            Vf_d,
                                            Rf_d,
                                            InertiaMoment,
                                            r1,
                                            r2,
                                            r3,
                                            R0,
                                            L0,
                                            Lf_d,
                                            ConstantOmegaSource});

  // Define simulation scenario
  String simName = "syncrhonousGeneratorTest";

  // Logger
  auto logger = DataLogger::make(simName);
  logger->logAttribute("omega", n18->attribute("v"));
  //   logger->logAttribute("Torque", InertiaMoment->attribute("i_intf"));
  logger->logAttribute("Ratio_d",
                       ElectroMechanicalConverter_d->attribute("flux"));
  logger->logAttribute("Ratio_q",
                       ElectroMechanicalConverter_q->attribute("flux"));
  //   logger->logAttribute("V_f_d", n10->attribute("v"));
  //   logger->logAttribute("I_f_d", Vf_d->attribute("i_intf"));
  logger->logAttribute("V_A", n1->attribute("v"));
  logger->logAttribute("V_B", n2->attribute("v"));
  logger->logAttribute("V_C", n3->attribute("v"));
  logger->logAttribute("I_A", r1->attribute("i_intf"));
  logger->logAttribute("I_B", r1->attribute("i_intf"));
  logger->logAttribute("I_C", r1->attribute("i_intf"));
  //   logger->logAttribute("V_Q", n11->attribute("v"));
  //   logger->logAttribute("V_D", n4->attribute("v"));
  //   logger->logAttribute("V_Q", n11->attribute("v"));
  //   logger->logAttribute("V_DS", n6->attribute("v"));
  //   logger->logAttribute("V_QS", n13->attribute("v"));
  //   logger->logAttribute("V_0", n17->attribute("v"));
  //   logger->logAttribute("I_sd", Rsd->attribute("i_intf"));
  //   logger->logAttribute("I_sq", Rsq->attribute("i_intf"));

  Simulation sim(simName);
  sim.setSystem(system);
  sim.setTimeStep(timeStep);
  sim.doSystemMatrixRecomputation(true);
  sim.setFinalTime(finalTime);
  sim.setDomain(Domain::EMT);
  sim.addLogger(logger);
  sim.doEigenvalueExtraction(doEigenvalueExtraction);
  sim.run();
}

void electromechanicalConverterTest(Real timeStep, Real finalTime,
                                    bool doEigenvalueExtraction) {

  auto n1 = SimNode::make("n1");
  auto n2 = SimNode::make("n2");
  auto n3 = SimNode::make("n3");
  auto n4 = SimNode::make("n4");

  auto v = VoltageSource::make("v", Logger::Level::debug);
  v->setParameters(Complex(10.0, 0.0), 0);
  v->connect({SimNode::GND, n1});

  auto l = Inductor::make("l", Logger::Level::debug);
  l->setParameters(0.1);
  l->connect({n1, n2});

  auto r = Resistor::make("r", Logger::Level::debug);
  r->setParameters(1.0);
  r->connect({n3, SimNode::GND});

  auto converter =
      ElectroMechanicalConverter::make("converter", Logger::Level::debug);
  converter->setVoltageReferenceNode(n1);
  converter->connect({n3, n2, n4, SimNode::GND});

  auto inertiaMoment =
      InertiaMoment::make("inertiaMoment", Logger::Level::debug);
  inertiaMoment->setParameters(1e-3);
  inertiaMoment->connect({n4, SimNode::GND});

  // Define system topology
  SystemTopology system(50, SystemNodeList{n1, n2, n3, n4},
                        SystemComponentList{v, l, r, converter, inertiaMoment});

  // Define simulation scenario
  String simName = "electromechanicalConverterTest";
  // Logger
  auto logger = DataLogger::make(simName);
  logger->logAttribute("Flux", converter->attribute("flux"));
  logger->logAttribute("V1", n1->attribute("v"));
  logger->logAttribute("V2", n2->attribute("v"));
  logger->logAttribute("V3", n3->attribute("v"));
  logger->logAttribute("V4", n4->attribute("v"));
  logger->logAttribute("I12", r->attribute("i_intf"));

  Simulation sim(simName);
  sim.setSystem(system);
  sim.setTimeStep(timeStep);
  sim.doSystemMatrixRecomputation(true);
  sim.setFinalTime(finalTime);
  sim.setDomain(Domain::EMT);
  sim.addLogger(logger);
  sim.doEigenvalueExtraction(doEigenvalueExtraction);
  sim.run();
}

void speedVoltageTermTest(Real timeStep, Real finalTime,
                          bool doEigenvalueExtraction) {

  // Paramters
  Real frequency = 50;
  Real omega = 2 * M_PI * frequency;
  Real l = 0.2;                 // Inductance
  Real r = 1.0;                 // Resistance
  Real voltageMagnitude = 10.0; // Voltage magnitude

  // Nodes
  auto n1 = SimNode::make("n1");
  auto n2 = SimNode::make("n2");
  auto n3 = SimNode::make("n3");
  auto n4 = SimNode::make("n4");
  auto n5 = SimNode::make("n5");
  auto n6 = SimNode::make("n6");
  auto n7 = SimNode::make("n7");
  auto n8 = SimNode::make("n8");

  // Components
  auto vD = VoltageSource::make("vD", Logger::Level::debug);
  vD->setParameters(Complex(voltageMagnitude, 0.0), frequency);
  vD->connect({SimNode::GND, n1});

  auto vQ = VoltageSource::make("vQ", Logger::Level::debug);
  vQ->setParameters(Complex(voltageMagnitude * cos(-M_PI / 2),
                            voltageMagnitude * sin(-M_PI / 2)),
                    frequency);
  vQ->connect({SimNode::GND, n5});

  auto resistorD = Resistor::make("resistorD", Logger::Level::debug);
  resistorD->setParameters(r);
  resistorD->connect({n1, n2});

  auto inductorD = Inductor::make("inductorD", Logger::Level::debug);
  inductorD->setParameters(l);
  inductorD->connect({n4, SimNode::GND});

  auto voltageSpeedTermD =
      VoltageSpeedTerm::make("voltageSpeedTermD", Logger::Level::debug);
  voltageSpeedTermD->setInitialOmega(omega);
  voltageSpeedTermD->setInductance(l);
  voltageSpeedTermD->setIsNegative(true);
  voltageSpeedTermD->setIsConstantSpeed(true);
  voltageSpeedTermD->connect({n7, n8, n2, n3});

  auto resistorQ = Resistor::make("resistorQ", Logger::Level::debug);
  resistorQ->setParameters(r);
  resistorQ->connect({n5, n6});

  auto inductorQ = Inductor::make("inductorQ", Logger::Level::debug);
  inductorQ->setParameters(l);
  inductorQ->connect({n8, SimNode::GND});

  auto voltageSpeedTermQ =
      VoltageSpeedTerm::make("voltageSpeedTermQ", Logger::Level::debug);
  voltageSpeedTermQ->setInitialOmega(omega);
  voltageSpeedTermQ->setInductance(l);
  voltageSpeedTermQ->setIsNegative(false);
  voltageSpeedTermQ->setIsConstantSpeed(true);
  voltageSpeedTermQ->connect({n3, n4, n6, n7}); // through other speed term

  // Define system topology
  SystemTopology system(50, SystemNodeList{n1, n2, n3, n4, n5, n6, n7, n8},
                        SystemComponentList{vD, vQ, resistorD, inductorD,
                                            voltageSpeedTermD, resistorQ,
                                            inductorQ, voltageSpeedTermQ});

  // Define simulation scenario
  String simName = "speedVoltageTermTest";
  // Logger
  auto logger = DataLogger::make(simName);
  logger->logAttribute("V_D", n1->attribute("v"));
  logger->logAttribute("V_Q", n2->attribute("v"));
  logger->logAttribute("I_D", resistorD->attribute("i_intf"));
  logger->logAttribute("I_Q", resistorQ->attribute("i_intf"));
  logger->logAttribute("V_SpeedTerm_D", n3->attribute("v"));
  logger->logAttribute("V_SpeedTerm_Q", n4->attribute("v"));
  logger->logAttribute("I_SpeedTerm_D", voltageSpeedTermD->attribute("i_intf"));
  logger->logAttribute("I_SpeedTerm_Q", voltageSpeedTermQ->attribute("i_intf"));

  Simulation sim(simName);
  sim.setSystem(system);
  sim.setTimeStep(timeStep);
  sim.doSystemMatrixRecomputation(true);
  sim.setFinalTime(finalTime);
  sim.setDomain(Domain::EMT);
  sim.addLogger(logger);
  sim.doEigenvalueExtraction(doEigenvalueExtraction);
  sim.run();
}

void timeLaggingVoltageSourceTest() {
  // Parameters
  Real frequency = 50;
  Real r = 1.0;                 // Resistance
  Real voltageMagnitude = 10.0; // Voltage magnitude

  // Nodes
  auto n1 = SimNode::make("n1");
  auto n2 = SimNode::make("n2");
  auto n3 = SimNode::make("n3");

  // Components
  auto vS = VoltageSource::make("vS", Logger::Level::debug);
  vS->setParameters(Complex(voltageMagnitude, 0.0), frequency);
  vS->connect({SimNode::GND, n1});

  auto r1 = Resistor::make("r1", Logger::Level::debug);
  r1->setParameters(r);
  r1->connect({n1, n2});

  auto r2 = Resistor::make("r2", Logger::Level::debug);
  r2->setParameters(r);
  r2->connect({n2, n3});

  auto tLVS = TimeLaggingVoltageSource::make("tLVS", Logger::Level::debug);
  tLVS->setVoltageReferenceNodes(n1, n2);
  tLVS->connect({n3, SimNode::GND});

  // Define system topology
  SystemTopology system(50, SystemNodeList{n1, n2, n3},
                        SystemComponentList{vS, r1, r2, tLVS});

  // Define simulation scenario
  String simName = "timeLaggingVoltageSourceTest";
  // Logger
  auto logger = DataLogger::make(simName);
  logger->logAttribute("V1", n1->attribute("v"));
  logger->logAttribute("V2", n2->attribute("v"));
  logger->logAttribute("V3", n3->attribute("v"));
  logger->logAttribute("I_r1", r1->attribute("i_intf"));
  logger->logAttribute("V_TLVS", tLVS->attribute("v"));

  Simulation sim(simName);
  sim.setSystem(system);
  sim.setTimeStep(1e-4);
  sim.doSystemMatrixRecomputation(false);
  sim.setFinalTime(0.1);
  sim.setDomain(Domain::EMT);
  sim.addLogger(logger);
  sim.doEigenvalueExtraction(false);
  sim.run();
}

void timeLaggingTorqueSourceTest() {
  // Parameters
  Real frequency = 50;
  Real r = 1.0;                 // Resistance
  Real voltageMagnitude = 10.0; // Voltage magnitude

  // Nodes
  auto n1 = SimNode::make("n1");
  auto n2 = SimNode::make("n2");
  auto n3 = SimNode::make("n3");

  // Components
  auto vS = VoltageSource::make("vS", Logger::Level::debug);
  vS->setParameters(Complex(voltageMagnitude, 0.0), frequency);
  vS->connect({SimNode::GND, n1});

  auto r1 = Resistor::make("r1", Logger::Level::debug);
  r1->setParameters(r);
  r1->connect({n2, SimNode::GND});

  auto r2 = Resistor::make("r2", Logger::Level::debug);
  r2->setParameters(r);
  r2->connect({n3, SimNode::GND});

  auto tLTS = TimeLaggingTorqueSource::make("tLTS", Logger::Level::debug);
  tLTS->connect({n1, n2, SimNode::GND, n3});

  // Define system topology
  SystemTopology system(50, SystemNodeList{n1, n2, n3},
                        SystemComponentList{vS, r1, r2, tLTS});

  // Define simulation scenario
  String simName = "timeLaggingTorqueSourceTest";
  // Logger
  auto logger = DataLogger::make(simName);
  logger->logAttribute("V1", n1->attribute("v"));
  logger->logAttribute("V2", n2->attribute("v"));
  logger->logAttribute("V3", n3->attribute("v"));
  logger->logAttribute("I_r1", r1->attribute("i_intf"));
  logger->logAttribute("I_r2", r2->attribute("i_intf"));
  logger->logAttribute("I_TLTS", tLTS->attribute("i"));

  Simulation sim(simName);
  sim.setSystem(system);
  sim.setTimeStep(1e-4);
  sim.doSystemMatrixRecomputation(false);
  sim.setFinalTime(0.1);
  sim.setDomain(Domain::EMT);
  sim.addLogger(logger);
  sim.doEigenvalueExtraction(false);
  sim.run();
}

void motorStartingTest(Real timeStep, Real finalTime,
                       bool doEigenvalueExtraction) {
  Real frequency = 50;
  Real omega = 2 * M_PI * frequency;
  Real voltageMagnitudeLL = 10e3;
  Real rS = 2.0995;
  Real xS = 6.9115;
  Real xM = 233.67;
  Real rR = 556.9 * 1e-3;
  Real xR = 6.9115;
  Real inertia = 10;

  Real voltageMagnitude = voltageMagnitudeLL * sqrt(2.0) / sqrt(3.0);
  Complex voltageL1 = Complex(voltageMagnitude, 0.0);
  Complex voltageL2 = Complex(voltageMagnitude * cos(-2 * M_PI / 3),
                              voltageMagnitude * sin(-2 * M_PI / 3));
  Complex voltageL3 = Complex(voltageMagnitude * cos(2 * M_PI / 3),
                              voltageMagnitude * sin(2 * M_PI / 3));
  Real lS = xS / omega;
  Real lM = xM / omega;
  Real lR = xR / omega;

  // Nodes
  auto n1 = SimNode::make("n1");
  auto n2 = SimNode::make("n2");
  auto n3 = SimNode::make("n3");
  auto n4 = SimNode::make("n4");
  auto n5 = SimNode::make("n5");
  auto n6 = SimNode::make("n6");
  auto n7 = SimNode::make("n7");
  auto n8 = SimNode::make("n8");
  auto n9 = SimNode::make("n9");
  auto n10 = SimNode::make("n10");
  auto n11 = SimNode::make("n11");
  auto n12 = SimNode::make("n12");
  auto n13 = SimNode::make("n13");
  auto n14 = SimNode::make("n14");
  auto n15 = SimNode::make("n15");
  auto n16 = SimNode::make("n16");
  auto n17 = SimNode::make("n17");
  auto n18 = SimNode::make("n18");
  auto n19 = SimNode::make("n19");
  auto n20 = SimNode::make("n20");
  auto n21 = SimNode::make("n21");
  auto n22 = SimNode::make("n22");
  auto n23 = SimNode::make("n23");
  auto n24 = SimNode::make("n24");
  auto n25 = SimNode::make("n25");
  auto n26 = SimNode::make("n26");
  auto n27 = SimNode::make("n27");
  auto n28 = SimNode::make("n28");
  auto n29 = SimNode::make("n29");
  auto n30 = SimNode::make("n30");
  auto n31 = SimNode::make("n31");
  auto n32 = SimNode::make("n32");
  auto n33 = SimNode::make("n33");
  auto n34 = SimNode::make("n34");
  auto n35 = SimNode::make("n35");
  auto n36 = SimNode::make("n36");
  auto n37 = SimNode::make("n37");
  auto n38 = SimNode::make("n38");
  auto n39 = SimNode::make("n39");
  auto n40 = SimNode::make("n40");
  auto n41 = SimNode::make("n41");
  auto n42 = SimNode::make("n42");
  auto n43 = SimNode::make("n43");
  auto n44 = SimNode::make("n44");
  auto n45 = SimNode::make("n45");
  auto n46 = SimNode::make("n46");
  auto n47 = SimNode::make("n47");

  // Components
  // Infeed
  auto v1 = VoltageSource::make("v1", Logger::Level::debug);
  v1->setParameters(voltageL1, frequency);
  v1->connect({SimNode::GND, n1});
  auto v2 = VoltageSource::make("v2", Logger::Level::debug);
  v2->setParameters(voltageL2, frequency);
  v2->connect({SimNode::GND, n2});
  auto v3 = VoltageSource::make("v3", Logger::Level::debug);
  v3->setParameters(voltageL3, frequency);
  v3->connect({SimNode::GND, n3});

  // Park transformer
  auto parkTrafo =
      ParkTransformer::make("ParkTransformer", Logger::Level::debug);
  parkTrafo->setIsOmegaConstant(true);
  parkTrafo->setInitialValues(omega, 0.0);
  parkTrafo->connect({n1, n2, n3, n4, n5, n6});

  // d-axis components
  auto rSd = Resistor::make("rSd", Logger::Level::debug);
  rSd->setParameters(rS);
  rSd->connect({n4, n7});

  auto ccvs1 = VoltageSpeedTerm::make("ccvs1", Logger::Level::debug);
  ccvs1->setInitialOmega(omega);
  ccvs1->setIsConstantSpeed(true);
  ccvs1->setInductance(lS);
  ccvs1->setIsNegative(false);
  ccvs1->connect({n13, n15, n10, n7});

  auto ccvs2 = VoltageSpeedTerm::make("ccvs2", Logger::Level::debug);
  ccvs2->setInitialOmega(omega);
  ccvs2->setIsConstantSpeed(true);
  ccvs2->setInductance(lM);
  ccvs2->setIsNegative(false);
  ccvs2->connect({n19, n21, n12, n10});

  auto lSd = Inductor::make("lSd", Logger::Level::debug);
  lSd->setParameters(lS);
  lSd->connect({n14, n16});

  auto lMd = Inductor::make("lMd", Logger::Level::debug);
  lMd->setParameters(lM);
  lMd->connect({n16, n18});

  auto lRd = Inductor::make("lRd", Logger::Level::debug);
  lRd->setParameters(lR);
  lRd->connect({n26, n16});

  auto ccvs5 = VoltageSpeedTerm::make("ccvs5", Logger::Level::debug);
  ccvs5->setInitialOmega(omega);
  ccvs5->setIsConstantSpeed(true);
  ccvs5->setInductance(lR);
  ccvs5->setIsNegative(false);
  ccvs5->connect({n39, n37, n26, n28});

  auto ccvs6 = VoltageSpeedTerm::make("ccvs6", Logger::Level::debug);
  ccvs6->setInitialOmega(omega);
  ccvs6->setIsConstantSpeed(true);
  ccvs6->setInductance(lM);
  ccvs6->setIsNegative(false);
  ccvs6->connect({n21, n23, n28, n30});

  auto rec1 = ElectroMechanicalConverter::make("rec1", Logger::Level::debug);
  rec1->setInitialFlux(0.0);
  rec1->setIsNegative(false);
  rec1->setVoltageReferenceNode(n27);
  rec1->connect({n32, n30, n46, SimNode::GND});

  auto ccvs9 = VoltageSpeedTerm::make("ccvs9", Logger::Level::debug);
  ccvs9->setInitialOmega(0);
  ccvs9->setIsConstantSpeed(false);
  ccvs9->setOmegaReferenceNode(n46);
  ccvs9->setInductance(lM);
  ccvs9->setIsNegative(false);
  ccvs9->connect({n23, n25, n34, n32});

  auto ccvs10 = VoltageSpeedTerm::make("ccvs10", Logger::Level::debug);
  ccvs10->setInitialOmega(0);
  ccvs10->setIsConstantSpeed(false);
  ccvs10->setOmegaReferenceNode(n46);
  ccvs10->setInductance(lR);
  ccvs10->setIsNegative(false);
  ccvs10->connect({n41, n39, n36, n34});

  auto rRd = Resistor::make("rRd", Logger::Level::debug);
  rRd->setParameters(rR);
  rRd->connect({n44, SimNode::GND});

  auto tlvs1 = TimeLaggingVoltageSource::make("tlvs1", Logger::Level::debug);
  tlvs1->setVoltageReferenceNodes(n44, n30);
  tlvs1->connect({n42, n44});

  // q-axis components
  auto rSq = Resistor::make("rSq", Logger::Level::debug);
  rSq->setParameters(rS);
  rSq->connect({n5, n8});

  auto ccvs3 = VoltageSpeedTerm::make("ccvs3", Logger::Level::debug);
  ccvs3->setInitialOmega(omega);
  ccvs3->setIsConstantSpeed(true);
  ccvs3->setInductance(lS);
  ccvs3->setIsNegative(false);
  ccvs3->connect({n12, n14, n8, n11});

  auto ccvs4 = VoltageSpeedTerm::make("ccvs4", Logger::Level::debug);
  ccvs4->setInitialOmega(omega);
  ccvs4->setIsConstantSpeed(true);
  ccvs4->setInductance(lM);
  ccvs4->setIsNegative(false);
  ccvs4->connect({n18, n20, n11, n13});

  auto lSq = Inductor::make("lSq", Logger::Level::debug);
  lSq->setParameters(lS);
  lSq->connect({n15, n17});

  auto lMq = Inductor::make("lMq", Logger::Level::debug);
  lMq->setParameters(lM);
  lMq->connect({n17, n19});

  auto lRq = Inductor::make("lRq", Logger::Level::debug);
  lRq->setParameters(lR);
  lRq->connect({n27, n17});

  auto ccvs7 = VoltageSpeedTerm::make("ccvs7", Logger::Level::debug);
  ccvs7->setInitialOmega(omega);
  ccvs7->setIsConstantSpeed(true);
  ccvs7->setInductance(lR);
  ccvs7->setIsNegative(false);
  ccvs7->connect({n38, n36, n29, n27});

  auto ccvs8 = VoltageSpeedTerm::make("ccvs8", Logger::Level::debug);
  ccvs8->setInitialOmega(omega);
  ccvs8->setIsConstantSpeed(true);
  ccvs8->setInductance(lM);
  ccvs8->setIsNegative(false);
  ccvs8->connect({n20, n22, n31, n29});

  auto rec2 = ElectroMechanicalConverter::make("rec2", Logger::Level::debug);
  rec2->setInitialFlux(0.0);
  rec2->setIsNegative(false);
  rec2->setVoltageReferenceNode(n26);
  rec2->connect({n31, n33, n46, SimNode::GND});

  auto ccvs11 = VoltageSpeedTerm::make("ccvs11", Logger::Level::debug);
  ccvs11->setInitialOmega(0);
  ccvs11->setIsConstantSpeed(false);
  ccvs11->setOmegaReferenceNode(n46);
  ccvs11->setInductance(lM);
  ccvs11->setIsNegative(false);
  ccvs11->connect({n22, n24, n33, n35});

  auto ccvs12 = VoltageSpeedTerm::make("ccvs12", Logger::Level::debug);
  ccvs12->setInitialOmega(0);
  ccvs12->setIsConstantSpeed(false);
  ccvs12->setOmegaReferenceNode(n46);
  ccvs12->setInductance(lR);
  ccvs12->setIsNegative(false);
  ccvs12->connect({n40, n38, n35, n37});

  auto rRq = Resistor::make("rRq", Logger::Level::debug);
  rRq->setParameters(rR);
  rRq->connect({n45, SimNode::GND});

  auto tlvs2 = TimeLaggingVoltageSource::make("tlvs2", Logger::Level::debug);
  tlvs2->setVoltageReferenceNodes(n31, n45);
  tlvs2->connect({n45, n43});

  // 0-axis components
  auto r0 = Resistor::make("r0", Logger::Level::debug);
  r0->setParameters(rS);
  r0->connect({n6, n9});

  auto l0 = Inductor::make("l0", Logger::Level::debug);
  l0->setParameters(lS);
  l0->connect({n9, SimNode::GND});

  // Mechanical components
  auto inertiaMoment =
      InertiaMoment::make("inertiaMoment", Logger::Level::debug);
  inertiaMoment->setParameters(inertia);
  inertiaMoment->connect({n47, SimNode::GND});

  auto ccts1 =
      CurrentControlledTorqueSource::make("ccts1", Logger::Level::debug);
  ccts1->setInductance(lM);
  ccts1->setInductor(lRd);
  ccts1->connect({n25, SimNode::GND, SimNode::GND, n46});

  auto ccts2 =
      CurrentControlledTorqueSource::make("ccts2", Logger::Level::debug);
  ccts2->setInductance(lR);
  ccts2->setInductor(lRd);
  ccts2->connect({n43, n41, SimNode::GND, n46});

  auto ccts3 =
      CurrentControlledTorqueSource::make("ccts3", Logger::Level::debug);
  ccts3->setInductance(lM);
  ccts3->setInductor(lRq);
  ccts3->connect({n24, SimNode::GND, n46, SimNode::GND});

  auto ccts4 =
      CurrentControlledTorqueSource::make("ccts4", Logger::Level::debug);
  ccts4->setInductance(lR);
  ccts4->setInductor(lRq);
  ccts4->connect({n42, n40, n46, SimNode::GND});

  auto tlts = TimeLaggingTorqueSource::make("tlts", Logger::Level::debug);
  tlts->connect({n46, n47, n46, SimNode::GND});

  // Define system topology
  SystemTopology system(
      50,
      SystemNodeList{n1,  n2,  n3,  n4,  n5,  n6,  n7,  n8,  n9,  n10, n11, n12,
                     n13, n14, n15, n16, n17, n18, n19, n20, n21, n22, n23, n24,
                     n25, n26, n27, n28, n29, n30, n31, n32, n33, n34, n35, n36,
                     n37, n38, n39, n40, n41, n42, n43, n44, n45, n46, n47},
      SystemComponentList{
          v1,     v2,     v3,    parkTrafo, rSd,           lSd,   lMd,
          lRd,    rRd,    ccvs1, ccvs2,     ccvs5,         ccvs6, rec1,
          ccvs9,  ccvs10, tlvs1, rSq,       lSq,           lMq,   lRq,
          rRq,    ccvs3,  ccvs4, ccvs7,     ccvs8,         rec2,  ccvs11,
          ccvs12, tlvs2,  r0,    l0,        inertiaMoment, ccts1, ccts2,
          ccts3,  ccts4,  tlts});

  // Define simulation scenario
  String simName = "motorStartingTest";
  // Logger
  auto logger = DataLogger::make(simName);
  logger->logAttribute("i_A", v1->attribute("i_intf"));
  logger->logAttribute("omega", n47->attribute("v"));
  logger->logAttribute("torque", inertiaMoment->attribute("i_intf"));

  Simulation sim(simName);
  sim.setSystem(system);
  sim.setTimeStep(timeStep);
  sim.doSystemMatrixRecomputation(true);
  sim.setFinalTime(finalTime);
  sim.setDomain(Domain::EMT);
  sim.addLogger(logger);
  sim.doEigenvalueExtraction(doEigenvalueExtraction);
  
//   inertiaMoment->setIntfVoltage(-omega);
  sim.run();
}

void synchronousGeneratorTest(Real timeStep, Real finalTime,
                       bool doEigenvalueExtraction) {
  Real frequency = 50;
  Real omega = 2 * M_PI * frequency;
  Real voltageMagnitudeLL = 24e3;
  Real rS = 0.00311;
  Real lS = 4.9553e-4;
  Real lMd = 0.00548;
  Real lMq = 0.00532;
  Real rFd = 6.227e-4;
  Real lFd = 5.451e-4;
  Real r1d = 0.02947;
  Real l1d = 5.659e-4;
  Real r1q = 0.00642;
  Real l1q = 0.00239;
  Real r2q = 0.02458;
  Real l2q = 4.129e-4;
//   Real rSource = 1.0;

  Real inertia = 2.25e4;
  Real shaftTorque = 3.1831e5;
  Real vExcitation = 7.0829; // 92.95;

  Real voltageMagnitude = voltageMagnitudeLL * sqrt(2.0) / sqrt(3.0);
  Complex voltageL1 = Complex(voltageMagnitude, 0.0);
  Complex voltageL2 = Complex(voltageMagnitude * cos(-2 * M_PI / 3),
                              voltageMagnitude * sin(-2 * M_PI / 3));
  Complex voltageL3 = Complex(voltageMagnitude * cos(2 * M_PI / 3),
                              voltageMagnitude * sin(2 * M_PI / 3));

  // Nodes
  auto n1 = SimNode::make("n1");
  auto n2 = SimNode::make("n2");
  auto n3 = SimNode::make("n3");
  auto n4 = SimNode::make("n4");
  auto n5 = SimNode::make("n5");
  auto n6 = SimNode::make("n6");
  auto n7 = SimNode::make("n7");
  auto n8 = SimNode::make("n8");
  auto n9 = SimNode::make("n9");
  auto n10 = SimNode::make("n10");
  auto n11 = SimNode::make("n11");
  auto n12 = SimNode::make("n12");
  auto n13 = SimNode::make("n13");
  auto n14 = SimNode::make("n14");
  auto n15 = SimNode::make("n15");
  auto n16 = SimNode::make("n16");
  auto n17 = SimNode::make("n17");
  auto n18 = SimNode::make("n18");
  auto n19 = SimNode::make("n19");
  auto n20 = SimNode::make("n20");
  auto n21 = SimNode::make("n21");
  auto n22 = SimNode::make("n22");
  auto n23 = SimNode::make("n23");
  auto n24 = SimNode::make("n24");
  auto n25 = SimNode::make("n25");
  auto n26 = SimNode::make("n26");
  auto n27 = SimNode::make("n27");
  auto n28 = SimNode::make("n28");
  auto n29 = SimNode::make("n29");
  auto n30 = SimNode::make("n30");
  auto n31 = SimNode::make("n31");
  auto n32 = SimNode::make("n32");
  auto n33 = SimNode::make("n33");
  auto n34 = SimNode::make("n34");
//   auto n35 = SimNode::make("n35");
//   auto n36 = SimNode::make("n36");
//   auto n37 = SimNode::make("n37");

  // Components
  // Infeed
  auto v1 = VoltageSource::make("v1", Logger::Level::debug);
  v1->setParameters(voltageL1, frequency);
  v1->connect({SimNode::GND, n32});
  auto v2 = VoltageSource::make("v2", Logger::Level::debug);
  v2->setParameters(voltageL2, frequency);
  v2->connect({SimNode::GND, n33});
  auto v3 = VoltageSource::make("v3", Logger::Level::debug);
  v3->setParameters(voltageL3, frequency);
  v3->connect({SimNode::GND, n34});

//   auto rSourceComp1 = Resistor::make("rSourceComp1", Logger::Level::debug);
//   rSourceComp1->setParameters(rSource);
//   rSourceComp1->connect({n32, n35});
//   auto rSourceComp2 = Resistor::make("rSourceComp2", Logger::Level::debug);
//   rSourceComp2->setParameters(rSource);
//   rSourceComp2->connect({n33, n36});
//   auto rSourceComp3 = Resistor::make("rSourceComp3", Logger::Level::debug);
//   rSourceComp3->setParameters(rSource);
//   rSourceComp3->connect({n34, n37});
  // Park transformer
  auto parkTrafo =
      ParkTransformer::make("ParkTransformer", Logger::Level::debug);
  parkTrafo->setIsOmegaConstant(false);
  parkTrafo->setOmegaReferenceNode(n30);
  parkTrafo->setInitialValues(omega, 0.0);
  parkTrafo->connect({n32, n33, n34, n1, n2, n3});

    // d-axis components
    auto rSd = Resistor::make("rSd", Logger::Level::debug);
    rSd->setParameters(rS);
    rSd->connect({n1, n4});

    auto tlvs1 = TimeLaggingVoltageSource::make("tlvs1", Logger::Level::debug);
    tlvs1->setVoltageReferenceNodes(n4, n13);
    tlvs1->connect({n7, n4});

    auto rec1 = ElectroMechanicalConverter::make("rec1", Logger::Level::debug);
    rec1->setInitialFlux(0.0);
    rec1->setIsNegative(false);
    rec1->setVoltageReferenceNode(n18);
    rec1->connect({n7, n9, n30, SimNode::GND});

    auto ccvs1 = VoltageSpeedTerm::make("ccvs1", Logger::Level::debug);
    ccvs1->setInitialOmega(0.0);
    ccvs1->setIsConstantSpeed(false);
    ccvs1->setOmegaReferenceNode(n30);
    ccvs1->setInductance(lMq);
    ccvs1->setIsNegative(false);
    ccvs1->connect({n24, SimNode::GND, n9, n11});

    auto ccvs2 = VoltageSpeedTerm::make("ccvs2", Logger::Level::debug);
    ccvs2->setInitialOmega(0.0);
    ccvs2->setIsConstantSpeed(false);
    ccvs2->setOmegaReferenceNode(n30);
    ccvs2->setInductance(lS);
    ccvs2->setIsNegative(false);
    ccvs2->connect({n16, n18, n11, n13});

    auto lSd = Inductor::make("lSd", Logger::Level::debug);
    lSd->setParameters(lS);
    lSd->connect({n17, n19});

    auto lMdComponent = Inductor::make("lMdComponent", Logger::Level::debug);
    lMdComponent->setParameters(lMd);
    lMdComponent->connect({n19, n21});

    auto lfdComponent = Inductor::make("lfdComponent", Logger::Level::debug);
    lfdComponent->setParameters(lFd);
    lfdComponent->connect({n19, n25});

    auto rfdComponent = Resistor::make("rfdComponent", Logger::Level::debug);
    rfdComponent->setParameters(rFd);
    rfdComponent->connect({n25, n27});

    auto vFd = VoltageSource::make("vF", Logger::Level::debug);
    vFd->setParameters(vExcitation, 0.0);
    vFd->connect({n27, SimNode::GND});

    auto l1dComponent = Inductor::make("l1dComponent", Logger::Level::debug);
    l1dComponent->setParameters(l1d);
    l1dComponent->connect({n19, n28});

    auto r1dComponent = Resistor::make("r1dComponent", Logger::Level::debug);
    r1dComponent->setParameters(r1d);
    r1dComponent->connect({n28, SimNode::GND});

    // q-axis components
    auto rSq = Resistor::make("rSq", Logger::Level::debug);
    rSq->setParameters(rS);
    rSq->connect({n2, n5});

    auto tlvs2 = TimeLaggingVoltageSource::make("tlvs2", Logger::Level::debug);
    tlvs2->setVoltageReferenceNodes(n14, n5);
    tlvs2->connect({n5, n8});

    auto rec2 = ElectroMechanicalConverter::make("rec2", Logger::Level::debug);
    rec2->setInitialFlux(0.0);
    rec2->setIsNegative(false);
    rec2->setVoltageReferenceNode(n17);
    rec2->connect({n10, n8, n30, SimNode::GND});

    auto ccvs3 = VoltageSpeedTerm::make("ccvs3", Logger::Level::debug);
    ccvs3->setInitialOmega(0.0);
    ccvs3->setIsConstantSpeed(false);
    ccvs3->setOmegaReferenceNode(n30);
    ccvs3->setInductance(lMd);
    ccvs3->setIsNegative(false);
    ccvs3->connect({n23, SimNode::GND, n12, n10});

    auto ccvs4 = VoltageSpeedTerm::make("ccvs4", Logger::Level::debug);
    ccvs4->setInitialOmega(0.0);
    ccvs4->setIsConstantSpeed(false);
    ccvs4->setOmegaReferenceNode(n30);
    ccvs4->setInductance(lS);
    ccvs4->setIsNegative(false);
    ccvs4->connect({n15, n17, n14, n12});

    auto lSq = Inductor::make("lSq", Logger::Level::debug);
    lSq->setParameters(lS);
    lSq->connect({n18, n20});

    auto lMqComponent = Inductor::make("lMqComponent", Logger::Level::debug);
    lMqComponent->setParameters(lMq);
    lMqComponent->connect({n20, n22});

    auto l1qComponent = Inductor::make("l1qComponent", Logger::Level::debug);
    l1qComponent->setParameters(l1q);
    l1qComponent->connect({n20, n26});

    auto r1qComponent = Resistor::make("r1qComponent", Logger::Level::debug);
    r1qComponent->setParameters(r1q);
    r1qComponent->connect({n26, SimNode::GND});

    auto l2qComponent = Inductor::make("l2qComponent", Logger::Level::debug);
    l2qComponent->setParameters(l2q);
    l2qComponent->connect({n20, n29});

    auto r2qComponent = Resistor::make("r2qComponent", Logger::Level::debug);
    r2qComponent->setParameters(r2q);
    r2qComponent->connect({n29, SimNode::GND});

    // 0-axis components
    auto r0 = Resistor::make("r0", Logger::Level::debug);
    r0->setParameters(rS);
    r0->connect({n3, n6});
    auto l0 = Inductor::make("l0", Logger::Level::debug);
    l0->setParameters(lS);
    l0->connect({n6, SimNode::GND});

    // Mechanical components
    auto inertiaMoment =
        InertiaMoment::make("inertiaMoment", Logger::Level::debug);
    inertiaMoment->setParameters(inertia);
    inertiaMoment->connect({n31, SimNode::GND});

    auto constantOmegaSource =
        VoltageSource::make("constantOmegaSource", Logger::Level::debug);
    constantOmegaSource->setParameters(
        omega, 0.0); // Set the constant angular frequency source
    constantOmegaSource->connect({SimNode::GND, n31});

    auto ccts1 =
        CurrentControlledTorqueSource::make("ccts1", Logger::Level::debug);
    ccts1->setInductance(lMq);
    ccts1->setInductor(lSq);
    ccts1->connect({n21, n23, SimNode::GND, n30});

    auto ccts2 =
        CurrentControlledTorqueSource::make("ccts2", Logger::Level::debug);
    ccts2->setInductance(lS);
    ccts2->setInductor(lSq);
    ccts2->connect({n13, n15, SimNode::GND, n30});

    auto ccts3 =
        CurrentControlledTorqueSource::make("ccts3", Logger::Level::debug);
    ccts3->setInductance(lMd);
    ccts3->setInductor(lSd);
    ccts3->connect({n22, n24, n30, SimNode::GND});

    auto ccts4 =
        CurrentControlledTorqueSource::make("ccts4", Logger::Level::debug);
    ccts4->setInductance(lS);
    ccts4->setInductor(lSd);
    ccts4->connect({n14, n16, n30, SimNode::GND});

    auto tlts = TimeLaggingTorqueSource::make("tlts", Logger::Level::debug);
    tlts->connect({n30, n31, n30, SimNode::GND});

    auto tMech = CurrentSource::make("tMech", Logger::Level::debug);
    tMech->setParameters(-shaftTorque, 0.0);
    tMech->connect({n31, SimNode::GND});

    // Define system topology
    SystemTopology system(
        50,
        SystemNodeList{n1,  n2,  n3,  n4,  n5,  n6,  n7,  n8,  n9,  n10, n11, n12,
                       n13, n14, n15, n16, n17, n18, n19, n20, n21, n22, n23, n24,
                       n25, n26, n27, n28, n29, n30, n31, n32, n33, n34,
                        // n35, n36, n37
                       },
        SystemComponentList{
            v1,        v2,        v3,        
            // rSourceComp1, rSourceComp2, rSourceComp3, 
            parkTrafo, rSd,       tlvs1,        rec1,
            ccvs1,     ccvs2,     lSd,       lMdComponent, lfdComponent,
            rfdComponent, vFd,      l1dComponent, r1dComponent, rSq,
            tlvs2,     rec2,      ccvs3,     ccvs4,       lSq,
            lMqComponent, l1qComponent, r1qComponent, l2qComponent, r2qComponent,
            r0,        l0,        inertiaMoment, ccts1,       ccts2,
            ccts3,     ccts4,     tlts,      tMech,
            constantOmegaSource 
            });

    // Define simulation scenario
    String simName = "synchronousGeneratorTest";
    // Logger
    auto logger = DataLogger::make(simName);
    logger->logAttribute("i_A", v1->attribute("i_intf"));
    logger->logAttribute("omega", n31->attribute("v"));
    logger->logAttribute("torqueDiff", inertiaMoment->attribute("i_intf"));
    logger->logAttribute("torque_EM_old", tlts->attribute("i"));
    logger->logAttribute("torque_Mech", tMech->attribute("i_intf"));
    logger->logAttribute("iSd", lSd->attribute("i_intf"));
    logger->logAttribute("vSd", lSd->attribute("v_intf"));
    logger->logAttribute("iSq", lSq->attribute("i_intf"));
    logger->logAttribute("vSq", lSq->attribute("v_intf"));
    logger->logAttribute("iFd", lfdComponent->attribute("i_intf"));
    logger->logAttribute("vFd", lfdComponent->attribute("v_intf"));
    logger->logAttribute("iMd", lMdComponent->attribute("i_intf"));
    logger->logAttribute("vMd", lMdComponent->attribute("v_intf"));
    logger->logAttribute("iMq", lMqComponent->attribute("i_intf"));
    logger->logAttribute("vMq", lMqComponent->attribute("v_intf"));

    Simulation sim(simName);
    sim.setSystem(system);
    sim.setTimeStep(timeStep);
    sim.doSystemMatrixRecomputation(true);
    sim.setFinalTime(finalTime);
    sim.setDomain(Domain::EMT);
    sim.addLogger(logger);
    sim.doEigenvalueExtraction(doEigenvalueExtraction);
    inertiaMoment->setIntfVoltage(-omega);
    sim.run();
}

int main(int argc, char *argv[]) {
  //   motorStartingTest(5e-5, 3.0, false);
  synchronousGeneratorTest(1e-4, 8, false);
  // timeLaggingTorqueSourceTest();
  // timeLaggingVoltageSourceTest();
  //  speedVoltageTermTest(1e-4, 1, true);
  //   syncronousGeneratorOldTest(1e-3, 70.0, false);
  // electromechanicalConverterTest(1e-4, 1.0, true);

  return 0;
}
