#include <DPsim.h>
#include <iostream>
#include <string>
#include <vector>

using namespace DPsim;
using namespace CPS::EMT;
using namespace CPS::EMT::Ph1;

void syncronousGeneratorTest(Real timeStep, Real finalTime,
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

void motorStartingTestRotorReference(Real timeStep, Real finalTime,
                                     bool doEigenvalueExtraction) {
  Real frequency = 50;
  Real omega = 2 * M_PI * frequency;
  Real voltageMagnitudeLL = 10e3;
  //   Real infeedResistance = 3.4e-3;
  //   Real infeedReactance = 9.4e-3;
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
  //   Real infeedInductance = infeedReactance / omega;
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
  //   auto rInfeed1 = Resistor::make("rInfeed", Logger::Level::debug);
  //   rInfeed1->setParameters(infeedResistance);
  //   rInfeed1->connect({n4, n1});
  //   auto rInfeed2 = Resistor::make("rInfeed2", Logger::Level::debug);
  //   rInfeed2->setParameters(infeedResistance);
  //   rInfeed2->connect({n5, n2});
  //   auto rInfeed3 = Resistor::make("rInfeed3", Logger::Level::debug);
  //   rInfeed3->setParameters(infeedResistance);
  //   rInfeed3->connect({n6, n3});
  //   auto lInfeed1 = Inductor::make("lInfeed", Logger::Level::debug);
  //   lInfeed1->setParameters(infeedInductance);
  //   lInfeed1->connect({n4, n7});
  //   auto lInfeed2 = Inductor::make("lInfeed2", Logger::Level::debug);
  //   lInfeed2->setParameters(infeedInductance);
  //   lInfeed2->connect({n5, n8});
  //   auto lInfeed3 = Inductor::make("lInfeed3", Logger::Level::debug);
  //   lInfeed3->setParameters(infeedInductance);
  //   lInfeed3->connect({n6, n9});

  // Park transformer
  auto parkTrafo =
      ParkTransformer::make("ParkTransformer", Logger::Level::debug);
  parkTrafo->setIsOmegaConstant(false);
  parkTrafo->setInitialValues(0.0, 0.0);
  //   parkTrafo->setInitialValues(omega, 0.0);
  //   parkTrafo->setIsOmegaConstant(true);
  parkTrafo->setOmegaReferenceNode(n22);
  //   parkTrafo->connect({n7, n8, n9, n10, n11, n12});
  parkTrafo->connect({n1, n2, n3, n10, n11, n12});

  // d-axis components
  auto rSd = Resistor::make("rSd", Logger::Level::debug);
  rSd->setParameters(rS);
  rSd->connect({n10, n13});
  auto lSd = Inductor::make("lSd", Logger::Level::debug);
  lSd->setParameters(lS);
  lSd->connect({n16, n18});
  auto lMd = Inductor::make("lMd", Logger::Level::debug);
  lMd->setParameters(lM);
  lMd->connect({n18, SimNode::GND});
  auto lRd = Inductor::make("lRd", Logger::Level::debug);
  lRd->setParameters(lR);
  lRd->connect({n18, n20});
  auto rRd = Resistor::make("rRd", Logger::Level::debug);
  rRd->setParameters(rR);
  rRd->connect({n20, SimNode::GND});
  auto electroMechanicalConverter_d = ElectroMechanicalConverter::make(
      "electroMechanicalConverter_d", Logger::Level::debug);
  electroMechanicalConverter_d->setVoltageReferenceNode(n17);
  electroMechanicalConverter_d->setIsNegative(true);
  electroMechanicalConverter_d->connect({n13, n16, n22, SimNode::GND});

  // q-axis components
  auto rSq = Resistor::make("rSq", Logger::Level::debug);
  rSq->setParameters(rS);
  rSq->connect({n11, n14});
  auto lSq = Inductor::make("lSq", Logger::Level::debug);
  lSq->setParameters(lS);
  lSq->connect({n17, n19});
  auto lMq = Inductor::make("lMq", Logger::Level::debug);
  lMq->setParameters(lM);
  lMq->connect({n19, SimNode::GND});
  auto lRq = Inductor::make("lRq", Logger::Level::debug);
  lRq->setParameters(lR);
  lRq->connect({n19, n21});
  auto rRq = Resistor::make("rRq", Logger::Level::debug);
  rRq->setParameters(rR);
  rRq->connect({n21, SimNode::GND});
  auto electroMechanicalConverter_q = ElectroMechanicalConverter::make(
      "electroMechanicalConverter_q", Logger::Level::debug);
  electroMechanicalConverter_q->setVoltageReferenceNode(n16);
  electroMechanicalConverter_q->connect({n14, n17, n22, SimNode::GND});

  // 0-axis components
  auto r0 = Resistor::make("r0", Logger::Level::debug);
  r0->setParameters(rS);
  r0->connect({n12, n15});
  auto l0 = Inductor::make("l0", Logger::Level::debug);
  l0->setParameters(lS);
  l0->connect({n15, SimNode::GND});

  // Mechanical components
  auto inertiaMoment =
      InertiaMoment::make("inertiaMoment", Logger::Level::debug);
  inertiaMoment->setParameters(inertia);
  inertiaMoment->connect({n22, SimNode::GND});

  // Define system topology
  SystemTopology system(
      50,
      SystemNodeList{n1, n2, n3,
                     // n4, n5, n6, n7, n8, n9,
                     n10, n11, n12, n13, n14, n15, n16, n17, n18, n19, n20, n21,
                     n22},
      SystemComponentList{v1, v2, v3,
                          // rInfeed1,
                          // rInfeed2,
                          // rInfeed3,
                          // lInfeed1,
                          // lInfeed2,
                          // lInfeed3,
                          parkTrafo, rSd, lSd, lMd, lRd, rRd,
                          electroMechanicalConverter_d, rSq, lSq, lMq, lRq, rRq,
                          electroMechanicalConverter_q, r0, l0, inertiaMoment});
  // Define simulation scenario
  String simName = "motorStartingTest";
  // Logger
  auto logger = DataLogger::make(simName);
  logger->logAttribute("omega", n22->attribute("v"));
  logger->logAttribute("Flux_d",
                       electroMechanicalConverter_d->attribute("flux"));
  logger->logAttribute("Flux_q",
                       electroMechanicalConverter_q->attribute("flux"));
  // logger->logAttribute("omega_Park", parkTrafo->attribute("omega"));
  // logger->logAttribute("V_A", n1->attribute("v"));
  // logger->logAttribute("V_B", n2->attribute("v"));
  // logger->logAttribute("V_C", n3->attribute("v"));
  // logger->logAttribute("I_A", rInfeed1->attribute("i_intf"));
  // logger->logAttribute("I_B", rInfeed2->attribute("i_intf"));
  // logger->logAttribute("I_C", rInfeed3->attribute("i_intf"));
  // logger->logAttribute("V_D", n10->attribute("v"));
  // logger->logAttribute("V_Q", n11->attribute("v"));
  // logger->logAttribute("V_0", n12->attribute("v"));
  logger->logAttribute("I_sd", rSd->attribute("i_intf"));
  logger->logAttribute("I_sq", rSq->attribute("i_intf"));
  logger->logAttribute("Torque_d",
                       electroMechanicalConverter_d->attribute("torque"));
  logger->logAttribute("Torque_q",
                       electroMechanicalConverter_q->attribute("torque"));

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

void motorStartingTestRotorReferenceSpeedVoltageTerms(
    Real timeStep, Real finalTime, bool doEigenvalueExtraction) {
  Real frequency = 50;
  Real omega = 2 * M_PI * frequency;
  Real voltageMagnitudeLL = 10e3;
  // Real infeedResistance = 3.4e-3;
  // Real infeedReactance = 9.4e-3;
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
  // Real infeedInductance = infeedReactance / omega;
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

  // Components
  // Infeed
  auto v1 = VoltageSource::make("v1", Logger::Level::debug);
  v1->setParameters(voltageL1, frequency);
  //   v1->connect({SimNode::GND, n25});
  v1->connect({SimNode::GND, n1});
  auto v2 = VoltageSource::make("v2", Logger::Level::debug);
  v2->setParameters(voltageL2, frequency);
  //   v2->connect({SimNode::GND, n26});
  v2->connect({SimNode::GND, n2});
  auto v3 = VoltageSource::make("v3", Logger::Level::debug);
  v3->setParameters(voltageL3, frequency);
  //   v3->connect({SimNode::GND, n27});
  v3->connect({SimNode::GND, n3});

  //   auto rInfeed1 = Resistor::make("rInfeed", Logger::Level::debug);
  //   rInfeed1->setParameters(infeedResistance);
  //   rInfeed1->connect({n28, n25});
  //   auto rInfeed2 = Resistor::make("rInfeed2", Logger::Level::debug);
  //   rInfeed2->setParameters(infeedResistance);
  //   rInfeed2->connect({n29, n26});
  //   auto rInfeed3 = Resistor::make("rInfeed3", Logger::Level::debug);
  //   rInfeed3->setParameters(infeedResistance);
  //   rInfeed3->connect({n30, n27});
  //   auto lInfeed1 = Inductor::make("lInfeed", Logger::Level::debug);
  //   lInfeed1->setParameters(infeedInductance);
  //   lInfeed1->connect({n28, n1});
  //   auto lInfeed2 = Inductor::make("lInfeed2", Logger::Level::debug);
  //   lInfeed2->setParameters(infeedInductance);
  //   lInfeed2->connect({n29, n2});
  //   auto lInfeed3 = Inductor::make("lInfeed3", Logger::Level::debug);
  //   lInfeed3->setParameters(infeedInductance);
  //   lInfeed3->connect({n30, n3});

  // Park transformer
  auto parkTrafo =
      ParkTransformer::make("ParkTransformer", Logger::Level::debug);
  parkTrafo->setIsOmegaConstant(false);
  parkTrafo->setInitialValues(0.0, 0.0);
  parkTrafo->setOmegaReferenceNode(n24);
  parkTrafo->connect({n1, n2, n3, n4, n5, n6});

  // d-axis components
  auto rSd = Resistor::make("rSd", Logger::Level::debug);
  rSd->setParameters(rS);
  rSd->connect({n4, n7});
  auto lSd = Inductor::make("lSd", Logger::Level::debug);
  lSd->setParameters(lS);
  lSd->connect({n16, n18});
  auto lMd = Inductor::make("lMd", Logger::Level::debug);
  lMd->setParameters(lM);
  lMd->connect({n18, n20});
  auto lRd = Inductor::make("lRd", Logger::Level::debug);
  lRd->setParameters(lR);
  lRd->connect({n18, n22});
  auto rRd = Resistor::make("rRd", Logger::Level::debug);
  rRd->setParameters(rR);
  rRd->connect({n22, SimNode::GND});
  auto voltageSpeedTermDFluxS =
      VoltageSpeedTerm::make("voltageSpeedTermDFluxS", Logger::Level::debug);
  voltageSpeedTermDFluxS->setIsNegative(true);
  voltageSpeedTermDFluxS->setInitialOmega(0.0);
  voltageSpeedTermDFluxS->setIsConstantSpeed(false);
  voltageSpeedTermDFluxS->setInductance(lS);
  voltageSpeedTermDFluxS->setOmegaReferenceNode(n24);
  voltageSpeedTermDFluxS->connect({n13, n15, n7, n10});
  auto voltageSpeedTermDFluxM =
      VoltageSpeedTerm::make("voltageSpeedTermDFluxM", Logger::Level::debug);
  voltageSpeedTermDFluxM->setIsNegative(true);
  voltageSpeedTermDFluxM->setInitialOmega(0.0);
  voltageSpeedTermDFluxM->setIsConstantSpeed(false);
  voltageSpeedTermDFluxM->setInductance(lM);
  voltageSpeedTermDFluxM->setOmegaReferenceNode(n24);
  voltageSpeedTermDFluxM->connect({n21, SimNode::GND, n10, n12});

  // q-axis components
  auto rSq = Resistor::make("rSq", Logger::Level::debug);
  rSq->setParameters(rS);
  rSq->connect({n5, n8});
  auto lSq = Inductor::make("lSq", Logger::Level::debug);
  lSq->setParameters(lS);
  lSq->connect({n17, n19});
  auto lMq = Inductor::make("lMq", Logger::Level::debug);
  lMq->setParameters(lM);
  lMq->connect({n19, n21});
  auto lRq = Inductor::make("lRq", Logger::Level::debug);
  lRq->setParameters(lR);
  lRq->connect({n19, n23});
  auto rRq = Resistor::make("rRq", Logger::Level::debug);
  rRq->setParameters(rR);
  rRq->connect({n23, SimNode::GND});
  auto voltageSpeedTermQFluxS =
      VoltageSpeedTerm::make("voltageSpeedTermQFluxS", Logger::Level::debug);
  voltageSpeedTermQFluxS->setInitialOmega(0.0);
  voltageSpeedTermQFluxS->setIsConstantSpeed(false);
  voltageSpeedTermQFluxS->setInductance(lS);
  voltageSpeedTermQFluxS->setOmegaReferenceNode(n24);
  voltageSpeedTermQFluxS->connect({n12, n14, n8, n11});

  auto voltageSpeedTermQFluxM =
      VoltageSpeedTerm::make("voltageSpeedTermQFluxM", Logger::Level::debug);
  voltageSpeedTermQFluxM->setInitialOmega(0.0);
  voltageSpeedTermQFluxM->setIsConstantSpeed(false);
  voltageSpeedTermQFluxM->setInductance(lM);
  voltageSpeedTermQFluxM->setOmegaReferenceNode(n24);
  voltageSpeedTermQFluxM->connect({n20, SimNode::GND, n11, n13});

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
  inertiaMoment->connect({n24, SimNode::GND});

  auto currentDControlledTorqueSource = CurrentControlledTorqueSource::make(
      "currentDControlledTorqueSource", Logger::Level::debug);
  currentDControlledTorqueSource->setInitialFlux(0.0);
  currentDControlledTorqueSource->setCoefficient(-1.0);
  currentDControlledTorqueSource->setVoltageReferenceNode(n17);
  currentDControlledTorqueSource->connect({n14, n16, SimNode::GND, n24});

  auto currentQControlledTorqueSource = CurrentControlledTorqueSource::make(
      "currentQControlledTorqueSource", Logger::Level::debug);
    currentQControlledTorqueSource->setInitialFlux(0.0);
  currentQControlledTorqueSource->setCoefficient(1.0);
    currentQControlledTorqueSource->setVoltageReferenceNode(n16);
  currentQControlledTorqueSource->connect({n15, n17, SimNode::GND, n24});

  // Define system topology
  SystemTopology system(50,
                        SystemNodeList{n1,  n2,  n3,  n4,  n5,  n6,  n7,  n8,
                                       n9,  n10, n11, n12, n13, n14, n15, n16,
                                       n17, n18, n19, n20, n21, n22, n23, n24,
                                    //    n25, n26, n27, n28, n29, n30
                                       },
                        SystemComponentList{
                            v1,
                            v2,
                            v3,
                            // rInfeed1,
                            // rInfeed2,
                            // rInfeed3,
                            // lInfeed1,
                            // lInfeed2,
                            // lInfeed3,
                            parkTrafo,
                            rSd,
                            lSd,
                            lMd,
                            lRd,
                            rRd,
                            voltageSpeedTermDFluxS,
                            voltageSpeedTermDFluxM,
                            rSq,
                            lSq,
                            lMq,
                            lRq,
                            rRq,
                            voltageSpeedTermQFluxS,
                            voltageSpeedTermQFluxM,
                            r0,
                            l0,
                            inertiaMoment,
                            currentDControlledTorqueSource,
                            currentQControlledTorqueSource,
                        });
  // Define simulation scenario
  String simName = "motorStartingTest";
  // Logger
  auto logger = DataLogger::make(simName);
  logger->logAttribute("omega", n24->attribute("v"));
  logger->logAttribute("omega_park", parkTrafo->attribute("omega"));
  logger->logAttribute("I_sd", rSd->attribute("i_intf"));
  logger->logAttribute("I_sq", rSq->attribute("i_intf"));
  logger->logAttribute("I_s0", r0->attribute("i_intf"));
//   logger->logAttribute("I_A", rInfeed1->attribute("i_intf"));
//   logger->logAttribute("I_B", rInfeed2->attribute("i_intf"));
//   logger->logAttribute("I_C", rInfeed3->attribute("i_intf"));

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

void motorStartingTestSynchReferenceSpeedVoltageTerms(
    Real timeStep, Real finalTime, bool doEigenvalueExtraction) {
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
  auto lSd = Inductor::make("lSd", Logger::Level::debug);
  lSd->setParameters(lS);
  lSd->connect({n16, n18});
  auto lMd = Inductor::make("lMd", Logger::Level::debug);
  lMd->setParameters(lM);
  lMd->connect({n18, n20});
  auto lRd = Inductor::make("lRd", Logger::Level::debug);
  lRd->setParameters(lR);
  lRd->connect({n18, n24});
  auto rRd = Resistor::make("rRd", Logger::Level::debug);
  rRd->setParameters(rR);
  rRd->connect({n30, SimNode::GND});

  auto voltageSpeedTermDFluxS =
      VoltageSpeedTerm::make("voltageSpeedTermDFluxS", Logger::Level::debug);
  voltageSpeedTermDFluxS->setIsNegative(true);
  voltageSpeedTermDFluxS->setInitialOmega(omega);
  voltageSpeedTermDFluxS->setIsConstantSpeed(true);
  voltageSpeedTermDFluxS->setInductance(lS);
  voltageSpeedTermDFluxS->connect({n13, n15, n7, n10});

  auto voltageSpeedTermDFluxM =
      VoltageSpeedTerm::make("voltageSpeedTermDFluxM", Logger::Level::debug);
  voltageSpeedTermDFluxM->setIsNegative(true);
  voltageSpeedTermDFluxM->setInitialOmega(omega);
  voltageSpeedTermDFluxM->setIsConstantSpeed(true);
  voltageSpeedTermDFluxM->setInductance(lM);
  voltageSpeedTermDFluxM->connect({n21, n23, n10, n12});

auto voltageSpeedTermDFluxR =
      VoltageSpeedTerm::make("voltageSpeedTermDFluxR", Logger::Level::debug);
    voltageSpeedTermDFluxR->setIsNegative(true);
    voltageSpeedTermDFluxR->setInitialOmega(0.0);
    voltageSpeedTermDFluxR->setIsConstantSpeed(false);
    voltageSpeedTermDFluxR->setInductance(lR);
    voltageSpeedTermDFluxR->setOmegaReferenceNode(n32);
    voltageSpeedTermDFluxR->setOmegaOffset(-omega);
    voltageSpeedTermDFluxR->connect({n31, n29, n24, n26});

    auto voltageSpeedTermDFluxRM =
        VoltageSpeedTerm::make("voltageSpeedTermDFluxRM", Logger::Level::debug);
    voltageSpeedTermDFluxRM->setIsNegative(true);
    voltageSpeedTermDFluxRM->setInitialOmega(0.0);
    voltageSpeedTermDFluxRM->setIsConstantSpeed(false);
    voltageSpeedTermDFluxRM->setInductance(lM);
    voltageSpeedTermDFluxRM->setOmegaReferenceNode(n32);
    voltageSpeedTermDFluxRM->setOmegaOffset(-omega);
    voltageSpeedTermDFluxRM->connect({n23, SimNode::GND, n26, n28});

    // q-axis components
    auto rSq = Resistor::make("rSq", Logger::Level::debug);
    rSq->setParameters(rS);
    rSq->connect({n5, n8});
    auto lSq = Inductor::make("lSq", Logger::Level::debug);
    lSq->setParameters(lS);
    lSq->connect({n17, n19});
    auto lMq = Inductor::make("lMq", Logger::Level::debug);
    lMq->setParameters(lM);
    lMq->connect({n19, n21});
    auto lRq = Inductor::make("lRq", Logger::Level::debug);
    lRq->setParameters(lR);
    lRq->connect({n19, n25});
    auto rRq = Resistor::make("rRq", Logger::Level::debug);
    rRq->setParameters(rR);
    rRq->connect({n31, SimNode::GND});
    
    auto voltageSpeedTermQFluxS =
        VoltageSpeedTerm::make("voltageSpeedTermQFluxS", Logger::Level::debug);
    voltageSpeedTermQFluxS->setInitialOmega(omega);
    voltageSpeedTermQFluxS->setIsConstantSpeed(true);
    voltageSpeedTermQFluxS->setInductance(lS);
    voltageSpeedTermQFluxS->connect({n12, n14, n8, n11});

    auto voltageSpeedTermQFluxM =
        VoltageSpeedTerm::make("voltageSpeedTermQFluxM", Logger::Level::debug);
    voltageSpeedTermQFluxM->setInitialOmega(omega);
    voltageSpeedTermQFluxM->setIsConstantSpeed(true);
    voltageSpeedTermQFluxM->setInductance(lM);
    voltageSpeedTermQFluxM->connect({n20, n22, n11, n13});

    auto voltageSpeedTermQFluxR =
        VoltageSpeedTerm::make("voltageSpeedTermQFluxR", Logger::Level::debug);
    voltageSpeedTermQFluxR->setInitialOmega(0.0);
    voltageSpeedTermQFluxR->setIsConstantSpeed(false);
    voltageSpeedTermQFluxR->setInductance(lR);
    voltageSpeedTermQFluxR->setOmegaReferenceNode(n32);
    voltageSpeedTermQFluxR->setOmegaOffset(-omega);
    voltageSpeedTermQFluxR->connect({n30, n28, n25, n27});

    auto voltageSpeedTermQFluxRM =
        VoltageSpeedTerm::make("voltageSpeedTermQFluxRM", Logger::Level::debug);
    voltageSpeedTermQFluxRM->setInitialOmega(0.0);
    voltageSpeedTermQFluxRM->setIsConstantSpeed(false);
    voltageSpeedTermQFluxRM->setInductance(lM);
    voltageSpeedTermQFluxRM->setOmegaReferenceNode(n32);
    voltageSpeedTermQFluxRM->setOmegaOffset(-omega);
    voltageSpeedTermQFluxRM->connect({n22, SimNode::GND, n27, n29});

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
    inertiaMoment->connect({n32, SimNode::GND});

    auto currentDControlledTorqueSource = CurrentControlledTorqueSource::make(
        "currentDControlledTorqueSource", Logger::Level::debug);
    currentDControlledTorqueSource->setInitialFlux(0.0);
    currentDControlledTorqueSource->setCoefficient(-1.0);
    currentDControlledTorqueSource->setVoltageReferenceNode(n17);
    currentDControlledTorqueSource->connect({n14, n16, SimNode::GND, n32});

    auto currentQControlledTorqueSource = CurrentControlledTorqueSource::make(
        "currentQControlledTorqueSource", Logger::Level::debug);
    currentQControlledTorqueSource->setInitialFlux(0.0);
    currentQControlledTorqueSource->setCoefficient(1.0);
    currentQControlledTorqueSource->setVoltageReferenceNode(n16);
    currentQControlledTorqueSource->connect({n15, n17, SimNode::GND, n32});

    // Define system topology
    SystemTopology system(
        50,
        SystemNodeList{n1,  n2,  n3,  n4,  n5,  n6,  n7,  n8,  n9,  n10, n11,
                       n12, n13, n14, n15, n16, n17, n18, n19, n20, n21, n22,
                       n23, n24, n25, n26, n27, n28, n29, n30, n31, n32},
        SystemComponentList{
            v1,
            v2,
            v3,
            parkTrafo,
            rSd,
            lSd,
            lMd,
            lRd,
            rRd,
            voltageSpeedTermDFluxS,
            voltageSpeedTermDFluxM,
            voltageSpeedTermDFluxR,
            voltageSpeedTermDFluxRM,
            rSq,
            lSq,
            lMq,
            lRq,
            rRq,
            voltageSpeedTermQFluxS,
            voltageSpeedTermQFluxM,
            voltageSpeedTermQFluxR,
            voltageSpeedTermQFluxRM,
            r0,
            l0,
            inertiaMoment,
            currentDControlledTorqueSource,
            currentQControlledTorqueSource,
        });
    // Define simulation scenario
    String simName = "motorStartingTest";
    // Logger
    auto logger = DataLogger::make(simName);
    logger->logAttribute("omega", n32->attribute("v"));
    logger->logAttribute("omega_park", parkTrafo->attribute("omega"));
    logger->logAttribute("I_sd", rSd->attribute("i_intf"));
    logger->logAttribute("I_sq", rSq->attribute("i_intf"));
    logger->logAttribute("I_s0", r0->attribute("i_intf"));

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

void motorStartingTestSynchReferenceSpeedVoltageTermsEMConverter(
    Real timeStep, Real finalTime, bool doEigenvalueExtraction) {
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
//   parkTrafo->setInitialValues(omega, 0.0);
  parkTrafo->setInitialValues(omega, PI / 4);
  parkTrafo->connect({n1, n2, n3, n4, n5, n6});

  // d-axis components
  auto rSd = Resistor::make("rSd", Logger::Level::debug);
  rSd->setParameters(rS);
  rSd->connect({n4, n7});
  auto lSd = Inductor::make("lSd", Logger::Level::debug);
  lSd->setParameters(lS);
  lSd->connect({n14, n16});
  auto lMd = Inductor::make("lMd", Logger::Level::debug);
  lMd->setParameters(lM);
  lMd->connect({n16, n18});
  auto lRd = Inductor::make("lRd", Logger::Level::debug);
  lRd->setParameters(lR);
  lRd->connect({n16, n22});
  auto rRd = Resistor::make("rRd", Logger::Level::debug);
  rRd->setParameters(rR);
  rRd->connect({n30, SimNode::GND});

  auto voltageSpeedTermDFluxS =
      VoltageSpeedTerm::make("voltageSpeedTermDFluxS", Logger::Level::debug);
  voltageSpeedTermDFluxS->setIsNegative(true);
  voltageSpeedTermDFluxS->setInitialOmega(omega);
  voltageSpeedTermDFluxS->setIsConstantSpeed(true);
  voltageSpeedTermDFluxS->setInductance(lS);
  voltageSpeedTermDFluxS->connect({n13, n15, n7, n10});

  auto voltageSpeedTermDFluxM =
      VoltageSpeedTerm::make("voltageSpeedTermDFluxM", Logger::Level::debug);
  voltageSpeedTermDFluxM->setIsNegative(true);
  voltageSpeedTermDFluxM->setInitialOmega(omega);
  voltageSpeedTermDFluxM->setIsConstantSpeed(true);
  voltageSpeedTermDFluxM->setInductance(lM);
  voltageSpeedTermDFluxM->connect({n19, n21, n10, n12});

  auto voltageSpeedTermDFluxR =
      VoltageSpeedTerm::make("voltageSpeedTermDFluxR", Logger::Level::debug);
  voltageSpeedTermDFluxR->setIsNegative(false);
  voltageSpeedTermDFluxR->setInitialOmega(omega);
  voltageSpeedTermDFluxR->setIsConstantSpeed(true);
  voltageSpeedTermDFluxR->setInductance(lR);
  voltageSpeedTermDFluxR->connect({n31, n29, n22, n24});

  auto voltageSpeedTermDFluxRM =
      VoltageSpeedTerm::make("voltageSpeedTermDFluxRM", Logger::Level::debug);
  voltageSpeedTermDFluxRM->setIsNegative(false);
  voltageSpeedTermDFluxRM->setInitialOmega(omega);
  voltageSpeedTermDFluxRM->setIsConstantSpeed(true);
  voltageSpeedTermDFluxRM->setInductance(lM);
  voltageSpeedTermDFluxRM->connect({n21, SimNode::GND, n24, n26});

  // q-axis components
  auto rSq = Resistor::make("rSq", Logger::Level::debug);
  rSq->setParameters(rS);
  rSq->connect({n5, n8});
  auto lSq = Inductor::make("lSq", Logger::Level::debug);
  lSq->setParameters(lS);
  lSq->connect({n15, n17});
  auto lMq = Inductor::make("lMq", Logger::Level::debug);
  lMq->setParameters(lM);
  lMq->connect({n17, n19});
  auto lRq = Inductor::make("lRq", Logger::Level::debug);
  lRq->setParameters(lR);
  lRq->connect({n17, n23});
  auto rRq = Resistor::make("rRq", Logger::Level::debug);
  rRq->setParameters(rR);
  rRq->connect({n31, SimNode::GND});

  auto voltageSpeedTermQFluxS =
      VoltageSpeedTerm::make("voltageSpeedTermQFluxS", Logger::Level::debug);
  voltageSpeedTermQFluxS->setInitialOmega(omega);
  voltageSpeedTermQFluxS->setIsConstantSpeed(true);
  voltageSpeedTermQFluxS->setInductance(lS);
  voltageSpeedTermQFluxS->connect({n12, n14, n8, n11});

  auto voltageSpeedTermQFluxM =
      VoltageSpeedTerm::make("voltageSpeedTermQFluxM", Logger::Level::debug);
  voltageSpeedTermQFluxM->setInitialOmega(omega);
  voltageSpeedTermQFluxM->setIsConstantSpeed(true);
  voltageSpeedTermQFluxM->setInductance(lM);
  voltageSpeedTermQFluxM->connect({n18, n20, n11, n13});

  auto voltageSpeedTermQFluxR =
      VoltageSpeedTerm::make("voltageSpeedTermQFluxR", Logger::Level::debug);
  voltageSpeedTermQFluxR->setIsNegative(true);    
  voltageSpeedTermQFluxR->setInitialOmega(omega);
  voltageSpeedTermQFluxR->setIsConstantSpeed(true);
  voltageSpeedTermQFluxR->setInductance(lR);
  voltageSpeedTermQFluxR->connect({n30, n28, n23, n25});

  auto voltageSpeedTermQFluxRM =
      VoltageSpeedTerm::make("voltageSpeedTermQFluxRM", Logger::Level::debug);
  voltageSpeedTermQFluxRM->setIsNegative(true);
  voltageSpeedTermQFluxRM->setInitialOmega(omega);
  voltageSpeedTermQFluxRM->setIsConstantSpeed(true);
  voltageSpeedTermQFluxRM->setInductance(lM);
  voltageSpeedTermQFluxRM->connect({n20, SimNode::GND, n25, n27});

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
  inertiaMoment->connect({n32, SimNode::GND});

  auto EMConverterD =
      ElectroMechanicalConverter::make("EMConverterD", Logger::Level::debug);
  EMConverterD->setInitialFlux(0.0);
  EMConverterD->setIsNegative(true);
  EMConverterD->setVoltageReferenceNode(n23);
  EMConverterD->connect({n26, n28, n32, SimNode::GND});

  auto EMConverterQ =
      ElectroMechanicalConverter::make("EMConverterQ", Logger::Level::debug);
  EMConverterQ->setInitialFlux(0.0);
  EMConverterQ->setIsNegative(false);
  EMConverterQ->setVoltageReferenceNode(n22);
  EMConverterQ->connect({n27, n29, n32, SimNode::GND});

  // Define system topology
  SystemTopology system(
      50, SystemNodeList{n1,  n2,  n3,  n4,  n5,  n6,  n7,  n8,  n9,  n10, n11,
                         n12, n13, n14, n15, n16, n17, n18, n19, n20, n21, n22,
                         n23, n24, n25, n26, n27, n28, n29, n30, n31, n32},
      SystemComponentList{
          v1,
          v2,
          v3,
          parkTrafo,
          rSd,
          lSd,
          lMd,
          lRd,
          rRd,
          voltageSpeedTermDFluxS,
          voltageSpeedTermDFluxM,
          voltageSpeedTermDFluxR,
          voltageSpeedTermDFluxRM,
          rSq,
          lSq,
          lMq,
          lRq,
          rRq,
          voltageSpeedTermQFluxS,
          voltageSpeedTermQFluxM,
          voltageSpeedTermQFluxR,
          voltageSpeedTermQFluxRM,
          r0,
          l0,
          inertiaMoment,
          EMConverterD,
          EMConverterQ,
      });
  // Define simulation scenario
  String simName = "motorStartingTest";
  // Logger
  auto logger = DataLogger::make(simName);
  logger->logAttribute("omega", n32->attribute("v"));
  logger->logAttribute("omega_park", parkTrafo->attribute("omega"));
  logger->logAttribute("I_sd", rSd->attribute("i_intf"));
  logger->logAttribute("I_sq", rSq->attribute("i_intf"));
  logger->logAttribute("I_s0", r0->attribute("i_intf"));
  logger->logAttribute("flux_q", EMConverterD->attribute("flux"));
  logger->logAttribute("flux_d", EMConverterQ->attribute("flux"));

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

void motorStartingTestStatorReference(Real timeStep, Real finalTime,
                                      bool doEigenvalueExtraction) {
  Real frequency = 50;
  Real omega = 2 * M_PI * frequency;
  Real voltageMagnitudeLL = 10e3;
  Real infeedResistance = 3.4e-3;
  Real infeedReactance = 9.4e-3;
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
  Real infeedInductance = infeedReactance / omega;
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
  auto rInfeed1 = Resistor::make("rInfeed", Logger::Level::debug);
  rInfeed1->setParameters(infeedResistance);
  rInfeed1->connect({n4, n1});
  auto rInfeed2 = Resistor::make("rInfeed2", Logger::Level::debug);
  rInfeed2->setParameters(infeedResistance);
  rInfeed2->connect({n5, n2});
  auto rInfeed3 = Resistor::make("rInfeed3", Logger::Level::debug);
  rInfeed3->setParameters(infeedResistance);
  rInfeed3->connect({n6, n3});
  auto lInfeed1 = Inductor::make("lInfeed", Logger::Level::debug);
  lInfeed1->setParameters(infeedInductance);
  lInfeed1->connect({n4, n7});
  auto lInfeed2 = Inductor::make("lInfeed2", Logger::Level::debug);
  lInfeed2->setParameters(infeedInductance);
  lInfeed2->connect({n5, n8});
  auto lInfeed3 = Inductor::make("lInfeed3", Logger::Level::debug);
  lInfeed3->setParameters(infeedInductance);
  lInfeed3->connect({n6, n9});

  // Park transformer
  auto parkTrafo =
      ParkTransformer::make("ParkTransformer", Logger::Level::debug);
  parkTrafo->setIsOmegaConstant(true);
  parkTrafo->setInitialValues(0.0, 0.0);
  parkTrafo->connect({n7, n8, n9, n10, n11, n12});

  // d-axis components
  auto rSd = Resistor::make("rSd", Logger::Level::debug);
  rSd->setParameters(rS);
  rSd->connect({n10, n13});
  auto lSd = Inductor::make("lSd", Logger::Level::debug);
  lSd->setParameters(lS);
  lSd->connect({n13, n16});
  auto lMd = Inductor::make("lMd", Logger::Level::debug);
  lMd->setParameters(lM);
  lMd->connect({n16, SimNode::GND});
  auto lRd = Inductor::make("lRd", Logger::Level::debug);
  lRd->setParameters(lR);
  lRd->connect({n16, n18});
  auto rRd = Resistor::make("rRd", Logger::Level::debug);
  rRd->setParameters(rR);
  rRd->connect({n20, SimNode::GND});
  auto electroMechanicalConverter_d = ElectroMechanicalConverter::make(
      "electroMechanicalConverter_d", Logger::Level::debug);
  electroMechanicalConverter_d->setVoltageReferenceNode(n19);
  electroMechanicalConverter_d->setIsNegative(true);
  electroMechanicalConverter_d->connect({n18, n20, n22, SimNode::GND});

  // q-axis components
  auto rSq = Resistor::make("rSq", Logger::Level::debug);
  rSq->setParameters(rS);
  rSq->connect({n11, n14});
  auto lSq = Inductor::make("lSq", Logger::Level::debug);
  lSq->setParameters(lS);
  lSq->connect({n14, n17});
  auto lMq = Inductor::make("lMq", Logger::Level::debug);
  lMq->setParameters(lM);
  lMq->connect({n17, SimNode::GND});
  auto lRq = Inductor::make("lRq", Logger::Level::debug);
  lRq->setParameters(lR);
  lRq->connect({n17, n19});
  auto rRq = Resistor::make("rRq", Logger::Level::debug);
  rRq->setParameters(rR);
  rRq->connect({n21, SimNode::GND});
  auto electroMechanicalConverter_q = ElectroMechanicalConverter::make(
      "electroMechanicalConverter_q", Logger::Level::debug);
  electroMechanicalConverter_q->setVoltageReferenceNode(n18);
  electroMechanicalConverter_q->connect({n19, n21, n22, SimNode::GND});

  // 0-axis components
  auto r0 = Resistor::make("r0", Logger::Level::debug);
  r0->setParameters(rS);
  r0->connect({n12, n15});
  auto l0 = Inductor::make("l0", Logger::Level::debug);
  l0->setParameters(lS);
  l0->connect({n15, SimNode::GND});

  // Mechanical components
  auto inertiaMoment =
      InertiaMoment::make("inertiaMoment", Logger::Level::debug);
  inertiaMoment->setParameters(inertia);
  inertiaMoment->connect({n22, SimNode::GND});

  // Define system topology
  SystemTopology system(
      50, SystemNodeList{n1,  n2,  n3,  n4,  n5,  n6,  n7,  n8,  n9,  n10, n11,
                         n12, n13, n14, n15, n16, n17, n18, n19, n20, n21, n22},
      SystemComponentList{v1,           v2,
                          v3,           rInfeed1,
                          rInfeed2,     rInfeed3,
                          lInfeed1,     lInfeed2,
                          lInfeed3,     parkTrafo,
                          rSd,          lSd,
                          lMd,          lRd,
                          rRd,          electroMechanicalConverter_d,
                          rSq,          lSq,
                          lMq,          lRq,
                          rRq,          electroMechanicalConverter_q,
                          r0,           l0,
                          inertiaMoment});
  // Define simulation scenario
  String simName = "motorStartingTestStatorReference";
  // Logger
  auto logger = DataLogger::make(simName);
  logger->logAttribute("omega", n22->attribute("v"));
  // logger->logAttribute("Flux_d",
  //                    electroMechanicalConverter_d->attribute("flux"));
  // logger->logAttribute("Flux_q",
  //                    electroMechanicalConverter_q->attribute("flux"));
  // logger->logAttribute("omega_Park", parkTrafo->attribute("omega"));
  // logger->logAttribute("V_A", n1->attribute("v"));
  // logger->logAttribute("V_B", n2->attribute("v"));
  // logger->logAttribute("V_C", n3->attribute("v"));
  logger->logAttribute("I_A", rInfeed1->attribute("i_intf"));
  logger->logAttribute("I_B", rInfeed2->attribute("i_intf"));
  logger->logAttribute("I_C", rInfeed3->attribute("i_intf"));
  // logger->logAttribute("V_D", n10->attribute("v"));
  // logger->logAttribute("V_Q", n11->attribute("v"));
  // logger->logAttribute("V_0", n12->attribute("v"));
  // logger->logAttribute("I_sd", rSd->attribute("i_intf"));
  // logger->logAttribute("I_sq", rSq->attribute("i_intf"));
  // logger->logAttribute("Torque_d",
  //                    electroMechanicalConverter_d->attribute("torque"));
  // logger->logAttribute("Torque_q",
  //                    electroMechanicalConverter_q->attribute("torque"));

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

int main(int argc, char *argv[]) {
    timeLaggingTorqueSourceTest();
    // timeLaggingVoltageSourceTest();
    // motorStartingTestSynchReferenceSpeedVoltageTermsEMConverter(1e-4, 3, true);
    // motorStartingTestSynchReferenceSpeedVoltageTerms(1e-4, 3, false);
//   motorStartingTestRotorReferenceSpeedVoltageTerms(1e-4, 3, false);
  //  speedVoltageTermTest(1e-4, 1, true);
  // motorStartingTestStatorReference(1e-5, 3.0, true);
  // motorStartingTestRotorReference(1e-4, 3.0, true);
  // syncronousGeneratorTest(1e-3, 70.0, false);
  // electromechanicalConverterTest(1e-4, 1.0, true);

  return 0;
}
