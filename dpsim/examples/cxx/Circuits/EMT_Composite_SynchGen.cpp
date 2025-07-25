#include <DPsim.h>
#include <iostream>
#include <string>
#include <vector>

using namespace DPsim;
using namespace CPS::EMT;
using namespace CPS::EMT::Ph1;

void CompositeSynchGen(Real timeStep, Real finalTime, bool doEigenvalueExtraction) {

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
  ParkTrafo->isOmegaConstant(
      true); // Set the Park transformer to use a constant angular frequency
  ParkTrafo->setParameters(2 * M_PI * 50,
                           M_PI/2); // Set the angular frequency and initial angle

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
                                       n11, n12, n13, n14, n15, n16, n17, n18, n19},
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
                                            ConstantOmegaSource
                                            });

  // Define simulation scenario
  String simName = "EMT_CompositeSynchGen" + std::to_string(timeStep);

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

void EMConverterEigenvaluesTest(Real timeStep, Real finalTime, bool doEigenvalueExtraction) {

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

  auto converter = ElectroMechanicalConverter::make("converter",
                                                         Logger::Level::debug);
  converter->setVoltageReferenceNode(n1);
  converter->connect({n3, n2, n4, SimNode::GND});

  auto inertiaMoment =
      InertiaMoment::make("inertiaMoment", Logger::Level::debug);
      inertiaMoment->setParameters(1e-3);
    inertiaMoment->connect({n4, SimNode::GND});

  // Define system topology
    SystemTopology system(50,
                            SystemNodeList{n1, n2, n3, n4},
                            SystemComponentList{v, l, r, converter, inertiaMoment});

  // Define simulation scenario
  String simName = "EMConverterEigenvaluesTest";
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

int main(int argc, char *argv[]) {

// CompositeSynchGen(1e-3, 10.0, false);
EMConverterEigenvaluesTest(1e-4, 1.0, true);

  return 0;
}
