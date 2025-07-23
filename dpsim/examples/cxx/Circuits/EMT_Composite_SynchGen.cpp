#include <DPsim.h>
#include <iostream>
#include <string>
#include <vector>

using namespace DPsim;
using namespace CPS::EMT;
using namespace CPS::EMT::Ph1;

void CompositeSynchGen() {

  Real timeStep = 5e-4;

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

  //Components
  auto v1 = VoltageSource::make("v1", Logger::Level::debug);
  v1->setParameters(Complex(10.0, 0.0), 50);
  auto v2 = VoltageSource::make("v2", Logger::Level::debug);
  v2->setParameters(Complex(-5.0, -5.0 * sqrt(3.0)), 50);
  auto v3 = VoltageSource::make("v3", Logger::Level::debug);
  v3->setParameters(Complex(-5.0, 5.0 * sqrt(3.0)), 50);

  auto R0 =
      Resistor::make("R0", Logger::Level::debug); // zero sequence resistor
  R0->setParameters(1e3);

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

  auto ElectroMechanicalConverter_d = IdealTransformerVariableRatio::make(
      "ElectroMechanicalConverter_d", Logger::Level::debug);
  ElectroMechanicalConverter_d->setStatorInductor(Ls_sigma_q);
  ElectroMechanicalConverter_d->setTimeStep(timeStep);

  auto ElectroMechanicalConverter_q = IdealTransformerVariableRatio::make(
      "ElectroMechanicalConverter_q", Logger::Level::debug);
  ElectroMechanicalConverter_q->setStatorInductor(Ls_sigma_d);
  ElectroMechanicalConverter_q->setTimeStep(timeStep);

  auto ParkTrafo =
      ParkTransformer::make("ParkTransformer", Logger::Level::debug);
  ParkTrafo->isOmegaConstant(
      true); // Set the Park transformer to use a constant angular frequency
  ParkTrafo->setParameters(2 * M_PI * 50,
                           0.0); // Set the angular frequency and initial angle

  //Mechanical components
  auto InertiaMoment =
      InertiaMoment::make("InertiaMoment", Logger::Level::debug);
  InertiaMoment->setParameters(0.0325); // Set the inertia moment value

  auto LoadTorque = CurrentSource::make("LoadTorque", Logger::Level::debug);
  LoadTorque->setParameters(Complex(400, 0.0));

  auto ConstantOmegaSource =
      VoltageSource::make("ConstantOmegaSource", Logger::Level::debug);
  ConstantOmegaSource->setParameters(
      Complex(2 * M_PI * 50, 0.0)); // Set the constant angular frequency source

  //Topology
  // abc side
  v1->connect({n1, SimNode::GND});
  v2->connect({n2, SimNode::GND});
  v3->connect({n3, SimNode::GND});
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
  R0->connect(
      {n17, SimNode::GND}); // zero sequence resistor connects to node n6

  //Mechanical side
  InertiaMoment->connect({n18, SimNode::GND});
  LoadTorque->connect({n18, SimNode::GND});
  ConstantOmegaSource->connect({SimNode::GND, n18});

  // Define system topology
  SystemTopology system(50,
                        SystemNodeList{n1, n2, n3, n4, n5, n6, n7, n8, n9, n10,
                                       n11, n12, n13, n14, n15, n16, n17, n18},
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
                                            v1,
                                            v2,
                                            v3,
                                            R0,
                                            Lf_d,
                                            LoadTorque,
                                            ConstantOmegaSource});

  // Define simulation scenario
  Real finalTime = 5.0;
  String simName = "EMT_CompositeSynchGen" + std::to_string(timeStep);

  // Logger
  auto logger = DataLogger::make(simName);
  logger->logAttribute("omega", n18->attribute("v"));
  logger->logAttribute("Torque", InertiaMoment->attribute("i_intf"));
  logger->logAttribute("Ratio_d",
                       ElectroMechanicalConverter_d->attribute("Ratio"));
  logger->logAttribute("Ratio_q",
                       ElectroMechanicalConverter_q->attribute("Ratio"));
  logger->logAttribute("V_lsd", Ls_sigma_d->attribute("v_intf"));
  logger->logAttribute("V_lsq", Ls_sigma_q->attribute("v_intf"));
  logger->logAttribute("V_m_d", Lmd->attribute("v_intf"));
  logger->logAttribute("V_m_q", Lmq->attribute("v_intf"));
  logger->logAttribute("V_f_d", Vf_d->attribute("v_intf"));
  logger->logAttribute("I_f_d", Vf_d->attribute("i_intf"));
  logger->logAttribute("V_D", n4->attribute("v"));
  logger->logAttribute("V_Q", n11->attribute("v"));
  logger->logAttribute("V_R0", R0->attribute("v_intf"));
  logger->logAttribute("VN7", n7->attribute("v"));
  logger->logAttribute("VN14", n14->attribute("v"));
  logger->logAttribute("i_Rfd", Rf_d->attribute("i_intf"));
  logger->logAttribute("v_Rfd", Rf_d->attribute("v_intf"));
  logger->logAttribute("v_Lfd", Lf_d->attribute("v_intf"));
  logger->logAttribute("Load", LoadTorque->attribute("i_intf"));
  logger->logAttribute("I_sd", Rsd->attribute("i_intf"));
  logger->logAttribute("I_sq", Rsq->attribute("i_intf"));

  Simulation sim(simName);
  sim.setSystem(system);
  sim.setTimeStep(timeStep);
  sim.doSystemMatrixRecomputation(true);
  sim.setFinalTime(finalTime);
  sim.setDomain(Domain::EMT);
  sim.addLogger(logger);
  sim.run();
}

int main(int argc, char *argv[]) {

  CompositeSynchGen();
  return 0;
}
