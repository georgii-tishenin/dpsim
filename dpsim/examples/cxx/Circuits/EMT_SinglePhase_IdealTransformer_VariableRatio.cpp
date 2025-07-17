/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <DPsim.h>

using namespace DPsim;
using namespace CPS::EMT;
using namespace CPS::EMT::Ph1;

/*
void Ideal_Transformer_VariableRatio() {


  // Nodes
  auto n1 = SimNode::make("n1");
  auto n2 = SimNode::make("n2");
  auto n3 = SimNode::make("n3");
  auto n4 = SimNode::make("n4");
 // auto n5 = SimNode::make("n5");

  // Components
  auto v1 = VoltageSource::make("v1", Logger::Level::debug);
  v1->setParameters(Complex(10, 0));
  
  auto t1 = IdealTransformerVariableRatio::make("t1", Logger::Level::debug);

  t1->setVoltageSource(v1);
  t1->setTimeStep(0.00005);
  //t1->setParameters(1);
  

  auto r1 = Resistor::make("r1", Logger::Level::debug);
  r1->setParameters(1);

  auto r_gnd = Resistor::make("r_gnd", Logger::Level::debug);
    r_gnd->setParameters(1);

  auto r_gnd2 = Resistor::make("r_gnd2", Logger::Level::debug);
    r_gnd2->setParameters(1);
  
 // auto l1 = Inductor::make("l1", Logger::Level::debug);
 // l1->setParameters(0.1);
  


  // Topology
  v1->connect({n1, n2});
  
  t1->connect({n1, n2, n3, n4});
  
  r1->connect({n3, n4});
  
 // l1->connect({n4, n5});

  //r_gnd->connect({SimNode::GND, n4});
  //r_gnd2->connect({SimNode::GND, n1});

  r_gnd->connect({n4, SimNode::GND});
  r_gnd2->connect({n1, SimNode::GND});

  // Define system topology
  SystemTopology system(50, SystemNodeList{
    n1, n2, n3, n4, SimNode::GND},
    SystemComponentList{
      v1, t1, r1, r_gnd, r_gnd2});


  // Define simulation scenario
  Real timeStep = 0.00005;
  Real finalTime = 0.1;
  String simName = "EMT_IdealTransformer_VariableRatio" + std::to_string(timeStep);


    // Logger
  auto logger = DataLogger::make(simName);
  logger->logAttribute("v_1", n2->attribute("v"));
  logger->logAttribute("v_r", r1->attribute("v_intf"));
  //logger->logAttribute("v_l", l1->attribute("v_intf"));
  logger->logAttribute("i_1", v1->attribute("i_intf"));
  logger->logAttribute("i_2", r1->attribute("i_intf"));
  logger->logAttribute("Ratio", t1->attribute("Ratio"));
  logger->logAttribute("v_2", n3->attribute("v"));
  logger->logAttribute("Test", n1->attribute("v"));
  logger->logAttribute("i_gnd", r_gnd->attribute("i_intf"));

  Simulation sim(simName);
  sim.setSystem(system);
  sim.setTimeStep(timeStep);
  sim.doSystemMatrixRecomputation(true);
  sim.setFinalTime(finalTime);
  sim.setDomain(Domain::EMT);
  sim.doEigenvalueExtraction(true);
  sim.addLogger(logger);
  sim.run();
}

void Ideal_Transformer_VariableRatio_Grounded() {


  // Nodes
  auto n1 = SimNode::make("n1");
 // auto n2 = SimNode::make("n2");
  auto n3 = SimNode::make("n3");
  auto n4 = SimNode::make("n4");
  auto n5 = SimNode::make("n5");

  // Components
  auto v1 = VoltageSource::make("v1", Logger::Level::debug);
  v1->setParameters(Complex(10, 0));
  
  auto t1 = IdealTransformerVariableRatio::make("t1", Logger::Level::debug);

  t1->setVoltageSource(v1);
  t1->setTimeStep(0.00005);
  //t1->setParameters(1);
  

  auto r1 = Resistor::make("r1", Logger::Level::debug);
  r1->setParameters(1);

 // auto r_gnd = Resistor::make("r_gnd", Logger::Level::debug);
 //   r_gnd->setParameters(100000);

 // auto r_gnd2 = Resistor::make("r_gnd2", Logger::Level::debug);
 // r_gnd2->setParameters(100000);
  
  auto l1 = Inductor::make("l1", Logger::Level::debug);
  l1->setParameters(0.1);
  


  // Topology
  v1->connect({n1, SimNode::GND});
  
  t1->connect({n1, SimNode::GND, n3, SimNode::GND});
  
  r1->connect({n3, n4});
  
  l1->connect({n4, SimNode::GND});

  //r_gnd->connect({n4, SimNode::GND});

  //r_gnd2->connect({n1, SimNode::GND});

  // Define system topology
  SystemTopology system(50, SystemNodeList{
    n1, n3, n4, SimNode::GND},
    SystemComponentList{
      v1, l1, t1, r1});


  // Define simulation scenario
  Real timeStep = 0.00005;
  Real finalTime = 0.1;
  String simName = "EMT_IdealTransformer_VariableRatio" + std::to_string(timeStep);


    // Logger
  auto logger = DataLogger::make(simName);
  logger->logAttribute("v_1", v1->attribute("v_intf"));
  logger->logAttribute("v_r", r1->attribute("v_intf"));
  logger->logAttribute("v_l", l1->attribute("v_intf"));
  logger->logAttribute("i_1", v1->attribute("i_intf"));
  logger->logAttribute("i_2", r1->attribute("i_intf"));
  logger->logAttribute("Ratio", t1->attribute("Ratio"));
  logger->logAttribute("Test", n3->attribute("v"));

  Simulation sim(simName);
  sim.setSystem(system);
  sim.setTimeStep(timeStep);
  sim.doSystemMatrixRecomputation(true);
  sim.setFinalTime(finalTime);
  sim.setDomain(Domain::EMT);
 // sim.doEigenvalueExtraction(true);
  sim.addLogger(logger);
  sim.run();
}


void Ideal_Transformer_VariableRatio_Eigenvalues() {


  // Nodes
  auto n1 = SimNode::make("n1");
  auto n2 = SimNode::make("n2");
  auto n3 = SimNode::make("n3");
  auto n4 = SimNode::make("n4");
  auto n5 = SimNode::make("n5");

  // Components
  auto v1 = VoltageSource::make("v1", Logger::Level::debug);
  v1->setParameters(Complex(10, 0));
  
  auto t1 = IdealTransformerVariableRatio::make("t1", Logger::Level::debug);

  t1->setVoltageSource(v1);
  t1->setTimeStep(0.00005);
  

  auto r1 = Resistor::make("r1", Logger::Level::debug);
  r1->setParameters(1);

  auto r_gnd = Resistor::make("r_gnd", Logger::Level::debug);
    r_gnd->setParameters(1);

  auto r_gnd2 = Resistor::make("r_gnd2", Logger::Level::debug);
    r_gnd2->setParameters(1);
  
  auto l1 = Inductor::make("l1", Logger::Level::debug);
  l1->setParameters(0.1);
  


  // Topology
  v1->connect({n1, n2});
  
  t1->connect({n1, n2, n3, n5});
  
  r1->connect({n3, n4});
  
  l1->connect({n4, n5});

  //r_gnd->connect({SimNode::GND, n4});
  //r_gnd2->connect({SimNode::GND, n1});

  r_gnd->connect({n4, SimNode::GND});
  r_gnd2->connect({n1, SimNode::GND});

  // Define system topology
  SystemTopology system(50, SystemNodeList{
    n1, n2, n3, n4, n5, SimNode::GND},
    SystemComponentList{
      v1, t1, r1, r_gnd, r_gnd2, l1});


  // Define simulation scenario
  Real timeStep = 0.00005;
  Real finalTime = 0.1;
  String simName = "EMT_IdealTransformer_VariableRatio" + std::to_string(timeStep);


    // Logger
  auto logger = DataLogger::make(simName);
  logger->logAttribute("v_1", n2->attribute("v"));
  logger->logAttribute("v_r", r1->attribute("v_intf"));
  //logger->logAttribute("v_l", l1->attribute("v_intf"));
  logger->logAttribute("i_1", v1->attribute("i_intf"));
  logger->logAttribute("i_2", r1->attribute("i_intf"));
  logger->logAttribute("Ratio", t1->attribute("Ratio"));
  logger->logAttribute("v_2", n3->attribute("v"));
  logger->logAttribute("Test", n1->attribute("v"));
  logger->logAttribute("i_gnd", r_gnd->attribute("i_intf"));

  Simulation sim(simName);
  sim.setSystem(system);
  sim.setTimeStep(timeStep);
  sim.doSystemMatrixRecomputation(true);
  sim.setFinalTime(finalTime);
  sim.setDomain(Domain::EMT);
  sim.doEigenvalueExtraction(true);
  sim.addLogger(logger);
  sim.run();
}


*/


int main(int argc, char *argv[]) {

// Ideal_Transformer_VariableRatio();

// Ideal_Transformer_VariableRatio_Grounded();

// Ideal_Transformer_VariableRatio_Eigenvalues();

    return 0;
}
