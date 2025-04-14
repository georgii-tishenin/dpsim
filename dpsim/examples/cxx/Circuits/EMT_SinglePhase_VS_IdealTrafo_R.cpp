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

int main(int argc, char *argv[]) {
  // Nodes
  auto n1 = SimNode::make("n1");
  auto n2 = SimNode::make("n2");
  auto n3 = SimNode::make("n3");
  auto n4 = SimNode::make("n4");

  // Components
  auto v1 = VoltageSource::make("v1", Logger::Level::debug);
  v1->setParameters(Complex(8.5, 0));
  //	auto l1 = Inductor::make("l_1");
  //	auto r2 = Resistor::make("r_2");
  
  
  //auto t1 = IdealTrafo::make("t1", Logger::Level::debug);
  //t1->setParameters(50000, 1000, 50, 0, 0, 0);

auto t1 = IdealTransformer::make("t1", Logger::Level::debug);
t1->setParameters(10);




  auto r1 = Resistor::make("r1", Logger::Level::debug);
  r1->setParameters(100);
  auto l1 = Inductor::make("l1", Logger::Level::debug);
  l1->setParameters(5);
  auto l2 = Inductor::make("l2", Logger::Level::debug);
  l2->setParameters(5);
  auto c1 = Capacitor::make("c1", Logger::Level::debug);
  c1->setParameters(250e-6);

  // Topology
  v1->connect({SimNode::GND, n1});
  //	l1->connect({ n1, n2 });
  //	r2->connect({ n2, SimNode::GND });
  t1->connect({n2, n3});
  r1->connect({n1, n2});
  l1->connect({n3, n4});
  l2->connect({n2, SimNode::GND});
  c1->connect({n4, SimNode::GND});


  // Define system topology
  SystemTopology system(50, SystemNodeList{n1, n2, n3, n4,  SimNode::GND},
                        SystemComponentList{v1, t1, r1, l1, c1, l2});

  // Define simulation scenario
  Real timeStep = 0.00005;
  Real finalTime = 0.4;
  String simName = "EMT_IdealTrafo_R_" + std::to_string(timeStep);

  //Simulation sim(simName, system, timeStep, finalTime);
  
   // Logger
  auto logger = DataLogger::make(simName);
  logger->logAttribute("v1", n1->attribute("v"));
  logger->logAttribute("v2", n2->attribute("v"));
  logger->logAttribute("i", r1->attribute("i_intf"));
  logger->logAttribute("virtual voltage", t1->virtualNode(0)->attribute("v"));
  logger->logAttribute("inductor current", l1->attribute("i_intf"));
  
  Simulation sim(simName);
  sim.setSystem(system);
  sim.setTimeStep(timeStep);
  sim.setFinalTime(finalTime);
  sim.setDomain(Domain::EMT);
  sim.doEigenvalueExtraction(true);
  sim.addLogger(logger);
  sim.run();

  return 0;
}
