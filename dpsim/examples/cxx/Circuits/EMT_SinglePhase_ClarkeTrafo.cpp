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

  // Components
  auto v1 = VoltageSource::make("v1", Logger::Level::debug);
    v1->setParameters(Complex(10, 0), 50);
/*
    auto v2 = VoltageSource::make("v2", Logger::Level::debug);
    v2->setParameters(Complex(10, 0), 50);
    auto v3 = VoltageSource::make("v3", Logger::Level::debug);
    v3->setParameters(Complex(10, 0), 50);
      */
  
  
    auto v2 = VoltageSource::make("v2", Logger::Level::debug);
    v2->setParameters(Complex(-5, -5*sqrt(3.0)), 50);
  auto v3 = VoltageSource::make("v3", Logger::Level::debug);
    v3->setParameters(Complex(-5, 5*sqrt(3.0)), 50);

  auto alphaBetaTrafo = ClarkeTransformer::make("ClarkeTrafo", Logger::Level::debug);
  //  alphaBetaTrafo->setParameters(true);
    
  auto r1 = Resistor::make("r1", Logger::Level::debug);
    r1->setParameters(1);
  auto r2 = Resistor::make("r2", Logger::Level::debug);
    r2->setParameters(1);
  auto r3 = Resistor::make("r3", Logger::Level::debug);
    r3->setParameters(1);

 
  auto r4 = Resistor::make("r4", Logger::Level::debug);
    r4->setParameters(50);
  auto r5 = Resistor::make("r5", Logger::Level::debug); 
    r5->setParameters(50);
  
  auto r6 = Resistor::make("r6", Logger::Level::debug);
    r6->setParameters(50);

    auto r_neutral = Resistor::make("r_neutral", Logger::Level::debug);
    r_neutral->setParameters(3*51);
  
  
  auto l1 = Inductor::make("l1", Logger::Level::debug);
    l1->setParameters(0.1);
  auto l2 = Inductor::make("l2", Logger::Level::debug);
    l2->setParameters(0.1);
  auto l3 = Inductor::make("l3", Logger::Level::debug);
   l3->setParameters(0.1);


  // Topology
  v1->connect({SimNode::GND, n1});
  v2->connect({SimNode::GND, n2});
  v3->connect({SimNode::GND, n3});

  alphaBetaTrafo->connect({n4, n5, n6, n7, n8, n9});


  r1->connect({n1, n4});
  r2->connect({n2, n5});
  r3->connect({n3, n6});
  l1->connect({n10, SimNode::GND});
  l2->connect({n11, SimNode::GND}); 
  l3->connect({n12, n13});

 r4->connect({n10, n7});
 r5->connect({n11, n8});
 r6->connect({n12, n9});
 r_neutral->connect({n13, SimNode::GND});





  // Define system topology
  SystemTopology system(50, SystemNodeList{n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, SimNode::GND},
                        SystemComponentList{v1, v2, v3, alphaBetaTrafo, r1, r2, r3, r4, r5, r6, l1, l2, l3, r_neutral});

  // Define simulation scenario
  Real timeStep = 0.00005;
  Real finalTime = 0.2;
  String simName = "EMT_alphaBetaTrafo_R_" + std::to_string(timeStep);

  //Simulation sim(simName, system, timeStep, finalTime);
  
   // Logger
  auto logger = DataLogger::make(simName);

  
  logger->logAttribute("va", n4->attribute("v"));
  logger->logAttribute("vb", n5->attribute("v"));
  logger->logAttribute("vc", n6->attribute("v"));
  logger->logAttribute("i_1", r1->attribute("i_intf"));
  logger->logAttribute("i_2", r2->attribute("i_intf"));
  logger->logAttribute("i_3", r3->attribute("i_intf"));
  logger->logAttribute("v_alpha", n7->attribute("v"));
  logger->logAttribute("v_beta", n8->attribute("v"));
  logger->logAttribute("v_0", n9->attribute("v"));
  logger->logAttribute("i_alpha", r4->attribute("i_intf"));
  logger->logAttribute("i_beta", r5->attribute("i_intf"));
  logger->logAttribute("i_0", r6->attribute("i_intf"));
  
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
