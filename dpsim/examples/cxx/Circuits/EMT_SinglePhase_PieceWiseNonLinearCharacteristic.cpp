#include <DPsim.h>
#include <iostream>
#include <string>
#include <vector>

using namespace DPsim;
using namespace CPS::EMT;
using namespace CPS::EMT::Ph1;


void testNonLinearCharacteristic() {
    // Define boundary points for the positive quadrant.
    std::vector<CPS::EMT::Ph1::PieceWiseNonLinearCharacteristic::Point> pts = {
        {0.0, 0.0},
        {0.2, 0.4},
        {0.4, 0.7},
        {1.0, 1.0}
    };

    // Create an instance of the characteristic with saturated inductance 0.001.
    CPS::EMT::Ph1::PieceWiseNonLinearCharacteristic characteristic(pts, 0.001);

    std::cout << "Testing PieceWiseNonLinearCharacteristic\n";
    std::cout << "Enter flux linkage values (or type 'exit' to quit):\n";

    std::string input;
    while (true) {
        std::cout << "Flux: ";
        std::getline(std::cin, input);
        if (input == "exit") {
            break;
        }
        try {
            double flux = std::stod(input);
            double inductance = characteristic.getInductance(flux);
            double current = characteristic.getCurrent(flux);
            std::cout << "For flux linkage " << flux << ":\n";
            std::cout << "  Inductance = " << inductance << "\n";
            std::cout << "  Current    = " << current << "\n";
        } catch (const std::exception &e) {
            std::cout << "Error: " << e.what() << "\n";
        }
    }
    std::cout << "Exiting test." << std::endl;
}


void testNonLinearInductor() {

    // Nodes
    auto n1 = SimNode::make("n1");

    //Components
    auto r1 = Resistor::make("r1", Logger::Level::debug);
    r1->setParameters(1.0);


    // Create a non-linear inductor with a piecewise characteristic.
    auto non_linear_inductor = CPS::EMT::Ph1::NonLinearInductor::make("NonLinearInductor", Logger::Level::debug);
    
    // Define the piecewise characteristic.
    std::vector<CPS::EMT::Ph1::PieceWiseNonLinearCharacteristic::Point> pts = {
        {0.0, 0.0},
        {0.2, 0.4},
        {0.4, 0.7},
        {1.0, 1.0}
    };

    // Create an instance of the characteristic with saturated inductance 0.001.
    auto characteristic = std::make_shared<CPS::EMT::Ph1::PieceWiseNonLinearCharacteristic>(pts, 0.001);
    non_linear_inductor->setPieceWiseCharacteristic(characteristic);

    
    // Set parameters for the inductor.
    non_linear_inductor->setParameters(Complex(0, 0));


    //Topology
    non_linear_inductor->connect({SimNode::GND, n1});
    r1->connect({n1, SimNode::GND});


      // Define system topology
  SystemTopology system(50, 
    SystemNodeList{n1, SimNode::GND},
    SystemComponentList{non_linear_inductor, r1});


  // Define simulation scenario
  Real timeStep = 0.00005;
  Real finalTime = 1.1;
  String simName = "EMT_NonLinearInductor" + std::to_string(timeStep);


    // Logger
  auto logger = DataLogger::make(simName);
  logger->logAttribute("v1", n1->attribute("v"));
  logger->logAttribute("v_r", r1->attribute("v_intf"));
  //logger->logAttribute("v_l", l1->attribute("v_intf"));
  logger->logAttribute("i_1", r1->attribute("i_intf"));
  logger->logAttribute("i_inductor", non_linear_inductor->attribute("i_intf"));

  Simulation sim(simName);
  sim.setSystem(system);
  sim.setTimeStep(timeStep);
  //sim.doSystemMatrixRecomputation(true);
  sim.setFinalTime(finalTime);
  sim.setDomain(Domain::EMT);
  //sim.doEigenvalueExtraction(true);
  sim.addLogger(logger);
  sim.run();
    
    
}




void test() {

    // Nodes
    auto n1 = SimNode::make("n1");
    //auto n2 = SimNode::make("n2");

    //Components
    auto r1 = Resistor::make("r1", Logger::Level::debug);
    r1->setParameters(1.0);
    //auto r2 = Resistor::make("r2", Logger::Level::debug);
    //r2->setParameters(1.0);
    //auto v1 = VoltageSource::make("v1", Logger::Level::debug);
    //v1->setParameters(Complex(10.0, 0.0), 0); 
    auto i1 = CurrentSource::make("i1", Logger::Level::debug);
    i1->setParameters(Complex(10.0, 0.0), 0); // Set current source to zero for now



    //Topology
    r1->connect({SimNode::GND, n1});
    //r2->connect({n2, n1});
    i1->connect({n1 ,SimNode::GND});


      // Define system topology
  SystemTopology system(50, 
    SystemNodeList{n1, SimNode::GND},
    SystemComponentList{i1, r1});


  // Define simulation scenario
  Real timeStep = 0.00005;
  Real finalTime = 1.1;
  String simName = "EMT_NonLinearInductor" + std::to_string(timeStep);


    // Logger
  auto logger = DataLogger::make(simName);
  logger->logAttribute("v1", n1->attribute("v"));
  logger->logAttribute("v_r1", r1->attribute("v_intf"));
  //logger->logAttribute("v_r2", r2->attribute("v_intf"));
  logger->logAttribute("i_1", r1->attribute("i_intf"));
  logger->logAttribute("v_source", i1->attribute("v_intf"));
  logger->logAttribute("i_source", i1->attribute("i_intf"));
  //logger->logAttribute("v2", n2->attribute("v"));
  //logger->logAttribute("i2", r2->attribute("i_intf"));

  Simulation sim(simName);
  sim.setSystem(system);
  sim.setTimeStep(timeStep);
  //sim.doSystemMatrixRecomputation(true);
  sim.setFinalTime(finalTime);
  sim.setDomain(Domain::EMT);
  //sim.doEigenvalueExtraction(true);
  sim.addLogger(logger);
  sim.run();
    
    
}

void testNonLinearInductor_complete() {

    // Nodes
    auto n1 = SimNode::make("n1");

    //Components
    auto v1 = VoltageSource::make("v1", Logger::Level::debug);
    v1->setParameters(Complex(400.0, 0.0), 50); // Set voltage source to 10V



    // Create a non-linear inductor with a piecewise characteristic.
    auto non_linear_inductor = CPS::EMT::Ph1::NonLinearInductor::make("NonLinearInductor", Logger::Level::debug);
    
    // Define the piecewise characteristic.
    std::vector<CPS::EMT::Ph1::PieceWiseNonLinearCharacteristic::Point> pts = {
        {0.0, 0.0},
        {0.01, 0.5},
        {0.1, 0.9},
        {100.0, 1.0}
    };

    // Create an instance of the characteristic with saturated inductance 0.001.
    auto characteristic = std::make_shared<CPS::EMT::Ph1::PieceWiseNonLinearCharacteristic>(pts, 0.001);
    non_linear_inductor->setPieceWiseCharacteristic(characteristic);

    
    // Set parameters for the inductor.
    non_linear_inductor->setParameters(Complex(0, 0));
    non_linear_inductor->setTimeStep(1e-6); // Set time step for the discrete integration


    //Topology
    v1->connect({SimNode::GND, n1});
    non_linear_inductor->connect({SimNode::GND, n1});


      // Define system topology
  SystemTopology system(50, 
    SystemNodeList{n1, SimNode::GND},
    SystemComponentList{non_linear_inductor, v1});


  // Define simulation scenario
  Real timeStep = 1e-6;
  Real finalTime = 0.4;
  String simName = "EMT_NonLinearInductor" + std::to_string(timeStep);


    // Logger
  auto logger = DataLogger::make(simName);
  logger->logAttribute("v1", n1->attribute("v"));
  logger->logAttribute("v_ind", non_linear_inductor->attribute("v_intf"));
  //logger->logAttribute("v_l", l1->attribute("v_intf"));
  logger->logAttribute("i_inductor", non_linear_inductor->attribute("i_intf"));
  logger->logAttribute("v_source", v1->attribute("v_intf"));
  logger->logAttribute("i_source", v1->attribute("i_intf"));

  Simulation sim(simName);
  sim.setSystem(system);
  sim.setTimeStep(timeStep);
  //sim.doSystemMatrixRecomputation(true);
  sim.setFinalTime(finalTime);
  sim.setDomain(Domain::EMT);
  //sim.doEigenvalueExtraction(true);
  sim.addLogger(logger);
  sim.run();
    
    
}



void testInductor() {

    // Nodes
    auto n1 = SimNode::make("n1");

    //Components
    auto v1 = VoltageSource::make("v1", Logger::Level::debug);
    v1->setParameters(Complex(1.0, 0.0), 0); // Set voltage source to 10V
    auto r1 = Resistor::make("r1", Logger::Level::debug);
    r1->setParameters(1.0);
    auto l1 = Inductor::make("l1", Logger::Level::debug);
    l1->setParameters(1); // Set inductance to 1 H


  

    //Topology
    v1->connect({SimNode::GND, n1});
    l1->connect({SimNode::GND, n1});
    r1->connect({SimNode::GND, n1});


      // Define system topology
  SystemTopology system(50, 
    SystemNodeList{n1, SimNode::GND},
    SystemComponentList{l1, r1, v1});


  // Define simulation scenario
  Real timeStep = 0.00005;
  Real finalTime = 1.1;
  String simName = "EMT_NonLinearInductor" + std::to_string(timeStep);


    // Logger
  auto logger = DataLogger::make(simName);
  logger->logAttribute("v1", n1->attribute("v"));
  logger->logAttribute("v_res", r1->attribute("v_intf"));
  logger->logAttribute("v_ind", l1->attribute("v_intf"));
  //logger->logAttribute("v_l", l1->attribute("v_intf"));
  logger->logAttribute("i_res", r1->attribute("i_intf"));
  logger->logAttribute("i_inductor", l1->attribute("i_intf"));
  logger->logAttribute("v_source", v1->attribute("v_intf"));
  logger->logAttribute("i_source", v1->attribute("i_intf"));

  Simulation sim(simName);
  sim.setSystem(system);
  sim.setTimeStep(timeStep);
  //sim.doSystemMatrixRecomputation(true);
  sim.setFinalTime(finalTime);
  sim.setDomain(Domain::EMT);
  //sim.doEigenvalueExtraction(true);
  sim.addLogger(logger);
  sim.run();
    
    
}


void testInrush() {


    // Nodes
    auto n1 = SimNode::make("n1");
    auto n2 = SimNode::make("n2");
    auto n3 = SimNode::make("n3");
    auto n4 = SimNode::make("n4");

    //Components
    auto v1 = VoltageSource::make("v1", Logger::Level::debug);
    v1->setParameters(Complex(400.0, 0.0), 50); // Set voltage source to 10V
    auto r1 = Resistor::make("r1", Logger::Level::debug);
    r1->setParameters(0.1);
    auto r2 = Resistor::make("r2", Logger::Level::debug);
    r2->setParameters(100);
    auto s1 = Switch::make("s", Logger::Level::info);
    s1->setParameters(1e8, 1e-4, false);
    auto l1 = Inductor::make("l1", Logger::Level::debug);
    l1->setParameters(0.001); // Set inductance to 1 mH







    // Create a non-linear inductor with a piecewise characteristic.
    auto non_linear_inductor = CPS::EMT::Ph1::NonLinearInductor::make("NonLinearInductor", Logger::Level::debug);
    
    // Define the piecewise characteristic.
    std::vector<CPS::EMT::Ph1::PieceWiseNonLinearCharacteristic::Point> pts = {
        {0.0, 0.0},
        {0.01, 0.5},
        {0.1, 0.9},
        {100.0, 1.0}
    };

    // Create an instance of the characteristic with saturated inductance 0.001.
    auto characteristic = std::make_shared<CPS::EMT::Ph1::PieceWiseNonLinearCharacteristic>(pts, 0.001);
    non_linear_inductor->setPieceWiseCharacteristic(characteristic);

    
    // Set parameters for the inductor.
    non_linear_inductor->setParameters(Complex(0, 0));
    non_linear_inductor->setTimeStep(1e-6); // Set time step for the discrete integration


    // Create a switch event to close the switch after 0.25 seconds
    auto sEvent = SwitchEvent::make(0.25, s1, true);



    //Topology
    v1->connect({SimNode::GND, n1});
    l1->connect({n3, n4});
    r1->connect({n2, n3});
    r2->connect({n1, SimNode::GND});
    non_linear_inductor->connect({n4, SimNode::GND});
    s1->connect({n1, n2});
    


      // Define system topology
  SystemTopology system(50, 
    SystemNodeList{n1, n2, n3, n4},
    SystemComponentList{r1, v1, r2, l1, non_linear_inductor, s1});


  // Define simulation scenario
  Real timeStep = 1e-6;
  Real finalTime = 0.4;
  String simName = "EMT_NonLinearInductor" + std::to_string(timeStep);


    // Logger
  auto logger = DataLogger::make(simName);
  logger->logAttribute("v_n2", n2->attribute("v"));
  logger->logAttribute("v_res", r1->attribute("v_intf"));
  logger->logAttribute("v_inrush", non_linear_inductor->attribute("v_intf"));
  logger->logAttribute("v_linearInd", l1->attribute("v_intf"));
  logger->logAttribute("i_res", r1->attribute("i_intf"));
  logger->logAttribute("i_inrush", non_linear_inductor->attribute("i_intf"));
  logger->logAttribute("v_source", v1->attribute("v_intf"));
  logger->logAttribute("i_source", v1->attribute("i_intf"));
  logger->logAttribute("Flux", non_linear_inductor->attribute("Flux"));
  logger->logAttribute("i_r2", r2->attribute("i_intf"));
  logger->logAttribute("v_n1", n1->attribute("v"));
  logger->logAttribute("v_n3", n3->attribute("v"));
  logger->logAttribute("v_n4", n4->attribute("v"));
  logger->logAttribute("v_sw", s1->attribute("v_intf"));
  logger->logAttribute("Inductance", non_linear_inductor->attribute("Inductance"));

  Simulation sim(simName);
  sim.setSystem(system);
  sim.setTimeStep(timeStep);
  sim.addEvent(sEvent);
  //sim.doSystemMatrixRecomputation(true);
  sim.setFinalTime(finalTime);
  sim.setDomain(Domain::EMT);
  //sim.doEigenvalueExtraction(true);
  sim.addLogger(logger);
  sim.run();
    
    







}




int main(int argc, char* argv[]) {

    //testNonLinearCharacteristic();

    //testNonLinearInductor();

    //test();

    //testNonLinearInductor_complete();

    //testInductor();

    testInrush();

    return 0;
}



//  cmake --build . --target EMT_SinglePhase_PieceWiseNonLinearCharacteristic -- -j10