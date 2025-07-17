/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/
#include <DPsim.h>
#include "../Examples.h"

using namespace DPsim;
using namespace CPS;
//using namespace CPS::CIM;

//Examples::Components::SynchronousGeneratorKundur::myMachineParameters
 //syngenKundur;
//Switch to trigger fault at generator terminal
Real SwitchOpen = 1e20;
Real SwitchClosed = 1e-8;


void scenario_B_step_B_EMT(String simName, Real timeStep, Real finalTime, Bool startFaultEvent, Bool endFaultEvent, Real startTimeFault, Real endTimeFault, Bool useVarResSwitch, Real cmdInertia_G1, Real cmdDamping_G1) {
	// ----- POWERFLOW FOR INITIALIZATION -----
	Real timeStepPF = finalTime;
	Real finalTimePF = finalTime+timeStepPF;
	String simNamePF = simName + "_PF";
	Logger::setLogDir("logs/" + simNamePF);

	// Components // EDITED ARMIN
	auto BUS_gas_PF = SimNode<Complex>::make("BUS_gas", PhaseType::Single);
	auto BUS_a_PF = SimNode<Complex>::make("BUS_a", PhaseType::Single);
	auto BUS_b_PF = SimNode<Complex>::make("BUS_b", PhaseType::Single);
	auto BUS_psha_PF = SimNode<Complex>::make("BUS_psha", PhaseType::Single);

	double gen_gas_power = 100e6;
	//Synchronous generator 1 // EDITED ARMIN
	auto GEN_gas_PF = SP::Ph1::SynchronGenerator::make("GEN_gas", Logger::Level::debug);
	// setPointVoltage is defined as the voltage at the transfomer primary side and should be transformed to network side
	GEN_gas_PF->setParameters(gen_gas_power, 10.5e3, 0.30880873930480557e6, 10.5e3, PowerflowBusType::VD);
	GEN_gas_PF->setBaseVoltage(10.5e3);
	
	
	double trafo_R = 0.00376196*(220000.0*220000.0)/gen_gas_power/2.0;
	double trafo_L = 0.1007298*(220000.0*220000.0)/gen_gas_power/(2*3.14159*50)/2.0;


	//Transformer
	auto TR_gas_PF = std::make_shared<SP::Ph1::Transformer>("TR_gas", Logger::Level::debug);
    TR_gas_PF->setParameters(10.5e3/*nomVoltageEnd1*/,  220e3/*nomVoltageEnd2*/, gen_gas_power /*ratedPower*/, 
				(10.5e3/220e3)/*ratioAbs*/, 0 /*ratioPhase*/, /*2*1.9129660649*/ /*3.64157728*/trafo_R*2 /*resistance*/, 
	/*2*0.163042774914*/ trafo_L*2/*0.15518646*2*/ /*inductance*/);
    //Real baseVolt = voltageNode1 >= voltageNode2 ? voltageNode1 : voltageNode2;
    TR_gas_PF->setBaseVoltage(220e3);

	auto TR_psh_PF = std::make_shared<SP::Ph1::Transformer>("TR_psh", Logger::Level::debug);
    TR_psh_PF->setParameters(18e3/*nomVoltageEnd1*/,  220e3/*nomVoltageEnd2*/, 200e6 /*ratedPower*/, 
				(18e3/220e3)/*ratioAbs*/, 0 /*ratioPhase*/, 0.41745*2 /*resistance*/, 
	0.0481260775594524*2 /*inductance*/);
    //Real baseVolt = voltageNode1 >= voltageNode2 ? voltageNode1 : voltageNode2;
    TR_psh_PF->setBaseVoltage(220e3);
	

	// shunt
	auto shunt_SR_bcb_PF = SP::Ph1::Shunt::make("shunt_SR_bcb", Logger::Level::debug);
	shunt_SR_bcb_PF->setParameters(3.0989e-06 /*conduntance*/, -3.0989e-04 /*susceptance*/);
	shunt_SR_bcb_PF->setBaseVoltage(220e3); // 2.0659e-07 - 2.0659e-05i


	

	auto shunt_SR_acb_PF = SP::Ph1::Shunt::make("shunt_SR_acb", Logger::Level::debug);
	shunt_SR_acb_PF->setParameters(3.0989e-06 /*conduntance*/, -3.0989e-04 /*susceptance*/);
	shunt_SR_acb_PF->setBaseVoltage(220e3); // 2.0659e-07 - 2.0659e-05i

	//Line1
	auto cable_PF = SP::Ph1::PiLine::make("cable", Logger::Level::debug);
	cable_PF->setParameters(0.032 * 10.5 /*R/km * km*/, 0.321493e-3 * 10.5 /*L/km * km*/, 0.21915e-6 * 10.5 /*Capacitance*/, 1e-15);
	cable_PF->setBaseVoltage(220e3);
	
	auto line_3_PF = SP::Ph1::PiLine::make("line_3", Logger::Level::debug);
	line_3_PF->setParameters(0.0749 * 5, 1.270693e-3 * 5, 0.00466961e-6 * 5, 1e-15);
	line_3_PF->setBaseVoltage(220e3);





	// Topology
	GEN_gas_PF->connect({ BUS_gas_PF });
	TR_gas_PF->connect({ BUS_gas_PF, BUS_b_PF});
	cable_PF->connect({ BUS_b_PF, BUS_a_PF});
	shunt_SR_bcb_PF->connect({ BUS_b_PF });
	shunt_SR_acb_PF->connect({ BUS_a_PF });
	line_3_PF->connect({ BUS_a_PF, BUS_psha_PF});




	// small power flow
	auto systemPF = SystemTopology(50, // das ist freq??
			SystemNodeList{BUS_gas_PF, BUS_b_PF, BUS_a_PF, BUS_psha_PF},
			SystemComponentList{GEN_gas_PF, TR_gas_PF, shunt_SR_bcb_PF, shunt_SR_acb_PF, cable_PF, line_3_PF});


	// Logging
	auto loggerPF = DataLogger::make(simNamePF);
	loggerPF->logAttribute("V_BUS_b_PF", BUS_b_PF->attribute("v"));
	loggerPF->logAttribute("V_BUS_gas_PF", BUS_gas_PF->attribute("v"));

	// Simulation
	Simulation simPF(simNamePF, Logger::Level::debug);
	simPF.setSystem(systemPF);
	simPF.setTimeStep(timeStepPF);
	simPF.setFinalTime(finalTimePF);
	simPF.setDomain(Domain::SP);
	simPF.setSolverType(Solver::Type::NRP);
	simPF.setSolverAndComponentBehaviour(Solver::Behaviour::Initialization);
	simPF.doInitFromNodesAndTerminals(true);
	simPF.addLogger(loggerPF);
	simPF.run();


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// ----- Dynamic simulation ------
	String simNameEMT = simName + "_EMT";
	Logger::setLogDir("logs/"+simNameEMT);

	// Nodes
	auto BUS_gas_EMT = SimNode<Real>::make("BUS_gas", PhaseType::ABC);
	auto BUS_a_EMT = SimNode<Real>::make("BUS_a", PhaseType::ABC);
	auto BUS_b_EMT = SimNode<Real>::make("BUS_b", PhaseType::ABC);
	auto BUS_psha_EMT = SimNode<Real>::make("BUS_psha", PhaseType::ABC);
	auto BUS_tr_psh_EMT = SimNode<Real>::make("BUS_tr_psh", PhaseType::ABC);
	auto BUS_psh_EMT = SimNode<Real>::make("BUS_psh", PhaseType::ABC);

	auto Bus_nl_ind_EMT = SimNode<Real>::make("Bus_nl_ind", PhaseType::ABC);
	auto Bus_resistor_EMT = SimNode<Real>::make("Bus_resistor", PhaseType::ABC);
	


	// Components

	// generator 6 order - working myb - odavde
	// auto GEN_gas_EMT =
    //   CPS::EMT::Ph3::SynchronGeneratorVBR::make("GEN_gas", Logger::Level::debug);
 	// 	 GEN_gas_EMT->setBaseAndOperationalPerUnitParameters(
    //   gen_gas_power/*nomPower*/, 10.5e3/*nomVoltage*/, 50/*nomFreq*/,
    //   2/*poleNum*/, 1300/*nomFieldCurr*/, 0.002/*Rs*/,
    //   2.4/*Ld*/, 1.33 /*Lq*/, 0.31/*Ld_t*/, 1.2/*Lq_t*/,
    //   0.24/*Ld_s*/, 0.35/*Lq_s*/, 0.135/*Ll*/, 1.45/*Td0_t*/,
    //   0.000001/*Tq0_t*/, 0.022/*Td0_s*/, 0.0095/*Tq0_s*/,
    //   5/*H*/);



	// // ovaj governor radi koliko-toliko
	// GEN_gas_EMT->addGovernor(0.3 /*Real Ta 0.3*/, 12.0 /*Real Tb 7.0*/, 0.5 /*Real Tc 0.5*/, 0.3 /*0.33 Real Fa 0.3*/,
    //                  0.3 /*0.065 Real Fb 0.3*/, 0.4 /*0.055 Real Fc 0.4*/, 20 /*Real K 20*/, 0.1 /*Real Tsr 0.1*/,
    //                  0.3 /*Real Tsm 0.3*/, 0.30880873930480557e6 / gen_gas_power /*Real Tm_init*/,
    //                  0.30880873930480557e6 / gen_gas_power /*Real PmRef*/);
	// 				//  (Real Ta, Real Tb, Real Tc, Real Fa,
    //                 //                       Real Fb, Real Fc, Real K, Real Tsr,
    //                 //                       Real Tsm, Real Tm_init, Real PmRef)



	// // second exciter parameters
	// GEN_gas_EMT->addExciter(0.01 /*Real Ta*/, 500 /*Real Ka*/, 0.46 /*Real Te*/, 0.0435 /*Real Ke*/,
    //                      1 /*Real Tf*/, 0.1 /*Real Kf*/, 0.01 /*Real Tr*/, 7.65 /*MaxVr*/, -6.5 /*MinVr*/);

	auto GEN_gas_EMT =
      CPS::EMT::Ph3::SynchronGenerator4OrderVBR::make("GEN_gas", Logger::Level::debug);
	
	GEN_gas_EMT->setOperationalParametersPerUnit(gen_gas_power,10.5e3,
                                       50, 5, 2.4, 1.33,
                                       0.135, 0.31, 1.2,
                                       1.45, 0.000001);

	

	// Trafo
	auto TR_gas_EMT = EMT::Ph3::Transformer::make("TR_gas", "TR_gas", Logger::Level::debug, true);
	TR_gas_EMT->setParameters(10.5e3, 220e3, gen_gas_power,
                     (10.5e3/220e3), 0, Math::singlePhaseParameterToThreePhase(trafo_R*2),
                     Math::singlePhaseParameterToThreePhase(trafo_L*2));

	auto TR_psh_EMT = EMT::Ph3::Transformer::make("TR_psh", "TR_psh", Logger::Level::debug, true);
	TR_psh_EMT->setParameters(18.0e3, 220e3, 200e6,
                     (18.0e3/220e3), 0, Math::singlePhaseParameterToThreePhase(0.41745*2),
                     Math::singlePhaseParameterToThreePhase(0.0481260775594524*2));
	 

	// Load
	auto load_shunt_bus_b_EMT = EMT::Ph3::RXLoad::make("shunt_SR_bcb", Logger::Level::debug);
	load_shunt_bus_b_EMT->setParameters(CPS::Math::singlePhasePowerToThreePhase((1.0*220000)*(1.0*220000)/(322698.9333)), 
	CPS::Math::singlePhasePowerToThreePhase((1.0*220000)*(1.0*220000)/(3226.9893)), 220e3);

	auto load_shunt_bus_a_EMT = EMT::Ph3::RXLoad::make("shunt_SR_acb", Logger::Level::debug);
	load_shunt_bus_a_EMT->setParameters(CPS::Math::singlePhasePowerToThreePhase((1.0*220000)*(1.0*220000)/(322698.9333)), 
	CPS::Math::singlePhasePowerToThreePhase((1.0*220000)*(1.0*220000)/(3226.9893)), 220e3);


	// Cable
	auto cable_EMT = EMT::Ph3::PiLine::make("cable", Logger::Level::debug);
	cable_EMT->setParameters(Math::singlePhaseParameterToThreePhase(0.032 * 10.5), 
	                      Math::singlePhaseParameterToThreePhase(0.321493e-3 * 10.5), 
					      Math::singlePhaseParameterToThreePhase(0.21915e-6 * 10.5),
						  Math::singlePhaseParameterToThreePhase(1e-15));


	// Line3
	auto line3_EMT = EMT::Ph3::PiLine::make("line_3", Logger::Level::debug);
	line3_EMT->setParameters(Math::singlePhaseParameterToThreePhase(0.0749 * 5), 
	                      Math::singlePhaseParameterToThreePhase(1.270693e-3 * 5), 
					      Math::singlePhaseParameterToThreePhase(0.00466961e-6 * 5),
						  Math::singlePhaseParameterToThreePhase(1e-15));


	auto Inductor = EMT::Ph3::Inductor::make("Inductor", Logger::Level::debug);
	Inductor->setParameters(Math::singlePhaseParameterToThreePhase(0.0001));

	auto Resistor = EMT::Ph3::Resistor::make("Resistor", Logger::Level::debug);
	Resistor->setParameters(Math::singlePhaseParameterToThreePhase(0.0001));
	

	// Switch

	auto breaker_psh_EMT = CPS::EMT::Ph3::Switch::make("breaker_psh", Logger::Level::debug);
	breaker_psh_EMT->setParameters(Math::singlePhaseParameterToThreePhase(SwitchOpen), 
							Math::singlePhaseParameterToThreePhase(SwitchClosed));
	breaker_psh_EMT->openSwitch();
	
	
	// Nonlinear inductor

Real V_LL_RMS = 220e3;

Real V_LG_Peak = sqrt(2.0/3.0) * V_LL_RMS;

Real S_three_phase = 200e6;

Real I_one_phase_RMS = S_three_phase / (V_LL_RMS * sqrt(3.0));

Real I_one_phase_Peak = sqrt(2.0) * I_one_phase_RMS;    // base value for characteristics - current

Real flux_LG_Peak = V_LG_Peak / (2 * M_PI * 50);        // base value for characteristics - flux

std::vector<CPS::EMT::Ph1::PieceWiseNonLinearCharacteristic::Point> Pts = {

{ 0.0 * I_one_phase_Peak,    0.0 * flux_LG_Peak  },

{ 1.1e-3 * I_one_phase_Peak, 1.1 * flux_LG_Peak },

};

auto characteristic = std::make_shared<

CPS::EMT::Ph1::PieceWiseNonLinearCharacteristic

>(Pts, 0.2 * flux_LG_Peak/I_one_phase_Peak);

	auto nl_ind = CPS::EMT::Ph3::NonLinearInductor::make("nl_ind", Logger::Level::info);
	nl_ind->setPieceWiseCharacteristic(characteristic);

	// Set time step for the nonlinear inductor
	nl_ind->setTimeStep(timeStep);
	

	// Topology
	GEN_gas_EMT->connect({BUS_gas_EMT});
	TR_gas_EMT->connect({BUS_gas_EMT, BUS_b_EMT});
	cable_EMT->connect({BUS_b_EMT, BUS_a_EMT});
	line3_EMT->connect({BUS_a_EMT, BUS_psha_EMT});
	load_shunt_bus_b_EMT->connect({BUS_b_EMT});
	load_shunt_bus_a_EMT->connect({BUS_a_EMT});
	
	TR_psh_EMT->connect({ BUS_psh_EMT, BUS_tr_psh_EMT});
	breaker_psh_EMT->connect({ BUS_psha_EMT, BUS_tr_psh_EMT});

	nl_ind->connect({ Bus_nl_ind_EMT, EMT::SimNode::GND });

	Inductor->connect({BUS_tr_psh_EMT, Bus_resistor_EMT});

	Resistor->connect({Bus_resistor_EMT, Bus_nl_ind_EMT});

	auto systemEMT = SystemTopology(50,
			SystemNodeList{BUS_gas_EMT, BUS_b_EMT, BUS_a_EMT, BUS_psha_EMT,
							BUS_psh_EMT, BUS_tr_psh_EMT, Bus_nl_ind_EMT, Bus_resistor_EMT},
			SystemComponentList{GEN_gas_EMT, TR_gas_EMT, cable_EMT, load_shunt_bus_b_EMT, 
								load_shunt_bus_a_EMT, line3_EMT,
								TR_psh_EMT, breaker_psh_EMT,
								nl_ind, Inductor, Resistor});

	// Initialization of dynamic topology
	systemEMT.initWithPowerflow(systemPF, Domain::EMT);

	for (auto attr : breaker_psh_EMT->attributes()) {
    	std::string name = attr.first;
    	std::cout << name << std::endl;
  	}


	// Logging
	auto loggerEMT = DataLogger::make(simNameEMT);
	loggerEMT->logAttribute("BUS_gas_EMT_v", BUS_gas_EMT->attribute("v"));
	loggerEMT->logAttribute("BUS_b_EMT_v", BUS_b_EMT->attribute("v"));
	loggerEMT->logAttribute("BUS_a_EMT_v", BUS_a_EMT->attribute("v"));
	loggerEMT->logAttribute("BUS_psha_EMT_v", BUS_psha_EMT->attribute("v"));
	loggerEMT->logAttribute("BUS_psh_EMT_v", BUS_psh_EMT->attribute("v"));
	loggerEMT->logAttribute("BUS_tr_psh_EMT_v", BUS_tr_psh_EMT->attribute("v"));

	loggerEMT->logAttribute("GEN_gas_EMT_wr", GEN_gas_EMT->attribute("w_r"));
	
	loggerEMT->logAttribute("TR_gas_EMT_i", TR_gas_EMT->attribute("i_intf"));
	loggerEMT->logAttribute("TR_psh_EMT_i", TR_psh_EMT->attribute("i_intf"));

	loggerEMT->logAttribute("breaker_psh_EMT_i", breaker_psh_EMT->attribute("i_intf"));

	loggerEMT->logAttribute("inrush", nl_ind->attribute("i_intf"));
	loggerEMT->logAttribute("inrush_v", nl_ind->attribute("v_intf"));

	Simulation simEMT(simNameEMT, Logger::Level::debug);
  	simEMT.doInitFromNodesAndTerminals(true);
  	simEMT.setSystem(systemEMT);
  	simEMT.setTimeStep(timeStep);
  	simEMT.setFinalTime(finalTime);
  	simEMT.setDomain(Domain::EMT);
  	simEMT.addLogger(loggerEMT);
	simEMT.doSystemMatrixRecomputation(true);

	if (useVarResSwitch == true) {
		simEMT.doSystemMatrixRecomputation(true);
	}

	if (startFaultEvent){
		auto sw1 = SwitchEvent3Ph::make(startTimeFault, breaker_psh_EMT, true);
		simEMT.addEvent(sw1);
	}

	simEMT.run();

}



int main(int argc, char* argv[]) {	
		
/*
	//Simultion parameters
	String simName="scenario_B_step_A";
	Real finalTime = 30;
	Real timeStep = 0.001;
	
	scenario_B_step_A(simName, timeStep, finalTime);
*/

//Simultion parameters
	String simName="scenario_B_step_B";
	Real finalTime = 1.00;
	Real timeStep = 1e-6;
	Bool startFaultEvent=true;
	Bool endFaultEvent=true;
	Bool useVarResSwitch=false;
	Real startTimeFault=0.4;
	Real endTimeFault=11100;
	Real cmdInertia_G1= 1.0;
	Real cmdDamping_G1= 1.0;


	CommandLineArgs args(argc, argv);
	if (argc > 1) {
		timeStep = args.timeStep;
		finalTime = args.duration;
		if (args.name != "dpsim")
			simName = args.name;
		if (args.options.find("SCALEINERTIA_G1") != args.options.end())
			cmdInertia_G1 = args.getOptionReal("SCALEINERTIA_G1");
		if (args.options.find("SCALEDAMPING_G1") != args.options.end())
			cmdDamping_G1 = args.getOptionReal("SCALEDAMPING_G1");
		if (args.options.find("STARTTIMEFAULT") != args.options.end())
			startTimeFault = args.getOptionReal("STARTTIMEFAULT");
		if (args.options.find("ENDTIMEFAULT") != args.options.end())
			endTimeFault = args.getOptionReal("ENDTIMEFAULT");
		// if (args.options.find("USEVARRESSWITCH") != args.options.end())
		// 	useVarResSwitch = args.options["USEVARRESSWITCH"];	
		// if (args.options.find("FAULTRESISTANCE") != args.options.end())
		// 	SwitchClosed = args.options["FAULTRESISTANCE"];	
	}
	
	scenario_B_step_B_EMT(simName, timeStep, finalTime, startFaultEvent, endFaultEvent, startTimeFault, endTimeFault, useVarResSwitch, cmdInertia_G1, cmdDamping_G1);

}


//  cmake --build . --target EMT_NonLinearInductorTest -- -j20