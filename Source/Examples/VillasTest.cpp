#include "VillasTest.h"
#include "../Simulation.h"
#include "../VillasInterface.h"
#include "../Utilities.h"

using namespace DPsim;

void DPsim::villasExample()
{
	// Very simple test circuit. Just 2 resistors and a current read from VILLASnode.
	Logger log, llog, rlog;
	std::vector<BaseComponent*> comps;

	ExternalCurrentSource *ecs = new ExternalCurrentSource("i1", 1, 0);
	comps.push_back(ecs);
	comps.push_back(new LinearResistor("r1", 1, 2, 1));
	LinearResistor *r2 = new LinearResistor("r2", 2, 0, 1);
	comps.push_back(r2);
	VillasInterface *villas = new VillasInterface("/villas1");
	villas->registerCurrentSource(ecs, 0, 1);
	villas->registerExportedVoltage(1, 0, 0, 1);
	villas->registerExportedCurrent(r2, 2, 3);

	// Set up simulation
	Real timeStep = 0.01;
	Simulation newSim(comps, 2.0*M_PI*50.0, timeStep, 0.1, log);
	newSim.addExternalInterface(villas);

	// Main Simulation Loop
	std::cout << "Start simulation." << std::endl;
	while (newSim.step(log, llog, rlog))
	{
		newSim.increaseByTimeStep();
		updateProgressBar(newSim.getTime(), newSim.getFinalTime());
	}
	std::cout << "Simulation finished." << std::endl;
	log.WriteLogToFile("output.log");
	rlog.WriteLogToFile("rvector.log");
	llog.WriteLogToFile("lvector.log");
	for (auto comp : comps) {
		delete comp;
	}
	delete villas;
}
