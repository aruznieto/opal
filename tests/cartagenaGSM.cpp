/***************************************************************/
//
//Copyright (c) 2021 Esteban Egea-Lopez http://girtel.upct.es/~simulations
//
/**************************************************************/
#include "cartagenaGSM.h"
#include "../timer.h"
#include "../Opal.h"
#include <memory>
#include <string>
#include <fstream>
#include <random>
#include "../curvedMeshSimulation.h"
#include "../curvedFlatMeshSimulation.h"
#include "../basicSimulation.h"
#include "../flatSimulation.h"
#include "../singleDiffraction.h"
#include "../rayDensityNormalizationSimulation.h"
#include "../util.h"

using namespace opal;
using namespace optix;
CartagenaGSM::CartagenaGSM(OpalSceneManager*   sceneManager, float sphereRadius, bool useDepolarization) : BasicTests(sceneManager, sphereRadius, useDepolarization) {
}
void  CartagenaGSM::loadScenarioCartagena(bool addRandom, bool useRDN) {
	//A few receivers
	float3 baseStation=make_float3(1153.33f,25.0f,936.0f);
	Timer timer;
	Timer timerRT;
	float freq = 1800e6f;
	std::cout<<"Load scenario Cartagena GSM"<<std::endl;
	timer.start();	
	//Init context before doing anything else
	//sceneManager->enableGenerateRaysOnLaunch();
	sceneManager->setMinEpsilon(1e-3);
	sceneManager->setUseAntennaGain(false);
	sceneManager->enableGenerateRaysOnLaunch();	
	ComputeMode mode=ComputeMode::VOLTAGE;
	if (useRDN) {
		RayDensityNormalizationSimulation* sim = new RayDensityNormalizationSimulation(sceneManager);
		sceneManager->setSimulation(sim);
	}  else {
		if (useDepolarization) {
			LPFlatMeshReflectionSimulation* sim = new LPFlatMeshReflectionSimulation(sceneManager);
			sceneManager->setSimulation(sim);
			sim->setComputeMode(mode);
			//sim->setEnableTraceLog(false);
			//sim->setEnableSimulation(false);
		} else {
			BasicFlatMeshReflectionSimulation* sim = new BasicFlatMeshReflectionSimulation(sceneManager);
			sceneManager->setSimulation(sim);
			sim->setComputeMode(mode);
		}
	}
//Add diffraction
	SingleDiffraction* simd= new SingleDiffraction(sceneManager);
	sceneManager->setSimulation(simd);
	simd->setComputeMode(mode);
	//simd->setEnableTraceLog(true);
	//simd->setEnableSimulation(false);
	
	sceneManager->initContext(freq);
//Exceptions
	//sceneManager->enableExceptions();	

//Load files here
	ScenarioLoader* loader=new ScenarioLoader(sceneManager);
        //std::string path("lora/cartagena.json");
        //std::string path("cartagena2.json");
        //std::string path("eldi.json");
        std::string path("cartagena-tfg.json");
	loader->loadJSONScenario(path);
        //std::string path("meshes/cartagena");
	//loader->loadMeshesFromFiles(path);
	//loader->loadEdgesFromFiles(path);

	optix::float3 polarization = make_float3(0.0f, 1.0f, 0.0f); 
//	optix::float3 posrx = make_float3(1547.14f, 40.69f, 620.8f);
//	float3 posrx = make_float3(1501.6, 25.1f, 609.9f);
	//float3 posrx=rx_e;
	//float3 posrx=rx3;
	std::string rx_file("rxn.txt");
	std::vector<float2> rx=loadReceiversFromFile(rx_file);
	std::default_random_engine gen;
  	std::uniform_real_distribution<float> df(-0.5f,0.5f);
	int i=0;
	for (auto p : rx) {
		float3 posrx = make_float3(p.x, 1.7f, p.y);
		if (addRandom) {
			posrx.x=posrx.x+df(gen);
			posrx.y=posrx.y+df(gen);
			posrx.z=posrx.z+df(gen);

		}
		sceneManager->addReceiver(i, posrx,polarization, sphereRadius, sceneManager->printPower);
		++i;
	}
	//sceneManager->addReceiver(0, posrx,polarization, sphereRadius, sceneManager->printPower);
	
	//AntennaGain gains=sceneManager->loadGainsFromFileIndBPower("lora/gananciasAntena.txt");
	//int gainId=sceneManager->registerAntennaGain(gains);
	//sceneManager->registerReceiverGain(0,gainId);
	//***Single ray transmit****
	//float3 mRay=normalize(make_float3(-0.004128737841, -0.9902680516, 0.1391119212));
	//sceneManager->createRaySphere2D(1,1,&mRay);
	

	if (useRDN) {
		sceneManager->finishSceneContext();
		int rayD=10000;
		sceneManager->setRayRange(0.0,180.0,0.0,360.0,rayD,rayD);
		RayDensityNormalizationSimulation* s=dynamic_cast<RayDensityNormalizationSimulation*>(sceneManager->getSimulation(0));
		s->setInitialDensity(((float)sceneManager->getRaySphere().rayCount)/(4*M_PIf));
		s->setFiltering(2u);

	} else {
		sceneManager->createRaySphere2D(0.0f,0.1,180.0f,0.0f,0.1,360.0f);
		sceneManager->finishSceneContext();
	}
//	//float3 postx = make_float3(1501.6, 25.1f, 609.9f);
//	//optix::float3 postx = make_float3(1547.14f, 40.69f, 620.8f);
//	optix::float3 postx = tx_f;
	optix::float3 postx = baseStation;
	timerRT.start();
	sceneManager->transmit(i+1, 0.01, postx, polarization);
	timer.stop();
	timerRT.stop();
	std::cout<<"Time\t"<<timer.getTime()<<"\t"<<timerRT.getTime()<<std::endl;
	delete loader;
}
std::vector<float2> CartagenaGSM::loadReceiversFromFile(std::string file) {
	ScenarioLoader* sl=new ScenarioLoader(sceneManager);
        std::ifstream infile(file);
	if (!infile.good()) {
		    std::cout<<"Error opening "<<file<<std::endl;
		    throw  opal::Exception("loadReceiversFromFile(): error opening file");
	}
    	std::cout<<"Loading receivers from  "<<file<<std::endl;
	std::string line;
	std::vector<float2> rx;
	while (std::getline(infile, line) ){
		optix::float2 v=sl->readFloat2(line);
		rx.push_back(v);
	}
	infile.close();
	delete sl;
	return rx;
}
void CartagenaGSM::runTests(std::string test) {
	std::vector<int> tokens=parseTestString(test);
	
	loadScenarioCartagena(tokens[0],tokens[1]);
}
