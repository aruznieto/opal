/***************************************************************/
//
//Copyright (c) 2021 Esteban Egea-Lopez http://girtel.upct.es/~simulations
//
/**************************************************************/
#include "parking.h"
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
Parking::Parking(OpalSceneManager*   sceneManager, float sphereRadius, bool useDepolarization) : BasicTests(sceneManager, sphereRadius, useDepolarization) {
	this->freq=2e9;
	this->baseStation=make_float3(1543.56f,37.6f,620.39f);
	this->gsufix=std::to_string(0);
	this->replications=1;
	this->useGain=true;
}
void  Parking::loadScenarioParkingMultifrequency(bool addRandom, bool useRDN, std::string path, std::vector<float>& frequencies, std::vector<float>& radius, std::string outputFile) {
	//A few receivers
	//float3 baseStation=make_float3(1543.56f,37.6f,620.39f);
	Timer timer;
	Timer timerRT;
	//float freq = 1800e6f;
	std::cout<<"Load parking with scenario "<<path<<"; freq="<<freq<<";tx="<<baseStation<<std::endl;
	timer.start();	
	//Init context before doing anything else
	//sceneManager->enableGenerateRaysOnLaunch();
	sceneManager->enableMultiChannel();
	sceneManager->setMinEpsilon(1e-3);
	sceneManager->setUseAntennaGain(useGain);
	sceneManager->enableGenerateRaysOnLaunch();
   	sceneManager->setPrintRecords(true);	
	ComputeMode mode=ComputeMode::VOLTAGE;
	if (useRDN) {
		RayDensityNormalizationSimulation* sim = new RayDensityNormalizationSimulation(sceneManager);
		sceneManager->setSimulation(sim);
		sim->setComputeMode(mode);
	}  else {
		if (useDepolarization) {
			LPFlatMeshReflectionSimulation* sim = new LPFlatMeshReflectionSimulation(sceneManager);
			sceneManager->setSimulation(sim);
			sim->setComputeMode(mode);
			sim->setEnableTraceLog(false);
			//sim->setEnableSimulation(false);
			//sim->setPrintHits(true);
		} else {
			BasicFlatMeshReflectionSimulation* sim = new BasicFlatMeshReflectionSimulation(sceneManager);
			sceneManager->setSimulation(sim);
			sim->setComputeMode(mode);
			//sim->setPrintHits(true);
		}
	}
//Add diffraction
	SingleDiffraction* simd= new SingleDiffraction(sceneManager);
	sceneManager->setSimulation(simd);
	simd->setComputeMode(mode);
	simd->setEnableTraceLog(false);
	//simd->setPrintHits(true);
	//simd->setEnableSimulation(false);
	
	sceneManager->initContext(freq);
//Exceptions
	//sceneManager->enableExceptions();	
	int gainIdTx=-1;
	int gainIdRx=-1;
	if (useGain) {	
		std::string gp("parking/april23/gainT");
		std::string gainPathT=gp+gsufix+".txt";
		AntennaGain gains=sceneManager->loadGainsFromFileIndBPower(gainPathT.c_str());
		gainIdTx=sceneManager->registerAntennaGain(gains);
		//Receiver gain
		std::string gpR("parking/april23/gainR");
		std::string gainPathR=gpR+gsufix+".txt";
		AntennaGain gainsRx=sceneManager->loadGainsFromFileIndBPower(gainPathR.c_str());
		gainIdRx=sceneManager->registerAntennaGain(gainsRx);
	}
//Load files here
	ScenarioLoader* loader=new ScenarioLoader(sceneManager);
        //std::string path("lora/cartagena.json");
        //std::string path("cartagena2.json");
        //std::string path("eldi.json");
        //std::string path("cartagena-tfg.json");
	loader->loadJSONScenario(path);
        //std::string path("meshes/cartagena");
	//loader->loadMeshesFromFiles(path);
	//loader->loadEdgesFromFiles(path);

	optix::float3 polarization = make_float3(0.0f, 1.0f, 0.0f); 


	std::string rx_file("parking/rx-ap23.txt");
	std::vector<float3> rx=loadReceiversFromFile(rx_file);
	std::default_random_engine gen;
  	std::uniform_real_distribution<float> df(-0.5f,0.5f);
  	std::uniform_real_distribution<float> dfy(-0.2f,0.2f);
	int i=0;
	//float3 posrx=make_float3(1554.8,22.9854,608.7);
	//float3 posrx=make_float3(1619.5,22.85484,580.6);
	//sceneManager->addReceiver(0, posrx,polarization, sphereRadius, sceneManager->printPower);
	//sceneManager->registerReceiverGain(0,gainIdRx);
        std::vector<float3> registeredReceivers;
	for (auto p : rx) {
		if (addRandom) {
			for (int j=0; j<replications; j++) {
				float3 posrx = make_float3(p.x,  p.y, p.z);
				posrx.x=posrx.x+df(gen);
				posrx.y=posrx.y+dfy(gen);
				posrx.z=posrx.z+df(gen);

				sceneManager->addReceiver(i, posrx,polarization, sphereRadius, sceneManager->printPower);
				if (useGain) {
					sceneManager->registerReceiverGain(i,gainIdRx);
				}
 				registeredReceivers.push_back(posrx);
				++i;
			}
		} else {
			float3 posrx = make_float3(p.x,  p.y, p.z);
			sceneManager->addReceiver(i, posrx,polarization, sphereRadius, sceneManager->printPower);
			if (useGain) {
				sceneManager->registerReceiverGain(i,gainIdRx);
			}
 			registeredReceivers.push_back(posrx);
			++i;
		}
	}
	//sceneManager->addReceiver(0, posrx,polarization, sphereRadius, sceneManager->printPower);
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
		sceneManager->finishSceneContext();
		sceneManager->createRaySphere2D(0.0f,0.1,180.0f,0.0f,0.1,360.0f);
	}
//	//float3 postx = make_float3(1501.6, 25.1f, 609.9f);
	//optix::float3 postx = make_float3(1547.14f, 40.69f, 620.8f);
	//optix::float3 postx = make_float3(1555.58f, 22.61f, 616.4f);
//	optix::float3 postx = tx_f;
	optix::float3 postx = baseStation;
	if (useGain) {
		sceneManager->registerTransmitterGain(i+1,gainIdTx);
	}
	ResultReport report;
	for (auto f: frequencies) {
	//for (int f=0; f<1;++f) {
		for (auto rad: radius) {
			for (int k=0; k<registeredReceivers.size(); ++k) {
				sceneManager->updateReceiver(k,registeredReceivers[k],rad);
			}
			timerRT.start();
			std::cout<<"Launching f="<<f<<"; rad="<<rad<<std::endl;
			ResultReport* rp=sceneManager->transmit(i+1, 1, postx, polarization, f, false, false);
			report.merge(*rp);
			timerRT.stop();
			std::cout<<"Time launch\t"<<timerRT.getTime()<<std::endl;
		}
	}
	timer.stop();
	std::cout<<"Time total\t"<<timer.getTime()<<std::endl;
	std::string rpath=outputFile+".csv";
	report.toCSV(rpath);
	delete loader;
}
void  Parking::loadScenarioParking(bool addRandom, bool useRDN, std::string path) {
	//A few receivers
	//float3 baseStation=make_float3(1543.56f,37.6f,620.39f);
	Timer timer;
	Timer timerRT;
	//float freq = 1800e6f;
	std::cout<<"Load parking with scenario "<<path<<"; freq="<<freq<<";tx="<<baseStation<<std::endl;
	timer.start();	
	//Init context before doing anything else
	//sceneManager->enableGenerateRaysOnLaunch();
	sceneManager->setMinEpsilon(1e-3);
	sceneManager->setUseAntennaGain(useGain);
	sceneManager->enableGenerateRaysOnLaunch();
   	sceneManager->setPrintRecords(true);	
	ComputeMode mode=ComputeMode::VOLTAGE;
	if (useRDN) {
		RayDensityNormalizationSimulation* sim = new RayDensityNormalizationSimulation(sceneManager);
		sceneManager->setSimulation(sim);
	}  else {
		if (useDepolarization) {
			LPFlatMeshReflectionSimulation* sim = new LPFlatMeshReflectionSimulation(sceneManager);
			sceneManager->setSimulation(sim);
			sim->setComputeMode(mode);
			sim->setEnableTraceLog(false);
			//sim->setEnableSimulation(false);
			//sim->setPrintHits(true);
		} else {
			BasicFlatMeshReflectionSimulation* sim = new BasicFlatMeshReflectionSimulation(sceneManager);
			sceneManager->setSimulation(sim);
			sim->setComputeMode(mode);
			//sim->setPrintHits(true);
		}
	}
//Add diffraction
	SingleDiffraction* simd= new SingleDiffraction(sceneManager);
	sceneManager->setSimulation(simd);
	simd->setComputeMode(mode);
	simd->setEnableTraceLog(false);
	//simd->setPrintHits(true);
	//simd->setEnableSimulation(false);
	
	sceneManager->initContext(freq);
//Exceptions
	//sceneManager->enableExceptions();	
	int gainIdTx=-1;
	int gainIdRx=-1;
	if (useGain) {	
		std::string gp("parking/april23/gainT");
		std::string gainPathT=gp+gsufix+".txt";
		AntennaGain gains=sceneManager->loadGainsFromFileIndBPower(gainPathT.c_str());
		gainIdTx=sceneManager->registerAntennaGain(gains);
		//Receiver gain
		std::string gpR("parking/april23/gainR");
		std::string gainPathR=gpR+gsufix+".txt";
		AntennaGain gainsRx=sceneManager->loadGainsFromFileIndBPower(gainPathR.c_str());
		gainIdRx=sceneManager->registerAntennaGain(gainsRx);
	}
//Load files here
	ScenarioLoader* loader=new ScenarioLoader(sceneManager);
        //std::string path("lora/cartagena.json");
        //std::string path("cartagena2.json");
        //std::string path("eldi.json");
        //std::string path("cartagena-tfg.json");
	loader->loadJSONScenario(path);
        //std::string path("meshes/cartagena");
	//loader->loadMeshesFromFiles(path);
	//loader->loadEdgesFromFiles(path);

	optix::float3 polarization = make_float3(0.0f, 1.0f, 0.0f); 


	std::string rx_file("parking/rx-vehicles.txt");
	std::vector<float3> rx=loadReceiversFromFile(rx_file);
	std::default_random_engine gen;
  	std::uniform_real_distribution<float> df(-0.5f,0.5f);
  	std::uniform_real_distribution<float> dfy(-0.2f,0.2f);
	int i=0;
	//float3 posrx=make_float3(1554.8,22.9854,608.7);
	//float3 posrx=make_float3(1619.5,22.85484,580.6);
	//sceneManager->addReceiver(0, posrx,polarization, sphereRadius, sceneManager->printPower);
	//sceneManager->registerReceiverGain(0,gainIdRx);
	for (auto p : rx) {
		if (addRandom) {
			for (int j=0; j<replications; j++) {
				float3 posrx = make_float3(p.x,  p.y, p.z);
				posrx.x=posrx.x+df(gen);
				posrx.y=posrx.y+dfy(gen);
				posrx.z=posrx.z+df(gen);

				sceneManager->addReceiver(i, posrx,polarization, sphereRadius, sceneManager->printPower);
				if (useGain) {
					sceneManager->registerReceiverGain(i,gainIdRx);
				}
				++i;
			}
		} else {
			float3 posrx = make_float3(p.x,  p.y, p.z);
			sceneManager->addReceiver(i, posrx,polarization, sphereRadius, sceneManager->printPower);
			if (useGain) {
				sceneManager->registerReceiverGain(i,gainIdRx);
			}
			++i;
		}
	}
	//sceneManager->addReceiver(0, posrx,polarization, sphereRadius, sceneManager->printPower);
	//sceneManager->registerReceiverGain(0,gainId);
	//***Single ray transmit****
	//float3 mRay=normalize(make_float3(-0.004128737841, -0.9902680516, 0.1391119212));
	//sceneManager->createRaySphere2D(1,1,&mRay);
	

	if (useRDN) {
		sceneManager->finishSceneContext();
		int rayD=30000;
		sceneManager->setRayRange(0.0,180.0,0.0,360.0,rayD,rayD);
		RayDensityNormalizationSimulation* s=dynamic_cast<RayDensityNormalizationSimulation*>(sceneManager->getSimulation(0));
		s->setInitialDensity(((float)sceneManager->getRaySphere().rayCount)/(4*M_PIf));
		s->setFiltering(2u);

	} else {
		sceneManager->finishSceneContext();
		sceneManager->createRaySphere2D(0.0f,0.05,180.0f,0.0f,0.05,360.0f);
	}
//	//float3 postx = make_float3(1501.6, 25.1f, 609.9f);
	//optix::float3 postx = make_float3(1547.14f, 40.69f, 620.8f);
	//optix::float3 postx = make_float3(1555.58f, 22.61f, 616.4f);
//	optix::float3 postx = tx_f;
	optix::float3 postx = baseStation;
	timerRT.start();
	if (useGain) {
		sceneManager->registerTransmitterGain(i+1,gainIdTx);
	}
	sceneManager->transmit(i+1, 1, postx, polarization);
	timer.stop();
	timerRT.stop();
	std::cout<<"Time\t"<<timer.getTime()<<"\t"<<timerRT.getTime()<<std::endl;
	delete loader;
}
std::vector<float3> Parking::loadReceiversFromFile(std::string file) {
	ScenarioLoader* sl=new ScenarioLoader(sceneManager);
        std::ifstream infile(file);
	if (!infile.good()) {
		    std::cout<<"Error opening "<<file<<std::endl;
		    throw  opal::Exception("loadReceiversFromFile(): error opening file");
	}
    	std::cout<<"Loading receivers from  "<<file<<std::endl;
	std::string line;
	std::vector<float3> rx;
	while (std::getline(infile, line) ){
		optix::float3 v=sl->readFloat3(line);
		rx.push_back(v);
	}
	infile.close();
	delete sl;
	return rx;
}
void Parking::runTests(std::string test, std::string outputFile) {

	std::vector<int> tokens=parseTestString(test);

	this->replications=20;
	bool useRDN=false;
	if (tokens[5]==1) {
		useRDN=true;
	}
	
	std::vector<float> frequencies;
	//frequencies.push_back(5e9);
	//frequencies.push_back(10e9);
	//frequencies.push_back(15e9);
	//frequencies.push_back(20e9);
	//frequencies.push_back(25e9);
	//frequencies.push_back(30e9);
	//frequencies.push_back(2e9);
	frequencies.push_back(2.5e9);
	frequencies.push_back(5e9);
	frequencies.push_back(7.5e9);
	frequencies.push_back(10e9);
	frequencies.push_back(12.5e9);
	frequencies.push_back(15e9);
	frequencies.push_back(17.5e9);
	frequencies.push_back(20e9);
	frequencies.push_back(22.5e9);
	frequencies.push_back(25e9);
	frequencies.push_back(27.5e9);
	std::vector<float> radius;
	if (useRDN) {
		radius.push_back(1.5);
		radius.push_back(1.0);
		radius.push_back(0.5);
		radius.push_back(0.25);
		radius.push_back(0.1);
		radius.push_back(0.05);
	} else {
		radius.push_back(3.5);
		radius.push_back(2.5);
		radius.push_back(1.0);
		radius.push_back(0.5);
		radius.push_back(0.1);
		radius.push_back(0.05);
	}
	freq=frequencies[tokens[0]];
	gsufix=std::to_string(tokens[0]);
	std::string path;
	if (tokens[1]==0) {
		//path="parking/parking-terrain-vehicles.json";
		//path="parking/lora-vehicles-improved.json";
		path="parking/parking-upct-smoothed.json";
	} else {
		//path="parking/parking-terrain-empty.json";
		//path="parking/lora-empty-improved.json";
		path="parking/parking-upct-smoothed-vehicles.json";
	}
	useGain=false;
	if (tokens[2]==1) {
		useGain=true;
	}
        this->sphereRadius=radius[tokens[3]];
        
        bool useRandom=false;
	if (tokens[4]==0) {
		useRandom=false;
	} else {
		useRandom=true;
	}

        std::cout<<"useRandom="<<useRandom<<"useGain="<<useGain<<"useRDN="<<useRDN<<"radius="<<sphereRadius<<std::endl;   
	//loadScenarioParking(useRandom,useRDN, path);
	loadScenarioParkingMultifrequency(useRandom,useRDN, path, frequencies, radius, outputFile);
}
void Parking::test() {
	this->replications=1;
	
	freq=5.9e9;
	std::string path("test-lora.json");;
	bool useRDN=false;
	useGain=false;
	loadScenarioParking(false,useRDN, path);
}

