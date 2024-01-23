/***************************************************************/
//
//Copyright (c) 2021 Esteban Egea-Lopez http://girtel.upct.es/~simulations
//
/**************************************************************/
#include "lora.h"
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
Lora::Lora(OpalSceneManager*   sceneManager, float sphereRadius, bool useDepolarization) : BasicTests(sceneManager, sphereRadius, useDepolarization) {
	this->freq=868e6;
	//Location of transmitter in our measurements
	//this->baseStation=make_float3(1543.56f,37.6f,620.39f);
	//this->baseStation=make_float3(1527.6,39.47,634.26);
	//this->baseStation=make_float3(1543.56f,37.2f,620.45f);

	//Both below are ELDI fixed Lora gateway
	//this->baseStation=make_float3(1463.48f,38.99f,617.73f);i
	this->baseStation=make_float3(1462.1f,38.31f,618.048f);
	this->gsufix=std::to_string(0);
	this->replications=1;
	this->useGain=true;
	this->rdnSim = nullptr;
	this->lpflatSim = nullptr;
	this->basicSim = nullptr;
	this->diffSim = nullptr;
	this->gIdRx=-1;
	this->gIdTx=-1;
	this->txPower=1;
	
}
void  Lora::runPoint(bool addRandom, bool useRDN, std::string path) {
	//A few receivers
	//float3 baseStation=make_float3(1543.56f,37.6f,620.39f);
	Timer timer;
	Timer timerRT;
	//float freq = 1800e6f;
	std::cout<<"Load LoRa with scenario "<<path<<"; freq="<<freq<<";tx="<<baseStation<<std::endl;
	timer.start();	
	//Init context before doing anything else
	//sceneManager->enableGenerateRaysOnLaunch();
	sceneManager->setMinEpsilon(1e-3);
	sceneManager->setUseAntennaGain(true);
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
			sim->setEnableTraceLog(true);
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
	simd->setEnableTraceLog(true);
	//simd->setPrintHits(true);
	simd->setEnableSimulation(true);

	sceneManager->initContext(freq);
	//Exceptions
	//sceneManager->enableExceptions();	

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
	int gainIdTx=-1;
	int gainIdRx=-1;
	if (useGain) {	
		//std::string gp("parking/gainT");
		//std::string gainPathT=gp+gsufix+".txt";
		//std::string gainPathT("lora/g11467.txt");
		std::string gainPathT("lora/dipole.txt");
		AntennaGain gains=sceneManager->loadGainsFromFileIndBPower(gainPathT.c_str());
		gainIdTx=sceneManager->registerAntennaGain(gains);
		//Receiver gain
		//std::string gpR("parking/gainR");
		//std::string gainPathR=gpR+gsufix+".txt";
		//std::string gainPathR("lora/g17514.txt");
		std::string gainPathR("lora/dipole.txt");
		AntennaGain gainsRx=sceneManager->loadGainsFromFileIndBPower(gainPathR.c_str());
		gainIdRx=sceneManager->registerAntennaGain(gainsRx);
	}

	optix::float3 polarization = make_float3(0.0f, 1.0f, 0.0f); 


	//Use antenna in ELDI as tx
	//float3 postx=make_float3(1527.6,39.47,634.26);
	float3 postx=baseStation;
	//std::string rx_file("lora/rx-medidas.txt");
	//std::vector<float3> rx=loadReceiversFromFile(rx_file);
	//std::default_random_engine gen;
	//std::uniform_real_distribution<float> df(-0.5f,0.5f);
	//std::uniform_real_distribution<float> dfy(-0.2f,0.2f);
	int i=0;
	//In LoRA the GW is a receiver
	//Location in ELDI of current GW
	//float3 posrx=make_float3(1527.6,39.47,634.26);
	//Nodes are transmitters
	//float3 postx = make_float3(1533.1f, 20.159f, 674.5f);
	//float3 postx = make_float3(1626.8f,28.123f,689.97f);

	////float3 posrx=make_float3(1543.5,37,620.45);
	float3 posrx=make_float3(1657.51001, 28.33124924, 626.5900269);
	//float3 posrx=make_float3(1619.5,22.85484,580.6);
	sceneManager->addReceiver(0, posrx,polarization, sphereRadius, sceneManager->printPower);
	if (useGain) {
		sceneManager->registerReceiverGain(0,gainIdRx);
	}
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
	//	optix::float3 postx = tx_f;
	//optix::float3 postx = baseStation;
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
void  Lora::loadScenarioLora(bool addRandom, bool useRDN, std::string path) {
	//A few receivers
	//float3 baseStation=make_float3(1543.56f,37.6f,620.39f);
	Timer timer;
	Timer timerRT;
	//float freq = 1800e6f;
	std::cout<<"Load LoRa with scenario "<<path<<"; freq="<<freq<<";tx="<<baseStation<<std::endl;
	timer.start();	
	//Init context before doing anything else
	//sceneManager->enableGenerateRaysOnLaunch();
	sceneManager->setMinEpsilon(1e-3);
	sceneManager->setUseAntennaGain(true);
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
			//sim->setEnableTraceLog(true);
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
	//simd->setEnableTraceLog(true);
	//simd->setPrintHits(true);
	simd->setEnableSimulation(true);

	sceneManager->initContext(freq);
	//Exceptions
	//sceneManager->enableExceptions();	

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
	int gainIdTx=-1;
	int gainIdRx=-1;
	if (useGain) {	
		//std::string gp("parking/gainT");
		//std::string gainPathT=gp+gsufix+".txt";
		//std::string gainPathT("lora/g11467.txt");
		std::string gainPathT("lora/dipole.txt");
		AntennaGain gains=sceneManager->loadGainsFromFileIndBPower(gainPathT.c_str());
		gainIdTx=sceneManager->registerAntennaGain(gains);
		//Receiver gain
		//std::string gpR("parking/gainR");
		//std::string gainPathR=gpR+gsufix+".txt";
		//std::string gainPathR("lora/g17514.txt");
		std::string gainPathR("lora/dipole.txt");
		AntennaGain gainsRx=sceneManager->loadGainsFromFileIndBPower(gainPathR.c_str());
		gainIdRx=sceneManager->registerAntennaGain(gainsRx);
	}

	optix::float3 polarization = make_float3(0.0f, 1.0f, 0.0f); 


	//Use antenna in ELDI as tx
	//float3 postx=make_float3(1527.6,39.47,634.26);
	float3 postx=baseStation;
	//std::string rx_file("lora/rx-tester-all.txt");
	//std::string rx_file("lora/rx-an-all.txt");
	std::string rx_file("lora/dec22/rx-dec22.txt");
	std::vector<float3> rx=loadReceiversFromFile(rx_file);
	std::default_random_engine gen;
	std::uniform_real_distribution<float> df(-0.5f,0.5f);
	std::uniform_real_distribution<float> dfy(-0.2f,0.2f);
	int i=0;
	//In LoRA the GW is a receiver
	//Location in ELDI of current GW
	//float3 posrx=make_float3(1527.6,39.47,634.26);
	//Nodes are transmitters
	//float3 postx = make_float3(1533.1f, 20.159f, 674.5f);
	//float3 postx = make_float3(1626.8f,28.123f,689.97f);

	////float3 posrx=make_float3(1543.5,37,620.45);
	//float3 posrx=make_float3(1620.13,23.03,611.18);
	////float3 posrx=make_float3(1619.5,22.85484,580.6);
	//sceneManager->addReceiver(0, posrx,polarization, sphereRadius, sceneManager->printPower);
	//if (useGain) {
	//	sceneManager->registerReceiverGain(0,gainIdRx);
	//}
	//sceneManager->addReceiver(0, posrx,polarization, sphereRadius, sceneManager->printPower);
	//sceneManager->registerReceiverGain(0,gainId);
	//***Single ray transmit****
	//float3 mRay=normalize(make_float3(-0.004128737841, -0.9902680516, 0.1391119212));
	//sceneManager->createRaySphere2D(1,1,&mRay);

	int totalRx = rx.size();
	int maxBatch=590;
	std::vector<float3> allReceivers; 
	if (addRandom) {
		totalRx = rx.size()*replications;
	}
	std::cout<<"totalRx="<<totalRx<<";maxBatch="<<maxBatch<<std::endl;
	if (totalRx <maxBatch) {

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

	} else {
		for (auto p : rx) {
			if (addRandom) {
				for (int j=0; j<replications; j++) {
					float3 posrx = make_float3(p.x,  p.y, p.z);
					posrx.x=posrx.x+df(gen);
					posrx.y=posrx.y+dfy(gen);
					posrx.z=posrx.z+df(gen);
					allReceivers.push_back(posrx);

				}
			} else {
				float3 posrx = make_float3(p.x,  p.y, p.z);
				allReceivers.push_back(posrx);
			}
		}
	}


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
		//sceneManager->createRaySphere2D(0.0f,0.01,180.0f,0.0f,0.01,360.0f);
	}
	//	//float3 postx = make_float3(1501.6, 25.1f, 609.9f);
	//optix::float3 postx = make_float3(1547.14f, 40.69f, 620.8f);
	//	optix::float3 postx = tx_f;
	//optix::float3 postx = baseStation;
	timerRT.start();
	if (useGain) {
		sceneManager->registerTransmitterGain(i+1,gainIdTx);
	}
	if (totalRx<maxBatch) {
		sceneManager->transmit(i+1, 0.0251188, postx, polarization);
	} else {
		int i=0;
		int j=0;
		for (auto p : allReceivers) {
			sceneManager->addReceiver(j, p,polarization, sphereRadius, sceneManager->printPower);
			if (useGain) {
				sceneManager->registerReceiverGain(j,gainIdRx);
			}
			++i;
			++j;
			if (i==maxBatch) {
				//sceneManager->transmit(j+1, 0.0251, postx, polarization);
				sceneManager->transmit(j+1, 1.0f, postx, polarization);
				i=0;
				std::cout<<"Clearing receivers"<<std::endl;
				sceneManager->clearReceivers();
			}
		}
		sceneManager->transmit(j+1, 0.0251188, postx, polarization);
	}
	timer.stop();
	timerRT.stop();
	std::cout<<"Time\t"<<timer.getTime()<<"\t"<<timerRT.getTime()<<std::endl;
	delete loader;
}
std::vector<float3> Lora::loadReceiversFromFile(std::string file) {
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
void Lora::runTests(std::string test) {
	std::vector<float> frequencies;
	std::vector<float> radius;
	radius.push_back(0.05);
	radius.push_back(0.5);
	radius.push_back(0.1);
	radius.push_back(1.0);
	radius.push_back(2.5);
	radius.push_back(3.5);

	frequencies.push_back(868e6);
	//frequencies.push_back(2e9);
	//frequencies.push_back(5e9);
	//frequencies.push_back(10e9);
	//frequencies.push_back(15e9);
	//frequencies.push_back(20e9);
	//frequencies.push_back(25e9);
	//frequencies.push_back(30e9);

	std::vector<int> tokens=parseTestString(test);

	this->replications=5;

	freq=frequencies[0];
	std::string path;
	//path="lora/patio.json";
	//path="lora/lora-empty-edges.json";
	path="lora/lora-empty-improved.json";
	if (tokens[0]==1) {
		//path="lora/lora-vehicles-edges.json";
		path="lora/lora-vehicles-improved.json";
	}
	bool useRandom=false;
	if (tokens[1]==0) {
		useRandom=false;
	} else {
		useRandom=true;
	}
	bool useRDN=false;
	if (tokens[2]==1) {
		useRDN=true;
	}
	useGain=true;
	if (tokens[3]==0) {
		useGain=false;
	}
	this->sphereRadius=radius[tokens[4]];
	std::cout<<"useRandom="<<useRandom<<"useGain="<<useGain<<"useRDN="<<useRDN<<"radius="<<sphereRadius<<std::endl;   
	loadScenarioLora(useRandom,useRDN, path);
	//runPoint(useRandom,useRDN, path);
}
void Lora::runAll(bool empty,std::string test, std::string folder) {
	std::vector<float> frequencies;
	std::vector<float> radius;
	radius.push_back(0.05);
	radius.push_back(0.5);
	radius.push_back(0.1);
	radius.push_back(1.0);
	radius.push_back(2.5);
	radius.push_back(3.5);
	//radius.push_back(4.5);
	//radius.push_back(5);

	frequencies.push_back(868e6);
	//frequencies.push_back(2e9);
	//frequencies.push_back(5e9);
	//frequencies.push_back(10e9);
	//frequencies.push_back(15e9);
	//frequencies.push_back(20e9);
	//frequencies.push_back(25e9);
	//frequencies.push_back(30e9);

	//this->txPower=1.0f;
	this->txPower=0.0251188f;
	this->replications=10;
	std::vector<int> tokens=parseTestString(test);

	freq=frequencies[0];
	//std::string path="lora/lora-empty-improved.json";
	std::string path="lora/lora-empty-meas-dec22.json";
	if (tokens[0]==1) {
		//path="lora/lora-vehicles-edges.json";
		path="lora/lora-parking-meas-dec22.json";
	}
	bool useRDN=false;
	if (tokens[1]==1) {
		useRDN=true;
	}
	buildAllScenario(useRDN, path );

	//empty or parking
		//random or not
		for (int r=0; r<1; r++) {
				for (int rad=4; rad<5; rad++) {
				//for (int rad=0; rad<radius.size(); rad++) {
				//for (int rad=7; rad<8; rad++) {
					std::string path;
					//path="lora/patio.json";
					//path="lora/lora-empty-edges.json";
					bool useRandom=false;
					if (r==0) {
						useRandom=false;
					} else {
						useRandom=true;
					}
					useGain=true;
					this->sphereRadius=radius[rad];
					std::cout<<"useRandom="<<useRandom<<"useGain="<<useGain<<"useRDN="<<useRDN<<"radius="<<sphereRadius<<std::endl;   
					runScenarioVariant(useRandom,useRDN, folder);
					std::cout<<"Clearing receivers"<<std::endl;
					sceneManager->clearReceivers();
				}
		}
	//runPoint(useRandom,useRDN, path);
}
void Lora::runAllMoving(bool empty,std::string test, std::string folder) {
	std::vector<float> frequencies;
	std::vector<float> radius;
	radius.push_back(0.05);
	radius.push_back(0.5);
	radius.push_back(0.1);
	radius.push_back(1.0);
	radius.push_back(2.5);
	radius.push_back(3.5);
	//radius.push_back(4.5);
	//radius.push_back(5);

	frequencies.push_back(868e6);
	//frequencies.push_back(2e9);
	//frequencies.push_back(5e9);
	//frequencies.push_back(10e9);
	//frequencies.push_back(15e9);
	//frequencies.push_back(20e9);
	//frequencies.push_back(25e9);
	//frequencies.push_back(30e9);

	//this->txPower=1.0f;
	this->txPower=0.0251188f;
	this->replications=10;
	std::vector<int> tokens=parseTestString(test);

	freq=frequencies[0];
	std::string path="lora/lora-empty-improved.json";
	//if (tokens[0]==1) {
	//	//path="lora/lora-vehicles-edges.json";
	//	path="lora/lora-vehicles-improved.json";
	//}
	bool useRDN=false;
	if (tokens[1]==1) {
		useRDN=true;
	}
	buildAllScenario(useRDN, path );

		//random or not
		for (int r=0; r<2; r++) {
				//for (int rad=0; rad<1; rad++) {
				for (int rad=0; rad<radius.size(); rad++) {
				//for (int rad=7; rad<8; rad++) {
					std::string path;
					//path="lora/patio.json";
					//path="lora/lora-empty-edges.json";
					bool useRandom=false;
					if (r==0) {
						useRandom=false;
					} else {
						useRandom=true;
					}
					useGain=true;
					this->sphereRadius=radius[rad];
					std::cout<<"useRandom="<<useRandom<<"useGain="<<useGain<<"useRDN="<<useRDN<<"radius="<<sphereRadius<<std::endl;   
					runScenarioVariantMovingTransmitter(useRandom,useRDN, folder);
					std::cout<<"Clearing receivers"<<std::endl;
					sceneManager->clearReceivers();
				}
		}
	//runPoint(useRandom,useRDN, path);
}
void  Lora::buildAllScenario(bool useRDN, std::string path ) {
	//A few receivers
	//float3 baseStation=make_float3(1543.56f,37.6f,620.39f);
	Timer timer;
	//float freq = 1800e6f;
	std::cout<<"Load LoRa with scenario "<<path<<"; freq="<<freq<<";tx="<<baseStation<<std::endl;
	timer.start();	
	//Init context before doing anything else
	//sceneManager->enableGenerateRaysOnLaunch();
	sceneManager->setMinEpsilon(1e-3);
	sceneManager->setUseAntennaGain(true);
	sceneManager->enableGenerateRaysOnLaunch();	
	ComputeMode mode=ComputeMode::FIELD;
	//ComputeMode mode=ComputeMode::VOLTAGE;
	if (useRDN) {
		RayDensityNormalizationSimulation* sim = new RayDensityNormalizationSimulation(sceneManager);
		rdnSim=sim;
		sim->setComputeMode(mode);
		sceneManager->setSimulation(sim);
	} else {
		if (useDepolarization) {
			LPFlatMeshReflectionSimulation* sim = new LPFlatMeshReflectionSimulation(sceneManager);
			sceneManager->setSimulation(sim);
			sim->setComputeMode(mode);
			lpflatSim=sim;
			//sim->setEnableTraceLog(true);
			//sim->setEnableSimulation(false);
			//sim->setPrintHits(true);
		} else {
			BasicFlatMeshReflectionSimulation* sim = new BasicFlatMeshReflectionSimulation(sceneManager);
			sceneManager->setSimulation(sim);
			sim->setComputeMode(mode);
			basicSim=sim;
			//sim->setPrintHits(true);
		}
	}
	//Add diffraction
	SingleDiffraction* simd= new SingleDiffraction(sceneManager);
	sceneManager->setSimulation(simd);
	simd->setComputeMode(mode);
	//simd->setEnableTraceLog(true);
	//simd->setPrintHits(true);
	simd->setEnableSimulation(true);
	diffSim=simd;

	sceneManager->initContext(freq);
	//Exceptions
	//sceneManager->enableExceptions();	

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
	int gainIdTx=-1;
	int gainIdRx=-1;
	if (useGain) {	
		//std::string gp("parking/gainT");
		//std::string gainPathT=gp+gsufix+".txt";
		//std::string gainPathT("lora/g11467.txt");
		std::string gainPathT("lora/dipole.txt");
		AntennaGain gains=sceneManager->loadGainsFromFileIndBPower(gainPathT.c_str());
		gainIdTx=sceneManager->registerAntennaGain(gains);
		gIdTx=gainIdTx;
		//Receiver gain
		//std::string gpR("parking/gainR");
		//std::string gainPathR=gpR+gsufix+".txt";
		//std::string gainPathR("lora/g17514.txt");
		std::string gainPathR("lora/dipole.txt");
		AntennaGain gainsRx=sceneManager->loadGainsFromFileIndBPower(gainPathR.c_str());
		gainIdRx=sceneManager->registerAntennaGain(gainsRx);
		gIdRx=gainIdRx;
	}
	sceneManager->finishSceneContext();
	timer.stop();
	std::cout<<"Time to build scenario\t"<<timer.getTime()<<std::endl;
	delete loader;

}
void  Lora::runScenarioVariantMovingTransmitter(bool addRandom, bool useRDN,  std::string folder) {
	
	Timer timerRT;
	Timer timer;
	//float freq = 1800e6f;
	timer.start();
	timerRT.start();	
	std::string type="rdn";
	std::string rand="norand";
	if (addRandom) {
		rand="random";
	}
	if (useRDN) {

	} else {
		type="dep";


	}
	int gainIdTx=gIdTx;
	int gainIdRx=gIdRx;
	optix::float3 polarization = make_float3(0.0f, 1.0f, 0.0f); 
	//Use antenna in ELDI as tx
	//float3 postx=make_float3(1527.6,39.47,634.26);
	float3 postx=baseStation;
	//std::string rx_file("lora/rx-tester-all.txt");
	//std::string rx_file("lora/rx-an-all.txt");
	std::string rx_file("lora/dec22/rx-dec22.txt");
	std::vector<float3> rx=loadReceiversFromFile(rx_file);
	std::default_random_engine gen;
	std::uniform_real_distribution<float> df(-0.5f,0.5f);
	std::uniform_real_distribution<float> dfy(-0.2f,0.2f);
	//In LoRA the GW is a receiver
	//Location in ELDI of current GW

	//Now ELDI is a receiver
	sceneManager->addReceiver(0, postx,polarization, sphereRadius, sceneManager->printPower);
	if (useGain) {
		sceneManager->registerReceiverGain(0,gainIdRx);
	}
	if (useRDN) {
		int rayD=10000;
		sceneManager->setRayRange(0.0,180.0,0.0,360.0,rayD,rayD);
		RayDensityNormalizationSimulation* s=dynamic_cast<RayDensityNormalizationSimulation*>(sceneManager->getSimulation(0));
		s->setInitialDensity(((float)sceneManager->getRaySphere().rayCount)/(4*M_PIf));
		s->setFiltering(2u);

	} else {
		sceneManager->createRaySphere2D(0.0f,0.1,180.0f,0.0f,0.1,360.0f);
		//sceneManager->createRaySphere2D(0.0f,0.01,180.0f,0.0f,0.01,360.0f);
	}

	int totalRx = rx.size();
	int maxBatch=590;
	std::vector<float3> allReceivers; 
	if (addRandom) {
		totalRx = rx.size()*replications;
	}
	std::cout<<"totalRx="<<totalRx<<";maxBatch="<<maxBatch<<std::endl;
	ResultReport report;
	int i=1;
	//if (totalRx <maxBatch) {

		for (auto p : rx) {
			if (addRandom) {
				for (int j=0; j<replications; j++) {
					float3 posrx = make_float3(p.x,  p.y, p.z);
					posrx.x=posrx.x+df(gen);
					posrx.y=posrx.y+dfy(gen);
					posrx.z=posrx.z+df(gen);
					if (useGain) {
						sceneManager->registerTransmitterGain(i,gainIdTx);
					}

					ResultReport* rp=sceneManager->transmit(i, txPower, posrx, polarization);
					report.merge(*rp);
					++i;
				}
			} else {
				float3 posrx = make_float3(p.x,  p.y, p.z);
				//float d=length(postx-posrx);
				float radius=this->sphereRadius;
				//float radius=3.5;
				//if (d<=250) {
				//	radius=2;
				//} else if (d<=350) {
				//	radius=2.5;
				//}
				//sceneManager->addReceiver(i, posrx,polarization, sphereRadius, sceneManager->printPower);
				if (useGain) {
					sceneManager->registerTransmitterGain(i,gainIdTx);
				}
				ResultReport* rp=sceneManager->transmit(i, txPower, posrx, polarization);
				report.merge(*rp);
				++i;
			}
		}

	//} 


	std::string savePath=folder+"/vol-"+type+"-"+rand+"-"+std::to_string(sphereRadius)+".csv";
	std::cout<<"Saving results to "<<savePath<<std::endl;
	report.toCSV(savePath);
	timer.stop();
	timerRT.stop();
	std::cout<<"Time\t"<<timer.getTime()<<"\t"<<timerRT.getTime()<<std::endl;
}


void  Lora::runScenarioVariant(bool addRandom, bool useRDN,  std::string folder) {
	
	Timer timerRT;
	Timer timer;
	//float freq = 1800e6f;
	timer.start();
	timerRT.start();	
	std::string type="rdn";
	std::string rand="norand";
	if (addRandom) {
		rand="random";
	}
	if (useRDN) {

	} else {
		type="dep";


	}
	int gainIdTx=gIdTx;
	int gainIdRx=gIdRx;
	optix::float3 polarization = make_float3(0.0f, 1.0f, 0.0f); 
	//Use antenna in ELDI as tx
	//float3 postx=make_float3(1527.6,39.47,634.26);
	float3 postx=baseStation;
	//std::string rx_file("lora/rx-tester-all.txt");
	//std::string rx_file("lora/rx-an-all.txt");
	std::string rx_file("lora/rx-dec22.txt");
	//std::string rx_file("lora/rx-dec22-one-receiver.txt");
	std::vector<float3> rx=loadReceiversFromFile(rx_file);
	std::default_random_engine gen;
	std::uniform_real_distribution<float> df(-0.5f,0.5f);
	std::uniform_real_distribution<float> dfy(-0.2f,0.2f);
	int i=0;
	//In LoRA the GW is a receiver
	//Location in ELDI of current GW
	//float3 posrx=make_float3(1527.6,39.47,634.26);
	//Nodes are transmitters
	//float3 postx = make_float3(1533.1f, 20.159f, 674.5f);
	//float3 postx = make_float3(1626.8f,28.123f,689.97f);

	////float3 posrx=make_float3(1543.5,37,620.45);
	//float3 posrx=make_float3(1620.13,23.03,611.18);
	////float3 posrx=make_float3(1619.5,22.85484,580.6);
	//sceneManager->addReceiver(0, posrx,polarization, sphereRadius, sceneManager->printPower);
	//if (useGain) {
	//	sceneManager->registerReceiverGain(0,gainIdRx);
	//}
	//sceneManager->addReceiver(0, posrx,polarization, sphereRadius, sceneManager->printPower);
	//sceneManager->registerReceiverGain(0,gainId);
	//***Single ray transmit****
	//float3 mRay=normalize(make_float3(-0.004128737841, -0.9902680516, 0.1391119212));
	//sceneManager->createRaySphere2D(1,1,&mRay);

	int totalRx = rx.size();
	int maxBatch=590;
	std::vector<float3> allReceivers; 
	if (addRandom) {
		totalRx = rx.size()*replications;
	}
	std::cout<<"totalRx="<<totalRx<<";maxBatch="<<maxBatch<<std::endl;
	if (totalRx <maxBatch) {

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
				//float d=length(postx-posrx);
				float radius=this->sphereRadius;
				//float radius=3.5;
				//if (d<=250) {
				//	radius=2;
				//} else if (d<=350) {
				//	radius=2.5;
				//}
				sceneManager->addReceiver(i, posrx,polarization, radius, sceneManager->printPower);
				//sceneManager->addReceiver(i, posrx,polarization, sphereRadius, sceneManager->printPower);
				if (useGain) {
					sceneManager->registerReceiverGain(i,gainIdRx);
				}
				++i;
			}
		}

	} else {
		for (auto p : rx) {
			if (addRandom) {
				for (int j=0; j<replications; j++) {
					float3 posrx = make_float3(p.x,  p.y, p.z);
					posrx.x=posrx.x+df(gen);
					posrx.y=posrx.y+dfy(gen);
					posrx.z=posrx.z+df(gen);
					allReceivers.push_back(posrx);

				}
			} else {
				float3 posrx = make_float3(p.x,  p.y, p.z);
				allReceivers.push_back(posrx);
			}
		}
	}


	if (useRDN) {
		int rayD=10000;
		sceneManager->setRayRange(0.0,180.0,0.0,360.0,rayD,rayD);
		RayDensityNormalizationSimulation* s=dynamic_cast<RayDensityNormalizationSimulation*>(sceneManager->getSimulation(0));
		s->setInitialDensity(((float)sceneManager->getRaySphere().rayCount)/(4*M_PIf));
		s->setFiltering(2u);

	} else {
		sceneManager->createRaySphere2D(0.0f,0.1,180.0f,0.0f,0.1,360.0f);
		//sceneManager->createRaySphere2D(0.0f,0.01,180.0f,0.0f,0.01,360.0f);
	}
	//	//float3 postx = make_float3(1501.6, 25.1f, 609.9f);
	//optix::float3 postx = make_float3(1547.14f, 40.69f, 620.8f);
	//	optix::float3 postx = tx_f;
	//optix::float3 postx = baseStation;
	if (useGain) {
		sceneManager->registerTransmitterGain(i+1,gainIdTx);
	}
	ResultReport report;
	if (totalRx<maxBatch) {
			ResultReport* rp=sceneManager->transmit(i+1, txPower, postx, polarization);
			report.merge(*rp);
	} else {
		int i=0;
		int j=0;
		for (auto p : allReceivers) {
			sceneManager->addReceiver(j, p,polarization, sphereRadius, sceneManager->printPower);
			if (useGain) {
				sceneManager->registerReceiverGain(j,gainIdRx);
			}
			++i;
			++j;
			if (i==maxBatch) {
				//sceneManager->transmit(j+1, 0.0251, postx, polarization);
				ResultReport* rp=sceneManager->transmit(j+1, txPower, postx, polarization);
				i=0;
				report.merge(*rp);
				std::cout<<"Clearing receivers"<<std::endl;
				sceneManager->clearReceivers();
			}
		}
		ResultReport* rp=sceneManager->transmit(j+1, txPower, postx, polarization);
		report.merge(*rp);
	}
	std::string savePath=folder+"/field-"+type+"-"+rand+"-"+std::to_string(sphereRadius)+".csv";
	std::cout<<"Saving results to "<<savePath<<std::endl;
	report.toCSV(savePath);
	timer.stop();
	timerRT.stop();
	std::cout<<"Time\t"<<timer.getTime()<<"\t"<<timerRT.getTime()<<std::endl;
}
void Lora::test() {
	this->replications=1;

	freq=868e6;
	std::string path("test-lora.json");;
	bool useRDN=false;
	useGain=false;
	loadScenarioLora(false,useRDN, path);
}

