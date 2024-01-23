/***************************************************************/
//
//Copyright (c) 2021 Esteban Egea-Lopez http://girtel.upct.es/~simulations
//
/**************************************************************/
#include "loraDoppler.h"
#include "../timer.h"
#include "../Opal.h"
#include "../configuration.h"
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
LoraDoppler::LoraDoppler(OpalSceneManager*   sceneManager, ConfigurationParser* config) : BasicTests(sceneManager, sphereRadius, useDepolarization) {
	// "Save" config file
	this->config = config;

	// Read some config parameters
	this->freq = 868e6;
	this->replications = (config != nullptr && !config->getKey("replications").is_null()) ? config->getKey("replications").get<uint>() : 1;
	this->useGain = (config != nullptr && !config->getKey("useGain").is_null()) ? config->getKey("useGain").get<bool>() : false;
	this->useDepolarization = (config != nullptr && !config->getKey("useDepolarization").is_null()) ? config->getKey("useDepolarization").get<bool>() : false;
	this->useRandom = (config != nullptr && !config->getKey("randomizePositions").is_null() )? config->getKey("randomizePositions").get<bool>() : false;
	this->txPower = (config != nullptr && !config->getKey("txPower").is_null()) ? config->getKey("txPower").get<float>() : 1;
	this->sphereRadius = 1;
	this->minEpsilon =(config != nullptr && !config->getKey("minEpsilon").is_null()) ? config->getKey("minEpsilon").get<float>() : 1e-3;


	float bSX = (config != nullptr && !config->getKey("baseStation").is_null()) ? config->getKey(config->getKey("baseStation"), "x").get<float>() : 1462.1f;
	float bSY = (config != nullptr && !config->getKey("baseStation").is_null()) ? config->getKey(config->getKey("baseStation"), "y").get<float>() : 38.31f;
	float bSZ = (config != nullptr && !config->getKey("baseStation").is_null()) ? config->getKey(config->getKey("baseStation"), "z").get<float>() : 618.048f;
	this->baseStation=make_float3(bSX,bSY,bSZ);

	this->gsufix=std::to_string(0);
	this->rdnSim = nullptr;
	this->lpflatSim = nullptr;
	this->basicSim = nullptr;
	this->diffSim = nullptr;
	this->gIdRx=-1;
	this->gIdTx=-1;

}
void  LoraDoppler::runPoint(bool addRandom, bool useRDN, std::string path) {

	Timer timer;
	Timer timerRT;
	std::cout<<"Load LoRa with scenario "<<path<<"; freq="<<freq<<";tx="<<baseStation<<std::endl;
	timer.start();
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
		//Transmitter gain
		std::string antenaGainTx = !config->getKey("gainPathTx").is_null() ? config->getKey("gainPathTx").get<std::string>() : "/lora/dipole.txt";
		AntennaGain gains = sceneManager->loadGainsFromFileIndBPower(antenaGainTx.c_str());
		gainIdTx=sceneManager->registerAntennaGain(gains);

        //Receiver gain
		std::string antenaGainRx = !config->getKey("gainPathRx").is_null() ? config->getKey("gainPathRx").get<std::string>() : "/lora/dipole.txt";
		AntennaGain gainsRx=sceneManager->loadGainsFromFileIndBPower(antenaGainRx.c_str());
		gainIdRx=sceneManager->registerAntennaGain(gainsRx);
	}

	optix::float3 polarization = make_float3(0.0f, 1.0f, 0.0f);


	//Use antenna in ELDI as tx
	float3 postx=baseStation;

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

void  LoraDoppler::loadScenarioLora(bool addRandom, bool useRDN, std::string path) {
	Timer timer;
	Timer timerRT;
	std::cout<<"Load LoRa with scenario "<<path<<"; freq="<<freq<<";tx="<<baseStation<<std::endl;
	timer.start();
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
		std::string antenaGainTx = !config->getKey("gainPathTx").is_null() ? config->getKey("gainPathTx").get<std::string>() : "/lora/dipole.txt";
		AntennaGain gains = sceneManager->loadGainsFromFileIndBPower(antenaGainTx.c_str());
		gainIdTx=sceneManager->registerAntennaGain(gains);
		//Receiver gain
		std::string antenaGainRx = !config->getKey("gainPathRx").is_null() ? config->getKey("gainPathRx").get<std::string>() : "/lora/dipole.txt";
		AntennaGain gainsRx=sceneManager->loadGainsFromFileIndBPower(antenaGainRx.c_str());
		gainIdRx=sceneManager->registerAntennaGain(gainsRx);
	}

	optix::float3 polarization = make_float3(0.0f, 1.0f, 0.0f);


	//Use antenna in ELDI as tx
	//float3 postx=make_float3(1527.6,39.47,634.26);
	float3 postx = baseStation;
	//std::string rx_file("lora/rx-tester-all.txt");
	//std::string rx_file("lora/rx-an-all.txt");

	std::vector<float3> rx{make_float3(0.0f, 0.0f, 0.0f)};
	if (!config->getKey("useNodesFile").is_null() && config->getKey("useNodesFile").get<bool>() == false) {
		if (!config->getKey("nodes").is_null()) {
			rx = loadReceiversFromJSON();
		}
	}
	else {
		std::string path = !config->getKey("nodesPath").is_null() ? config->getKey("nodesPath").get<std::string>() : "lora/dec22/rx-dec22.txt";
		rx = loadReceiversFromFile(path);
	}



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
	if (totalRx < maxBatch) {

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

std::vector<float3> LoraDoppler::loadReceiversFromFile(std::string file) {
	// Load Receivers from File
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

std::vector<float3> LoraDoppler::loadReceiversFromJSON() {
	// Load Receivers from JSON (when we use configuration file)
	std::vector<nlohmann::json> nodes = config->getKey("nodes").get<std::vector<nlohmann::json>>();
	std::vector<float3> rx;
	std::cout<<"Loading receivers from JSON"<<std::endl;
	for (nlohmann::json n : nodes) {
		float x = config->getKey(n, "x").get<float>();
		float y = config->getKey(n, "y").get<float>();
		float z = config->getKey(n, "z").get<float>();
		optix::float3 nodePosition = make_float3(x, y, z);
		rx.push_back(nodePosition);
	}
	return rx;
}

void LoraDoppler::runTests(std::string test) {
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

void LoraDoppler::runAll(bool empty, std::string test, std::string folder) {

	// Radius vector
	std::vector<float> radii = (config != nullptr && !config->getKey("sphereRadius").is_null()) ? config->getKey("sphereRadius").get<std::vector<float>>() : std::vector<float>{1};

	// Freq vector
	std::vector<float> frequencies = (config != nullptr && !config->getKey("freq").is_null()) ? config->getKey("freq").get<std::vector<float>>() : std::vector<float>{868e8};

	// Path for scenario
	std::string path = (config != nullptr && !config->getKey("scenarioPath").is_null()) ? config->getKey("scenarioPath").get<std::string>() : "lora/lora-empty-meas-dec22.json";

	bool useRDN = (config != nullptr && !config->getKey("useRDN").is_null()) ? config->getKey("useRDN").get<bool>() : false;

	buildAllScenario(useRDN, path);
	// We combine all frequencies and radii to perform all the simulations
	for (float frequency: frequencies) {
		this->freq = frequency;
		for (float radius : radii) {
			this->sphereRadius = radius;
			std::cout<<"Running simulation: "<< "\n- Freq: " << frequency << "\n- Radius: " << radius << "\n- useRandom: "<<this->useRandom<<"\n- useRDN: "<<useRDN <<std::endl;
			runScenarioVariant(this->useRandom, useRDN, folder);
			std::cout<<"Clearing receivers"<<std::endl;
			sceneManager->clearReceivers();
		}
	}
}

void LoraDoppler::runAllMoving(bool empty,std::string test, std::string folder) {
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

void  LoraDoppler::buildAllScenario(bool useRDN, std::string path ) {
		Timer timer;
		std::cout<<"Load LoRa with scenario "<<path<<"; freq="<<freq<<";tx="<<baseStation<<std::endl;
		timer.start();

		sceneManager->setMinEpsilon(this->minEpsilon);
		sceneManager->setUseAntennaGain(this->useGain);
		sceneManager->enableGenerateRaysOnLaunch();
		sceneManager->enableMultiChannel();
		ComputeMode mode=ComputeMode::FIELD;
		// ComputeMode mode=ComputeMode::VOLTAGE;

		if (useRDN) {
			RayDensityNormalizationSimulation* sim = new RayDensityNormalizationSimulation(sceneManager);
			rdnSim=sim;
			sim->setComputeMode(mode);
			sceneManager->setSimulation(sim);
		} else if (useDepolarization){
			LPFlatMeshReflectionSimulation* sim = new LPFlatMeshReflectionSimulation(sceneManager);
			sceneManager->setSimulation(sim);
			sim->setComputeMode(mode);
			lpflatSim=sim;
		} else {
			BasicFlatMeshReflectionSimulation* sim = new BasicFlatMeshReflectionSimulation(sceneManager);
			sceneManager->setSimulation(sim);
			sim->setComputeMode(mode);
			basicSim=sim;
		}

		//Add diffraction
		SingleDiffraction* simd= new SingleDiffraction(sceneManager);
		sceneManager->setSimulation(simd);
		simd->setComputeMode(mode);
		simd->setEnableSimulation(true);
		diffSim=simd;

		sceneManager->initContext(freq);
		sceneManager->enableExceptions();
		sceneManager->setPrintEnabled(1023*1023*1203, make_uint3(1179, 2123,0));

		//Load files here
		ScenarioLoader* loader=new ScenarioLoader(sceneManager);
		loader->loadJSONScenario(path);

		int gainIdTx=-1;
		int gainIdRx=-1;
		if (useGain) {
			// Transmitter gain
			std::string antenaGainTx = !config->getKey("gainPathTx").is_null() ? config->getKey("gainPathTx").get<std::string>() : "/lora/dipole.txt";
			AntennaGain gains = sceneManager->loadGainsFromFileIndBPower(antenaGainTx.c_str());
			gainIdTx=sceneManager->registerAntennaGain(gains);
			gIdTx=gainIdTx;

			// Receiver gain
			std::string antenaGainRx = !config->getKey("gainPathRx").is_null() ? config->getKey("gainPathRx").get<std::string>() : "/lora/dipole.txt";
			AntennaGain gainsRx=sceneManager->loadGainsFromFileIndBPower(antenaGainRx.c_str());
			gainIdRx=sceneManager->registerAntennaGain(gainsRx);
			gIdRx=gainIdRx;
		}

		sceneManager->finishSceneContext();
		// RDN use
		if (useRDN) {
			int rayD=30000;
			sceneManager->setRayRange(0.0,180.0,0.0,360.0,rayD,rayD);
			RayDensityNormalizationSimulation* s=dynamic_cast<RayDensityNormalizationSimulation*>(sceneManager->getSimulation(0));
			s->setInitialDensity(((float)sceneManager->getRaySphere().rayCount)/(4*M_PIf));
			s->setFiltering(2u);
		} else {
			sceneManager->createRaySphere2D(0.0f,0.1,180.0f,0.0f,0.1,360.0f);
		}

		timer.stop();
		std::cout<<"Time to build scenario\t"<<timer.getTime()<<std::endl;
		delete loader;
	}

void  LoraDoppler::runScenarioVariantMovingTransmitter(bool addRandom, bool useRDN,  std::string folder) {

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

		std::vector<float3> rx{make_float3(0.0f, 0.0f, 0.0f)};
		if (config != nullptr && !config->getKey("useNodesFile").is_null() && config->getKey("useNodesFile").get<bool>() == false) {
			if (!config->getKey("nodes").is_null()) {
				rx = loadReceiversFromJSON();
			}
		}
		else {
			std::string path = !config->getKey("nodesPath").is_null() ? config->getKey("nodesPath").get<std::string>() : "lora/dec22/rx-dec22.txt";
			rx = loadReceiversFromFile(path);
		}

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

void  LoraDoppler::runScenarioVariant(bool addRandom, bool useRDN,  std::string folder) {
        std::string baseName = (config != nullptr && !config->getKey("baseName").is_null()) ? config->getKey("baseName").get<std::string>() : "simulation";
		Timer timerRT;
		Timer timer;
		timer.start();
		timerRT.start();
		std::string type="rdn";
		std::string rand="norand";

		if (addRandom) {
			rand="random";
		}
		if (useRDN) {
			type="rdn";
		} else {
			type="dep";
		}

		int gainIdTx=gIdTx;
		int gainIdRx=gIdRx;

		float3 postx = baseStation;

		std::vector<float3> rx{make_float3(0.0f, 0.0f, 0.0f)};
		if (!config->getKey("useNodesFile").is_null() && config->getKey("useNodesFile").get<bool>() == false) {
			if (!config->getKey("nodes").is_null()) {
				rx = loadReceiversFromJSON();
			}
		}
		else {
			std::string path = !config->getKey("nodesPath").is_null() ? config->getKey("nodesPath").get<std::string>() : "lora/dec22/rx-dec22.txt";
			rx = loadReceiversFromFile(path);
		}

		unsigned int maxReceivers=rx.size();
		int totalRx = rx.size();
		if (!config->getKey("batchSize").is_null() ) {
			maxReceivers=config->getKey("batchSize").get<unsigned int>();

		}
		// Random movement
		std::default_random_engine gen;
		std::uniform_real_distribution<float> df(-0.5f,0.5f);
		std::uniform_real_distribution<float> dfy(-0.2f,0.2f);


		if (addRandom) {
			totalRx = rx.size()*replications;
		}

		int numBatches=ceil(totalRx/(float) maxReceivers);
		std::cout <<"Number of required batches:"<<numBatches<<std::endl;
		ResultReport report;
		int i=0;

        optix::float3 polarization = make_float3(0.0f, 1.0f, 0.0f);

		for (int b=0; b < numBatches; ++b) {
			std::cout<<"Running batch "<<b<<"/"<<numBatches<<". Total receivers="<<totalRx<<"; Max receivers per batch="<<maxReceivers<<std::endl;
			std::vector<float3> aux;
			for (int r=0; r <maxReceivers; r++) {
				int index=(b*maxReceivers)+r;
				if (index<rx.size()) {
					aux.push_back(rx[index]);
				}
			}
			for (auto p : aux) {
				if (addRandom) {
					for (int j=0; j<replications; j++) {
						float3 posrx=make_float3(0,0,0);
						posrx.x=p.x+df(gen);
						posrx.y=p.y+dfy(gen);
						posrx.z=p.z+df(gen);
						sceneManager->addReceiver(i, posrx,polarization, sphereRadius, sceneManager->printPower);
						if (useGain) {
							sceneManager->registerReceiverGain(i,gainIdRx);
						}
						++i;
					}

				} else {
					float radius=this->sphereRadius;
					sceneManager->addReceiver(i, p,polarization, radius, sceneManager->printPower);
					if (useGain) {
						sceneManager->registerReceiverGain(i,gainIdRx);
					}
				}
				++i;
			}

			// Gain use
			if (useGain) {
				sceneManager->registerTransmitterGain(totalRx+1,gainIdTx);
			}

			ResultReport* rp=sceneManager->transmit(totalRx+1, txPower, postx, polarization,freq,false);
			report.merge(*rp);
			sceneManager->clearReceivers();
		}

		std::string savePath=folder+"/"+baseName+"-"+std::to_string(rx.size())+"-"+type+"-"+rand+"-"+std::to_string(sphereRadius)+".csv";
		std::cout<<"Saving results to "<<savePath<<std::endl;
		report.toCSV(savePath);
		timer.stop();
		timerRT.stop();
		std::cout<<"Time: \t"<<timer.getTime()<<"\t"<<timerRT.getTime()<<std::endl;
	}

void LoraDoppler::test() {
		this->replications=1;

		freq=868e6;
		std::string path("test-lora.json");;
		bool useRDN=false;
		useGain=false;
		loadScenarioLora(false,useRDN, path);
	}

