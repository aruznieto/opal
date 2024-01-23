/***************************************************************/
//
//Copyright (c) 2021 Esteban Egea-Lopez http://girtel.upct.es/~simulations
//
/**************************************************************/
#include "lille.h"
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
#include <optixu/optixu_quaternion_namespace.h> 
#include <sstream>

using namespace opal;
using namespace optix;

std::string Basestation::toString() const {
	std::stringstream os;
	os<<id<<"\t"<<antenna_id<<"\t"<<postx<<"\t"<<azimuth<<"\t"<<freq<<"\t"<<freq_end<<"\t"<<site<<"\t"<<tech<<"\t"<<lat<<"\t"<<lon<<"\t"<<antenna_index<<std::endl;
	return os.str();
}
Lille::Lille(OpalSceneManager*   sceneManager, float sphereRadius, bool useDepolarization) : BasicTests(sceneManager, sphereRadius, useDepolarization) {
	this->freq=2e9;
	//this->baseStation= make_float3(365.4184f, 7.0f, 346.6396f);;
	this->baseStation= make_float3(359.674f, 2.10f, 340.9827f);

	this->gsufix=std::to_string(0);
	this->replications=1;
	this->useGain=false;
	this->txPower=0.001f;
	this->delta_rx=0.075f;
	this->delta_tx_x=0.025f;
	this->delta_tx_y=0.025f;
	this->mimotx=8;
	this->mimorx=1;
}
void Lille::runBasestations(std::string test, std::string outputFile) {
	std::vector<int> tokens = parseTestString(test);

	this->replications = 1;

	std::string path;
	path = "emf_5G/waz_scene2.json";
	//this->rx_file = "emf_5G/rx-waz.txt";
	this->rx_file = "emf_5G/waz_50sens.txt";
	std::string txPath="emf_5G/wazemmes_tx_info_unity.csv";
	if (tokens[0]==1) {
		txPath="emf_5G/euratech_tx_info_unity.csv";
	}
	txPower = 1.0f; //1 W


	bool useRDN = false;

	useGain = false;
	if (tokens[1] == 1) {
		useGain = true;
	}
	//true here adds a random displacement to receivers positions
	bool addRandom = false;
	executeBasestation(addRandom, useRDN, path, outputFile, txPath);
	//executeSingleBasestation(addRandom, useRDN, path, outputFile, txPath);
}

void Lille::testAntennaOrientation() {
	//orientateAntennaPattern(14.0,-5.0);
	float freq = 5.9e9f;
	sceneManager->enableMultiChannel();
	sceneManager->setMinEpsilon(1e-3);
	sceneManager->setUseAntennaGain(true);
	//sceneManager->enableGenerateRaysOnLaunch();
	//Use field or induced voltage
	ComputeMode mode = ComputeMode::FIELD;
	LPFlatMeshReflectionSimulation* sim = new LPFlatMeshReflectionSimulation(sceneManager);
	sceneManager->setSimulation(sim);
	sim->setComputeMode(mode);
	sceneManager->initContext(freq);


	optix::float3 polarization = make_float3(0.0f, 1.0f, 0.0f);
	std::string gp("emf_5G/diagram.txt");
	//std::string gainPathT = gp + gsufix + ".txt";
	AntennaGain gains = sceneManager->loadGainsFromFileIndBPower(gp.c_str());
	int gainIdTx = sceneManager->registerAntennaGain(gains);
	optix::float3 posrx = make_float3(10.0f, 10.0f, 0.0f);
	sceneManager->addReceiver(1, posrx, polarization, sphereRadius, sceneManager->printPower);
	sceneManager->registerTransmitterGain(0,gainIdTx);


	optix::float3 postx = make_float3(0.0f, 10.0f, 0.0f);
	float3 mobile=make_float3(5.0,4.0,3.0);
	//The zero azimuth orientation and 90 degrees elevation have to point to the mobile receiver
	float3 zeroAzimuth=normalize(mobile-postx);
	//***Single ray transmit****
	float3 mRay=zeroAzimuth;
	//float3 mRay=normalize(make_float3(1,0,0));
	sceneManager->createRaySphere2D(1,1,&mRay);


	sceneManager->finishSceneContext();
	Matrix<4,4> t=sceneManager->orientateAntennaPattern(90.0,0.0);
	//Test passing axis
	//float3 azimuthZero=make_float3(1,0,0);
	//float3 azimuthRotation = make_float3(0,1,0);
	//Matrix<4,4> t=sceneManager->orientateAntennaPattern(90.0,0.0,azimuthRotation, azimuthZero);


	///Normal to Y (up) and zeroAzimuth: is the normal vector to the plane Y-zeroAzimuth and defines the azimuth direction
//	float3 x=normalize(cross(make_float3(0,1,0),zeroAzimuth));
//	//Now, a normal to zeroAzimuth and x defines the rotationVector for the azimuth
//	float3 azimuthRotation = normalize(cross(zeroAzimuth, x));
//	Matrix<4,4> t=sceneManager->orientateAntennaPattern(0.0,0.0,azimuthRotation, zeroAzimuth);
//
	sceneManager->getActiveTransmitters()[0]->transformToPolarization=t;
	sceneManager->transmit(0, 1.0f, postx, polarization, false,true);
	return;
}

//Example of how to create a directed beam from basestation to mobile phone
void  Lille::executeBasestationBeamforming(bool addRandom, bool useRDN, std::string scenarioPath, std::string outputFile, std::string txPath) {
	Timer timer;
	Timer timerRT;
	std::cout << "Load executeBasestationBeamforming with scenario " << scenarioPath << "; txPath=" << txPath <<  std::endl;
	timer.start();
	//Init context before doing anything else
	//Enable multichannel since each transmitter use different frequencies
	sceneManager->enableMultiChannel();
	sceneManager->setMinEpsilon(1e-3);
	sceneManager->setUseAntennaGain(useGain);
	sceneManager->enableGenerateRaysOnLaunch();
	//Use field or induced voltage
	ComputeMode mode = ComputeMode::FIELD;
	//ComputeMode mode=ComputeMode::VOLTAGE;
	if (useRDN) {
		RayDensityNormalizationSimulation* sim = new RayDensityNormalizationSimulation(sceneManager);
		sceneManager->setSimulation(sim);
		sim->setComputeMode(mode);
	}
	else {
		if (useDepolarization) {
			LPFlatMeshReflectionSimulation* sim = new LPFlatMeshReflectionSimulation(sceneManager);
			sceneManager->setSimulation(sim);
			sim->setComputeMode(mode);
			//Create log trace (only for one tx and one rx)
			//sim->setEnableTraceLog(false);
			//sim->setEnableSimulation(false);
			//sim->setPrintHits(true);
		}
		else {
			BasicFlatMeshReflectionSimulation* sim = new BasicFlatMeshReflectionSimulation(sceneManager);
			sceneManager->setSimulation(sim);
			sim->setComputeMode(mode);
			//sim->setPrintHits(true);
		}
	}
	//Add diffraction
	SingleDiffraction* simd = new SingleDiffraction(sceneManager);
	sceneManager->setSimulation(simd);
	simd->setComputeMode(mode);
	//Create log trace (only for one tx and one rx)
	//simd->setEnableTraceLog(false);
	//simd->setPrintHits(true);
	//Disable or enable diffraction
	//simd->setEnableSimulation(false);

	sceneManager->initContext(freq);
	//Exceptions
	//sceneManager->enableExceptions();	

	//Load files here

	std::vector<Basestation> tx=loadTransmittersFromFile(txPath);
	for (auto bs : tx) {
		std::cout<<bs.toString();
	}
	int num_tx=tx.size();
	std::cout<<"Loaded a total of "<<num_tx<<" transmitters"<<std::endl;
	ScenarioLoader* loader = new ScenarioLoader(sceneManager);
	loader->loadJSONScenario(scenarioPath);
	std::vector<int> gainIdTx(7);
	int gainIdRx = -1;
	if (useGain) {
		for (int i=0; i<7; ++i) {
			std::string gainPathT ="emf_5G/diagram"  + std::to_string(i) + ".txt";
			std::string gp(gainPathT);
			AntennaGain gains = sceneManager->loadGainsFromFileIndBPower(gp.c_str());
			gainIdTx[i] = sceneManager->registerAntennaGain(gains);
		}
		//No gain for receivers. Uncomment if necessary
		//Receiver gain
		//std::string gpR("emf_5G/diagram.txt");
		//std::string gainPathR = gpR + gsufix + ".txt";
		//AntennaGain gainsRx = sceneManager->loadGainsFromFileIndBPower(gpR.c_str());
		// gainIdRx = sceneManager->registerAntennaGain(gainsRx);
	}

	optix::float3 polarization = make_float3(0.0f, 1.0f, 0.0f);

	if (useRDN) {
		sceneManager->finishSceneContext();
		int rayD = 10000;
		sceneManager->setRayRange(0.0, 180.0, 0.0, 360.0, rayD, rayD);
		RayDensityNormalizationSimulation* s = dynamic_cast<RayDensityNormalizationSimulation*>(sceneManager->getSimulation(0));
		s->setInitialDensity(((float)sceneManager->getRaySphere().rayCount) / (4 * M_PIf));
		s->setFiltering(2u);

	}
	else {
		sceneManager->finishSceneContext();
		sceneManager->createRaySphere2D(0.0f, 0.1, 180.0f, 0.0f, 0.1, 360.0f);
	}

	std::vector<float3> rx = loadReceiversFromFile(rx_file);
	Basestation aux=tx[tx.size()-1];
	int i = aux.id+1;
	int maxReceivers = 500;
	int j = 0;
	int maxReportIndex=50;
	int reportIndex=0;




	//Multiple transmitters or launches. full report
	ResultReport report;
	if (maxReceivers >=rx.size()) {
		for (auto p : rx) {

			float3 posrx = make_float3(p.x, p.y, p.z);

			sceneManager->addReceiver(i, posrx, polarization, sphereRadius, sceneManager->printPower);
			//No gain for receivers. Uncomment if necessary
			//if (useGain) {
			//	sceneManager->registerReceiverGain(i, gainIdRx);
			//}
			++i;
		}
	}


	//Create a mobile receiver (it is not actually registered as receiver, but it could be registered, with that position)
	float3 mobile = make_float3(773.624, 1.55, 1500.55); 

	//Loop over all basestations
	for (auto bs : tx) {
		//Loop over all receivers in case we have to split them due to memory constraints
		if (maxReceivers <rx.size()) {
			for (auto p : rx) {

				float3 posrx = make_float3(p.x, p.y, p.z);

				sceneManager->addReceiver(i, posrx, polarization, sphereRadius, sceneManager->printPower);
				//No gain for receivers. Uncomment if necessary
				//if (useGain) {
				//	sceneManager->registerReceiverGain(i, gainIdRx);
				//}
				++i;
				++j;
				if (j == maxReceivers) {

					//Transmitter position
					optix::float3 postx = bs.postx;

					timerRT.start();
					if (useGain) {
						sceneManager->registerTransmitterGain(bs.id, gainIdTx[bs.antenna_index]);
					}
					//Orientate antenna pattern
					//The zero azimuth orientation and 90 degrees elevation have to point to the mobile receiver
//					float3 zeroAzimuth=normalize(mobile-postx);
//					//x=Normal to Y (up) and zeroAzimuth: is the normal vector to the plane Y-zeroAzimuth. 
//					//The plane zeroAzimuth-x defines the azimuth direction
//					float3 x=normalize(cross(make_float3(0,1,0),zeroAzimuth));
//					//Now, a normal to zeroAzimuth and x defines the rotationVector for the azimuth
//					float3 azimuthRotation = normalize(cross(zeroAzimuth, x));
//
					//Matrix<4,4> t=sceneManager->orientateAntennaPattern(0.0,0.0,azimuthRotation, zeroAzimuth);
					Matrix<4,4> t=sceneManager->pointAntennaPatternTo(postx, mobile);
					sceneManager->getActiveTransmitters()[0]->transformToPolarization=t;
					ResultReport* rp = sceneManager->transmit(bs.id, txPower, postx, polarization, bs.freq, false, true);
					//Multiple transmitters. Merge report
					report.merge(*rp);
					timerRT.stop();
					std::cout <<"Launch time="<< timerRT.getTime() << std::endl;
					//Start for new receivers
					j = 0;
					sceneManager->clearReceivers();
				}
			}
			//rp->toCSV(outputFile);


		}

		optix::float3 postx = bs.postx;
		timerRT.start();
		if (useGain) {
			sceneManager->registerTransmitterGain(bs.id, gainIdTx[bs.antenna_index]);
		}
		//Orientate antenna pattern
		//The zero azimuth orientation and 90 degrees elevation have to point to the mobile receiver
//		float3 zeroAzimuth=normalize(mobile-postx);
//		//x=Normal to Y (up) and zeroAzimuth: is the normal vector to the plane Y-zeroAzimuth. 
//		//The plane zeroAzimuth-x defines the azimuth direction
//		float3 x=normalize(cross(make_float3(0,1,0),zeroAzimuth));
//		//Now, a normal to zeroAzimuth and x defines the rotationVector for the azimuth
//		float3 azimuthRotation = normalize(cross(zeroAzimuth, x));
//
		//Matrix<4,4> t=sceneManager->orientateAntennaPattern(0.0,0.0,azimuthRotation, zeroAzimuth);
		Matrix<4,4> t=sceneManager->pointAntennaPatternTo(postx,mobile);
		sceneManager->getActiveTransmitters()[0]->transformToPolarization=t;
		ResultReport* rp = sceneManager->transmit(bs.id, txPower, postx, polarization, bs.freq, false, true);
		//Multiple transmitters. Merge report
		report.merge(*rp);

		++reportIndex;
		if (reportIndex==maxReportIndex) {
			std::string rpath=outputFile+"-"+std::to_string(bs.id)+".csv";
			report.toCSV(rpath);
			reportIndex=0;
			report.clear();
		}

		timerRT.stop();
		std::cout << "Launch time=" << timerRT.getTime() << std::endl;
		if (maxReceivers <rx.size()) {
			i = aux.id+1;
			//Start for new receivers
			j = 0;
			sceneManager->clearReceivers();
		} else {
			//No need to clear the receivers
		}

	}


	//Multiple transmitters. save full report csv
	report.toCSV(outputFile);
	timer.stop();

	std::cout << "Total time\t" << timer.getTime() << "\t" << timerRT.getTime() << std::endl;
	delete loader;
}
void  Lille::executeSingleBasestation(bool addRandom, bool useRDN, std::string scenarioPath, std::string outputFile, std::string txPath) {
	Timer timer;
	Timer timerRT;
	std::cout << "Load executeBasestation with scenario " << scenarioPath << "; txPath=" << txPath <<  std::endl;
	timer.start();
	//Init context before doing anything else
	//Enable multichannel since each transmitter use different frequencies
	sceneManager->enableMultiChannel();
	sceneManager->setMinEpsilon(1e-3);
	sceneManager->setUseAntennaGain(useGain);
	sceneManager->enableGenerateRaysOnLaunch();
	//Use field or induced voltage
	ComputeMode mode = ComputeMode::FIELD;
	//ComputeMode mode=ComputeMode::VOLTAGE;
	if (useRDN) {
		RayDensityNormalizationSimulation* sim = new RayDensityNormalizationSimulation(sceneManager);
		sceneManager->setSimulation(sim);
		sim->setComputeMode(mode);
	}
	else {
		if (useDepolarization) {
			LPFlatMeshReflectionSimulation* sim = new LPFlatMeshReflectionSimulation(sceneManager);
			sceneManager->setSimulation(sim);
			sim->setComputeMode(mode);
			//sim->setEnableSimulation(false);
			//Create log trace (only for one tx and one rx)
			//sim->setEnableTraceLog(true);
			//sim->setPrintHits(true);
		}
		else {
			BasicFlatMeshReflectionSimulation* sim = new BasicFlatMeshReflectionSimulation(sceneManager);
			sceneManager->setSimulation(sim);
			sim->setComputeMode(mode);
			//sim->setPrintHits(true);
		}
	}
	//Add diffraction
	SingleDiffraction* simd = new SingleDiffraction(sceneManager);
	sceneManager->setSimulation(simd);
	simd->setComputeMode(mode);
	//Create log trace (only for one tx and one rx)
	//simd->setEnableTraceLog(true);
	//simd->setPrintHits(true);
	//Disable or enable diffraction
	//simd->setEnableSimulation(false);

	sceneManager->initContext(freq);
	//Exceptions
	//sceneManager->enableExceptions();	

	//Load files here

	std::vector<Basestation> tx=loadTransmittersFromFile(txPath);
	for (auto bs : tx) {
		std::cout<<bs.toString();
	}
	int num_tx=tx.size();
	std::cout<<"Loaded a total of "<<num_tx<<" transmitters"<<std::endl;
	ScenarioLoader* loader = new ScenarioLoader(sceneManager);
	loader->loadJSONScenario(scenarioPath);
	std::vector<int> gainIdTx(7);
	int gainIdRx = -1;
	if (useGain) {
		for (int i=0; i<7; ++i) {
			std::string gainPathT ="emf_5G/diagram"  + std::to_string(i) + ".txt";
			std::string gp(gainPathT);
			AntennaGain gains = sceneManager->loadGainsFromFileIndBPower(gp.c_str());
			gainIdTx[i] = sceneManager->registerAntennaGain(gains);
		}
		//No gain for receivers. Uncomment if necessary
		//Receiver gain
		//std::string gpR("emf_5G/diagram.txt");
		//std::string gainPathR = gpR + gsufix + ".txt";
		//AntennaGain gainsRx = sceneManager->loadGainsFromFileIndBPower(gpR.c_str());
		// gainIdRx = sceneManager->registerAntennaGain(gainsRx);
	}

	optix::float3 polarization = make_float3(0.0f, 1.0f, 0.0f);

	if (useRDN) {
		sceneManager->finishSceneContext();
		int rayD = 10000;
		sceneManager->setRayRange(0.0, 180.0, 0.0, 360.0, rayD, rayD);
		RayDensityNormalizationSimulation* s = dynamic_cast<RayDensityNormalizationSimulation*>(sceneManager->getSimulation(0));
		s->setInitialDensity(((float)sceneManager->getRaySphere().rayCount) / (4 * M_PIf));
		s->setFiltering(2u);

	}
	else {
		sceneManager->finishSceneContext();
		sceneManager->createRaySphere2D(0.0f, 0.1, 180.0f, 0.0f, 0.1, 360.0f);
	}

	std::vector<float3> rx = loadReceiversFromFile(rx_file);
	Basestation aux=tx[tx.size()-1];
	int i = aux.id+1;
	int maxReceivers = 500;
	int j = 0;
	int maxReportIndex=50;
	int reportIndex=0;
	//Multiple transmitters or launches. full report
	ResultReport report;
	float3 p=rx[1];
	float3 posrx = make_float3(p.x, p.y, p.z);

	sceneManager->addReceiver(i, posrx, polarization, sphereRadius, sceneManager->printPower);
	//No gain for receivers. Uncomment if necessary
	//if (useGain) {
	//	sceneManager->registerReceiverGain(i, gainIdRx);
	//}
	++i;
	//Loop over all receivers in case we have to split them due to memory constraints
	timerRT.start();

	Basestation bs=tx[0];
	optix::float3 postx = bs.postx;
	timerRT.start();
	std::cout<<"antenna index"<<bs.antenna_index<<std::endl;
	if (useGain) {
		sceneManager->registerTransmitterGain(bs.id, gainIdTx[bs.antenna_index]);
	}
	//Orientate antenna pattern
	//No tilt, we used antenna index up to 6, only 7 had a tilt of -5 (upward)
	Matrix<4,4> t=sceneManager->orientateAntennaPattern(bs.azimuth,0.0);
	sceneManager->getActiveTransmitters()[0]->transformToPolarization=t;
	ResultReport* rp = sceneManager->transmit(bs.id, txPower, postx, polarization, bs.freq, false, true);
	//ResultReport* rp = sceneManager->transmit(rx.size() + 1, txPower, postx, polarization);
	//Multiple transmitters. Merge report
	report.merge(*rp);

	timerRT.stop();
	std::cout << "Launch time=" << timerRT.getTime() << std::endl;




	//Multiple transmitters. save full report csv
	report.toCSV(outputFile);
	timer.stop();

	std::cout << "Total time\t" << timer.getTime() << "\t" << timerRT.getTime() << std::endl;
	delete loader;
}

void  Lille::executeBasestation(bool addRandom, bool useRDN, std::string scenarioPath, std::string outputFile, std::string txPath) {
	Timer timer;
	Timer timerRT;
	std::cout << "Load executeBasestation with scenario " << scenarioPath << "; txPath=" << txPath <<  std::endl;
	timer.start();
	//Init context before doing anything else
	//Enable multichannel since each transmitter use different frequencies
	sceneManager->enableMultiChannel();
	sceneManager->setMinEpsilon(1e-3);
	sceneManager->setUseAntennaGain(useGain);
	sceneManager->enableGenerateRaysOnLaunch();
	//Use field or induced voltage
	ComputeMode mode = ComputeMode::FIELD;
	//ComputeMode mode=ComputeMode::VOLTAGE;
	if (useRDN) {
		RayDensityNormalizationSimulation* sim = new RayDensityNormalizationSimulation(sceneManager);
		sceneManager->setSimulation(sim);
		sim->setComputeMode(mode);
	}
	else {
		if (useDepolarization) {
			LPFlatMeshReflectionSimulation* sim = new LPFlatMeshReflectionSimulation(sceneManager);
			sceneManager->setSimulation(sim);
			sim->setComputeMode(mode);
			//Create log trace (only for one tx and one rx)
			//sim->setEnableTraceLog(false);
			//sim->setEnableSimulation(false);
			//sim->setPrintHits(true);
		}
		else {
			BasicFlatMeshReflectionSimulation* sim = new BasicFlatMeshReflectionSimulation(sceneManager);
			sceneManager->setSimulation(sim);
			sim->setComputeMode(mode);
			//sim->setPrintHits(true);
		}
	}
	//Add diffraction
	SingleDiffraction* simd = new SingleDiffraction(sceneManager);
	sceneManager->setSimulation(simd);
	simd->setComputeMode(mode);
	//Create log trace (only for one tx and one rx)
	//simd->setEnableTraceLog(false);
	//simd->setPrintHits(true);
	//Disable or enable diffraction
	//simd->setEnableSimulation(false);

	sceneManager->initContext(freq);
	//Exceptions
	//sceneManager->enableExceptions();	

	//Load files here

	std::vector<Basestation> tx=loadTransmittersFromFile(txPath);
	for (auto bs : tx) {
		std::cout<<bs.toString();
	}
	int num_tx=tx.size();
	std::cout<<"Loaded a total of "<<num_tx<<" transmitters"<<std::endl;
	ScenarioLoader* loader = new ScenarioLoader(sceneManager);
	loader->loadJSONScenario(scenarioPath);
	std::vector<int> gainIdTx(7);
	int gainIdRx = -1;
	if (useGain) {
		for (int i=0; i<7; ++i) {
			std::string gainPathT ="emf_5G/diagram"  + std::to_string(i) + ".txt";
			std::string gp(gainPathT);
			AntennaGain gains = sceneManager->loadGainsFromFileIndBPower(gp.c_str());
			gainIdTx[i] = sceneManager->registerAntennaGain(gains);
		}
		//No gain for receivers. Uncomment if necessary
		//Receiver gain
		//std::string gpR("emf_5G/diagram.txt");
		//std::string gainPathR = gpR + gsufix + ".txt";
		//AntennaGain gainsRx = sceneManager->loadGainsFromFileIndBPower(gpR.c_str());
		// gainIdRx = sceneManager->registerAntennaGain(gainsRx);
	}

	optix::float3 polarization = make_float3(0.0f, 1.0f, 0.0f);

	if (useRDN) {
		sceneManager->finishSceneContext();
		int rayD = 10000;
		sceneManager->setRayRange(0.0, 180.0, 0.0, 360.0, rayD, rayD);
		RayDensityNormalizationSimulation* s = dynamic_cast<RayDensityNormalizationSimulation*>(sceneManager->getSimulation(0));
		s->setInitialDensity(((float)sceneManager->getRaySphere().rayCount) / (4 * M_PIf));
		s->setFiltering(2u);

	}
	else {
		sceneManager->finishSceneContext();
		sceneManager->createRaySphere2D(0.0f, 0.1, 180.0f, 0.0f, 0.1, 360.0f);
	}

	std::vector<float3> rx = loadReceiversFromFile(rx_file);
	Basestation aux=tx[tx.size()-1];
	int i = aux.id+1;
	int maxReceivers = 500;
	int j = 0;
	int maxReportIndex=50;
	int reportIndex=0;
	//Multiple transmitters or launches. full report
	ResultReport report;
	if (maxReceivers >=rx.size()) {
		for (auto p : rx) {

			float3 posrx = make_float3(p.x, p.y, p.z);

			sceneManager->addReceiver(i, posrx, polarization, sphereRadius, sceneManager->printPower);
			//No gain for receivers. Uncomment if necessary
			//if (useGain) {
			//	sceneManager->registerReceiverGain(i, gainIdRx);
			//}
			++i;
		}
	}
	//Loop over all basestations
	for (auto bs : tx) {
		//Loop over all receivers in case we have to split them due to memory constraints
		timerRT.start();
		if (maxReceivers <rx.size()) {
			for (auto p : rx) {

				float3 posrx = make_float3(p.x, p.y, p.z);

				sceneManager->addReceiver(i, posrx, polarization, sphereRadius, sceneManager->printPower);
				//No gain for receivers. Uncomment if necessary
				//if (useGain) {
				//	sceneManager->registerReceiverGain(i, gainIdRx);
				//}
				++i;
				++j;
				if (j == maxReceivers) {

					//Transmitter position
					optix::float3 postx = bs.postx;

					timerRT.start();
					if (useGain) {
						sceneManager->registerTransmitterGain(bs.id, gainIdTx[bs.antenna_index]);
					}
					//Orientate antenna pattern
					//No tilt, we used antenna index up to 6, only 7 had a tilt of -5 (upward)
					Matrix<4,4> t=sceneManager->orientateAntennaPattern(bs.azimuth,0.0);
					sceneManager->getActiveTransmitters()[0]->transformToPolarization=t;
					ResultReport* rp = sceneManager->transmit(bs.id, txPower, postx, polarization, bs.freq, false, true);
					//Multiple transmitters. Merge report
					report.merge(*rp);
					timerRT.stop();
					std::cout <<"Launch time="<< timerRT.getTime() << std::endl;
					//Start for new receivers
					j = 0;
					sceneManager->clearReceivers();
				}
			}
			//rp->toCSV(outputFile);


		}

		optix::float3 postx = bs.postx;
		timerRT.start();
		if (useGain) {
			sceneManager->registerTransmitterGain(bs.id, gainIdTx[bs.antenna_index]);
		}
		//Orientate antenna pattern
		//No tilt, we used antenna index up to 6, only 7 had a tilt of -5 (upward)
		Matrix<4,4> t=sceneManager->orientateAntennaPattern(bs.azimuth,0.0);
		sceneManager->getActiveTransmitters()[0]->transformToPolarization=t;
		ResultReport* rp = sceneManager->transmit(bs.id, txPower, postx, polarization, bs.freq, false, true);
		//Multiple transmitters. Merge report
		report.merge(*rp);

		++reportIndex;
		if (reportIndex==maxReportIndex) {
			std::string rpath=outputFile+"-"+std::to_string(bs.id)+".csv";
			report.toCSV(rpath);
			reportIndex=0;
			report.clear();
		}

		timerRT.stop();
		std::cout << "Launch time=" << timerRT.getTime() << std::endl;
		if (maxReceivers <rx.size()) {
			i = aux.id+1;
			//Start for new receivers
			j = 0;
			sceneManager->clearReceivers();
		} else {
			//No need to clear the receivers
		}

	}


	//Multiple transmitters. save full report csv
	std::string rpath=outputFile+"-last.csv";
	report.toCSV(rpath);
	timer.stop();

	std::cout << "Total time\t" << timer.getTime() << "\t" << timerRT.getTime() << std::endl;
	delete loader;
}

void Lille::runWaz(std::string test, std::string outputFile) {
	std::vector<float> frequencies;
	frequencies.push_back(1860e6);
	frequencies.push_back(791e6);
	//From 5.85 GHz to 5.93, 818 point
	/*for (int i = 0; i < 818; i++) {
	  float f = 5.85e9 + i * 9.78e4;
	  frequencies.push_back(f);

	  }*/
	std::vector<float3> txPositions;
	txPositions.push_back(make_float3(1050.638f, 15, 161.0417f));
	txPositions.push_back(make_float3(989.638f, 15, 365.0417f));
	std::vector<int> tokens = parseTestString(test);

	this->replications = 1;

	freq = frequencies[tokens[0]];
	baseStation = txPositions[tokens[0]];
	std::string path;
	path = "waz_scene.json";
	this->rx_file = "rx-waz.txt";
	txPower = 1.0f; //1 W


	bool useRDN = false;

	useGain = false;
	if (tokens[1] == 1) {
		useGain = true;
	}
	//true here adds a random displacement to receivers positions
	bool addRandom = false;
	loadWazem(addRandom, useRDN, path, outputFile);
}

void  Lille::loadWazem(bool addRandom, bool useRDN, std::string path, std::string outputFile) {
	Timer timer;
	Timer timerRT;
	std::cout << "Load Wazem with scenario " << path << "; freq=" << freq << ";tx=" << baseStation << std::endl;
	timer.start();
	//Init context before doing anything else
	//sceneManager->enableGenerateRaysOnLaunch();
	sceneManager->setMinEpsilon(1e-3);
	sceneManager->setUseAntennaGain(useGain);
	sceneManager->enableGenerateRaysOnLaunch();
	//Use field or induced voltage
	ComputeMode mode = ComputeMode::FIELD;
	//ComputeMode mode=ComputeMode::VOLTAGE;
	if (useRDN) {
		RayDensityNormalizationSimulation* sim = new RayDensityNormalizationSimulation(sceneManager);
		sceneManager->setSimulation(sim);
		sim->setComputeMode(mode);
	}
	else {
		if (useDepolarization) {
			LPFlatMeshReflectionSimulation* sim = new LPFlatMeshReflectionSimulation(sceneManager);
			sceneManager->setSimulation(sim);
			sim->setComputeMode(mode);
			//Create log trace (only for one tx and one rx)
			//sim->setEnableTraceLog(false);
			//sim->setEnableSimulation(false);
			//sim->setPrintHits(true);
		}
		else {
			BasicFlatMeshReflectionSimulation* sim = new BasicFlatMeshReflectionSimulation(sceneManager);
			sceneManager->setSimulation(sim);
			sim->setComputeMode(mode);
			//sim->setPrintHits(true);
		}
	}
	//Add diffraction
	SingleDiffraction* simd = new SingleDiffraction(sceneManager);
	sceneManager->setSimulation(simd);
	simd->setComputeMode(mode);
	//Create log trace (only for one tx and one rx)
	//simd->setEnableTraceLog(false);
	//simd->setPrintHits(true);
	//Disable or enable diffraction
	//simd->setEnableSimulation(false);

	sceneManager->initContext(freq);
	//Exceptions
	//sceneManager->enableExceptions();	

	//Load files here
	ScenarioLoader* loader = new ScenarioLoader(sceneManager);
	loader->loadJSONScenario(path);
	int gainIdTx = -1;
	int gainIdRx = -1;
	if (useGain) {
		//Change antenna paths	
		std::string gp("dipole.txt");
		//std::string gainPathT = gp + gsufix + ".txt";
		AntennaGain gains = sceneManager->loadGainsFromFileIndBPower(gp.c_str());
		int gainIdTx = sceneManager->registerAntennaGain(gains);
		//Receiver gain
		std::string gpR("dipole.txt");
		//std::string gainPathR = gpR + gsufix + ".txt";
		AntennaGain gainsRx = sceneManager->loadGainsFromFileIndBPower(gpR.c_str());
		int gainIdRx = sceneManager->registerAntennaGain(gainsRx);
	}

	optix::float3 polarization = make_float3(0.0f, 1.0f, 0.0f);

	if (useRDN) {
		sceneManager->finishSceneContext();
		int rayD = 10000;
		sceneManager->setRayRange(0.0, 180.0, 0.0, 360.0, rayD, rayD);
		RayDensityNormalizationSimulation* s = dynamic_cast<RayDensityNormalizationSimulation*>(sceneManager->getSimulation(0));
		s->setInitialDensity(((float)sceneManager->getRaySphere().rayCount) / (4 * M_PIf));
		s->setFiltering(2u);

	}
	else {
		sceneManager->finishSceneContext();
		sceneManager->createRaySphere2D(0.0f, 0.1, 180.0f, 0.0f, 0.1, 360.0f);
	}

	std::vector<float3> rx = loadReceiversFromFile(rx_file);
	int i = 0;
	int maxReceivers = 6000;
	int j = 0;
	//Multiple transmitters or launches. full report
	ResultReport report;
	for (auto p : rx) {

		float3 posrx = make_float3(p.x, p.y, p.z);

		sceneManager->addReceiver(i, posrx, polarization, sphereRadius, sceneManager->printPower);
		if (useGain) {
			sceneManager->registerReceiverGain(i, gainIdRx);
		}
		++i;
		++j;
		if (j == maxReceivers) {

			//Transmitter position
			optix::float3 postx = baseStation;
			timerRT.start();
			if (useGain) {
				sceneManager->registerTransmitterGain(rx.size() + 1, gainIdTx);
			}
			ResultReport* rp = sceneManager->transmit(rx.size() + 1, txPower, postx, polarization);
			//Multiple transmitters. Merge report
			report.merge(*rp);

			//Multiple transmitters. New transmitter
			/*postx = make_float3(1024.0, 47.5, 512.1);
			  if (useGain) {
			  sceneManager->registerTransmitterGain(rx.size() + 2, gainIdTx);
			  }
			  rp = sceneManager->transmit(rx.size() + 2, txPower, postx, polarization);
			  report.merge(*rp);
			  */

			timerRT.stop();
			std::cout <<"Launch time"<< timerRT.getTime() << std::endl;
			//Start for new receivers
			j = 0;
			sceneManager->clearReceivers();
		}
		//rp->toCSV(outputFile);


	}

	optix::float3 postx = baseStation;
	timerRT.start();
	if (useGain) {
		sceneManager->registerTransmitterGain(rx.size() + 1, gainIdTx);
	}
	ResultReport* rp = sceneManager->transmit(rx.size() + 1, txPower, postx, polarization);
	//Multiple transmitters. Merge report
	report.merge(*rp);

	//Multiple transmitters. New transmitter
	/*postx = make_float3(1024.0, 47.5, 512.1);
	  if (useGain) {
	  sceneManager->registerTransmitterGain(i + 2, gainIdTx);
	  }
	  rp = sceneManager->transmit(i + 2, txPower, postx, polarization);
	  report.merge(*rp);
	  */

	timerRT.stop();
	std::cout << "Launch time" << timerRT.getTime() << std::endl;
	//Start for new receivers
	j = 0;
	sceneManager->clearReceivers();


	//Multiple transmitters. save full report csv
	report.toCSV(outputFile);
	timer.stop();

	std::cout << "Total time\t" << timer.getTime() << "\t" << timerRT.getTime() << std::endl;
	delete loader;
}

//Callback functions
void printAoD(int txId, int rxId, HitRecord info) {
	//std::cout<<std::setprecision(15)<<"AoD\t"<<txId<<"\t"<<rxId<<"\t"<<E.y<<"\t"<<doad.x<<"\t"<<doad.y<<"\t"<<doad.z<<"\t"<<doad.w<<"\t"<<receivers[index]->externalId<<"\t"<< tx->externalId<<"\t"<<dod.x<<"\t"<<dod.y<<"\t"<<dod.z<<std::endl;
	float2 E = info.E;
	
	float3 dod=info.directionOfDeparture;
	float3 doad=info.directionOfArrival;
	//Departure
	float2 a=GeometryUtils::getAngles(dod);
	if (a.y<0) {
		a.y=(2*M_PIf) +a.y;
	}
	//Now, we have to remove 45 degrees to the azimuth, because the antenna was actually rotated 45 degrees with respect to the global X (see below)
	//So we get the azimuth in the antenna local frame
	a.y=a.y-(M_PIf/4);
	
	float2 b=GeometryUtils::getAngles(doad);
	//std::cout<<"AoD,"<<txId<<","<<rxId<<","<<E.x<<","<<E.y<<","<<a.x<<","<<a.y<<","<<b.x<<","<<b.y<<","<<doad.w<<std::endl;
}
//void printAoD(int txId, int rxId, HitRecord info) {
//	//std::cout<<std::setprecision(15)<<"AoD\t"<<txId<<"\t"<<rxId<<"\t"<<E.y<<"\t"<<doad.x<<"\t"<<doad.y<<"\t"<<doad.z<<"\t"<<doad.w<<"\t"<<receivers[index]->externalId<<"\t"<< tx->externalId<<"\t"<<dod.x<<"\t"<<dod.y<<"\t"<<dod.z<<std::endl;
//	float2 E = make_float2(info->EEx.x,info->EEx.y);
//	
//	float4 dod=info->rdud;
//	float4 doad=info->doaD;
//	//Departure
//	float2 a=GeometryUtils::getAngles(make_float3(dod.x,dod.y,dod.z));
//	if (a.y<0) {
//		a.y=(2*M_PIf) +a.y;
//	}
//	//Now, we have to remove 45 degrees to the azimuth, because the antenna was actually rotated 45 degrees with respect to the global X (see below)
//	//So we get the azimuth in the antenna local frame
//	a.y=a.y-(M_PIf/4);
//	
//	float2 b=GeometryUtils::getAngles(make_float3(doad.x,doad.y,doad.z));
//	//std::cout<<"AoD,"<<txId<<","<<rxId<<","<<E.x<<","<<E.y<<","<<a.x<<","<<a.y<<","<<b.x<<","<<b.y<<","<<doad.w<<std::endl;
//}
void printAoDDiffraction(int txId, int rxId, HitRecord info) {
	//std::cout<<std::setprecision(15)<<"AoD\t"<<txId<<"\t"<<rxId<<"\t"<<E.y<<"\t"<<doad.x<<"\t"<<doad.y<<"\t"<<doad.z<<"\t"<<doad.w<<"\t"<<receivers[index]->externalId<<"\t"<< tx->externalId<<"\t"<<dod.x<<"\t"<<dod.y<<"\t"<<dod.z<<std::endl;
	float2 E = info.E;
	
	float3 dod=info.directionOfDeparture;
	float3 doad=info.directionOfArrival;
	//Departure
	float2 a=GeometryUtils::getAngles(dod);
	if (a.y<0) {
		a.y=(2*M_PIf) +a.y;
	}
	//Now, we have to remove 45 degrees to the azimuth, because the antenna was actually rotated 45 degrees with respect to the global X (see below)
	//So we get the azimuth in the antenna local frame
	a.y=a.y-(M_PIf/4);
	
	float2 b=GeometryUtils::getAngles(doad);
	//std::cout<<"AoDi,"<<txId<<","<<rxId<<","<<E.x<<","<<E.y<<","<<a.x<<","<<a.y<<","<<b.x<<","<<b.y<<","<<doad.w<<std::endl;
}
//void printAoDDiffraction(int txId, int rxId, RDNHit* info) {
//	//std::cout<<std::setprecision(15)<<"AoD\t"<<txId<<"\t"<<rxId<<"\t"<<E.y<<"\t"<<doad.x<<"\t"<<doad.y<<"\t"<<doad.z<<"\t"<<doad.w<<"\t"<<receivers[index]->externalId<<"\t"<< tx->externalId<<"\t"<<dod.x<<"\t"<<dod.y<<"\t"<<dod.z<<std::endl;
//	float2 E = make_float2(info->EEx.x,info->EEx.y);
//	
//	float4 dod=info->doDu;
//	float4 doad=info->doaD;
//	//Departure
//	float2 a=GeometryUtils::getAngles(make_float3(dod.x,dod.y,dod.z));
//	if (a.y<0) {
//		a.y=(2*M_PIf) +a.y;
//	}
//	//Now, we have to remove 45 degrees to the azimuth, because the antenna was actually rotated 45 degrees with respect to the global X (see below)
//	//So we get the azimuth in the antenna local frame
//	a.y=a.y-(M_PIf/4);
//	
//	float2 b=GeometryUtils::getAngles(make_float3(doad.x,doad.y,doad.z));
//	//std::cout<<"AoDi,"<<txId<<","<<rxId<<","<<E.x<<","<<E.y<<","<<a.x<<","<<a.y<<","<<b.x<<","<<b.y<<","<<doad.w<<std::endl;
//}
void  Lille::angleOfDeparture(bool addRandom, bool useRDN, std::string path, std::string outputFile) {


//	float3 ray=make_float3(-1.0,1.0,0.0);
//	float rad2deg=180.0f/M_PIf;
//	//
//	const float3 up=make_float3(0.0,1.0,0.0);
//	const Quaternion a(up,45.0);
//	float3 rr=a*ray;
//	float2 an=GeometryUtils::getAngles(normalize(ray))*rad2deg;
//	float2 anr=GeometryUtils::getAngles(normalize(rr))*rad2deg;
//	if (an.y<0) {
//		an.y=360 +an.y;
//	}
//	if (anr.y<0) {
//		anr.y=360 +anr.y;
//	}
//	std::cout<<"ray"<<ray<<"rr="<<rr<<"an="<<an<<"anr="<<anr<<std::endl;
//	return;
//
	Timer timer;
	Timer timerRT;
	ResultReport report;
	std::cout << "Load Lille with scenario " << path << "; freq=" << freq << ";tx=" << baseStation << std::endl;
	timer.start();
	std::vector<float> frequencies;
	//From 5.85 GHz to 5.93, 818 point
	for (int i = 0; i < 818; i++) {
		float f = 5.85e9 + i * 9.78e4;
		frequencies.push_back(f);

	}
	//Init context before doing anything else
	//Enable multichannel since each transmitter use different frequencies
	sceneManager->enableMultiChannel();
	sceneManager->setMinEpsilon(1e-3);
	sceneManager->setUseAntennaGain(useGain);
	sceneManager->enableGenerateRaysOnLaunch();
	//Use field or induced voltage
	//ComputeMode mode = ComputeMode::FIELD;
	ComputeMode mode = ComputeMode::VOLTAGE;
	if (useRDN) {
		RayDensityNormalizationSimulation* sim = new RayDensityNormalizationSimulation(sceneManager);
		sceneManager->setSimulation(sim);
		sim->setComputeMode(mode);
	}
	else {
		if (useDepolarization) {
			LPFlatMeshReflectionSimulation* sim = new LPFlatMeshReflectionSimulation(sceneManager);
			sceneManager->setSimulation(sim);
			sim->setComputeMode(mode);
			//To get the angles
			sim->setCallback(printAoD);
			//Create log trace (only for one tx and one rx)
			//sim->setEnableTraceLog(false);
			//sim->setEnableSimulation(false);
			//sim->setPrintHits(true);
		}
		else {
			BasicFlatMeshReflectionSimulation* sim = new BasicFlatMeshReflectionSimulation(sceneManager);
			sceneManager->setSimulation(sim);
			sim->setComputeMode(mode);
			sim->setCallback(printAoD);
			//sim->setPrintHits(true);
		}
	}
	//Add diffraction
	SingleDiffraction* simd = new SingleDiffraction(sceneManager);
	sceneManager->setSimulation(simd);
	simd->setComputeMode(mode);
	simd->setCallback(printAoDDiffraction);
	//simd->setAdjustToMemorySize(true);
	//simd->setMaxMemoryFraction(0.1);
	//Create log trace (only for one tx and one rx)
	//simd->setEnableTraceLog(false);
	//simd->setPrintHits(true);
	//Disable or enable diffraction
	//simd->setEnableSimulation(false);

	sceneManager->initContext(freq);
	//Exceptions
	//sceneManager->enableExceptions();	

	//Load files here
	ScenarioLoader* loader = new ScenarioLoader(sceneManager);
	loader->loadJSONScenario(path);
	int gainIdTx = -1;
	int gainIdRx = -1;
	if (useGain) {
		//Change antenna paths	
		std::string gp("lille/gainT");
		std::string gainPathT = gp + gsufix + ".txt";
		AntennaGain gains = sceneManager->loadGainsFromFileIndBPower(gainPathT.c_str());
		int gainIdTx = sceneManager->registerAntennaGain(gains);
		//Receiver gain
		std::string gpR("lille/gainR");
		std::string gainPathR = gpR + gsufix + ".txt";
		AntennaGain gainsRx = sceneManager->loadGainsFromFileIndBPower(gainPathR.c_str());
		int gainIdRx = sceneManager->registerAntennaGain(gainsRx);
	}

	optix::float3 polarization = make_float3(0.0f, 1.0f, 0.0f);






	std::vector<float3> rx = loadReceiversFromFile(rx_file);
	std::cout<<"Loaded "<<rx.size()<<" receivers from "<<rx_file<<std::endl;
	//std::vector<float3> rx = loadReceiversFromFile(rxfiles[0]);
	//std::vector<float3> rx;
	std::default_random_engine gen;
	std::uniform_real_distribution<float> df(-0.5f, 0.5f);



	int i = 0;
	int maxRx = 1500;
	//for (auto p : rx) {
	for (int p=0; p< rx.size();p +=1 ) {
		//if (i == maxRx) {
		//	break;
		//}
		for (int j = 0; j < mimorx; j++) {
			//for (int j=0; j<replications; j++) {
			float3 posrx = make_float3(rx[p].x, rx[p].y, rx[p].z);
			if (addRandom) {
				posrx.x = posrx.x + df(gen);
				//posrx.y=posrx.y+df(gen);
				posrx.z = posrx.z + df(gen);

			}
			posrx = make_float3(rx[p].x + (-3.5 + j) * delta_rx, rx[p].y, rx[p].z);
			sceneManager->addReceiver(i, posrx, polarization, sphereRadius, sceneManager->printPower);


			if (useGain) {
				sceneManager->registerReceiverGain(i, gainIdRx);
			}
			++i;
		}
	}


		if (useRDN) {
			sceneManager->finishSceneContext();
			int rayD = 10000;
			sceneManager->setRayRange(0.0, 180.0, 0.0, 360.0, rayD, rayD);
			RayDensityNormalizationSimulation* s = dynamic_cast<RayDensityNormalizationSimulation*>(sceneManager->getSimulation(0));
			s->setInitialDensity(((float)sceneManager->getRaySphere().rayCount) / (4 * M_PIf));
			s->setFiltering(2u);

		}
		else {
			sceneManager->finishSceneContext();
			sceneManager->createRaySphere2D(0.0f, 0.1, 180.0f, 0.0f, 0.1, 360.0f);
		}

		//for (int f = 0; f < frequencies.size(); f++) {



		//std::vector<float3> rx = loadReceiversFromFile(rx_file);

		optix::float3 postx = baseStation;
		timerRT.start();
		int indexTx = i + 1;
		int launchNumber=0;
		//ResultReport report;
		mimotx=1;
		//for (int f = 0; f < frequencies.size(); f++) {
			//for (int f = 0; f <1; f++) {
			for (int x = 0; x < mimotx; ++x) {
				//float3 tx_aux = postx + make_float3((-3.5 + x) * delta_tx_x, 0.0, 0.0);
				//float aux_x = postx.x + (-3.5 + x) * delta_tx_x;
				float aux_x = (-3.5 + x) * delta_tx_x;
				for (int y = 0; y < mimotx; ++y) {
                                        //float aux_y=postx.y + (-3.5 + y) * delta_tx_y; 
                                        float aux_y= (-3.5 + y) * delta_tx_y; 
					//float3 tx_aux = make_float3(aux_x, aux_y, postx.z);
					float3 tx_aux = make_float3(aux_x, aux_y,0.0);
					float3 up=make_float3(0.0,1.0,0.0);
					float azimuthDegrees=45.0f; //Rotation of the mimosa panel with respect to x
					const Quaternion ar(up,azimuthDegrees);
					float3 rotatedTx=ar*tx_aux;
					float3 translatedTx=postx+rotatedTx;
					//std::cout<<"tx_aux="<<tx_aux<<"rotatedTx="<<rotatedTx<<"translatedTx="<<translatedTx<<std::endl;
					if (useGain) {
						sceneManager->registerTransmitterGain(indexTx, gainIdTx);
					}
					/*int k = 0;
					  while (k<maxRx) {
					  float3 prx = rx[k];
					  sceneManager->updateReceiver(k, prx, sphereRadius);
					  ++k;
					  }*/
					std::cout << launchNumber<<";f="<<frequencies[frequencies.size()-1]<<";tx="<<tx_aux << std::endl;

					Timer timerLaunch;
					timerLaunch.start();
					ResultReport* rp = sceneManager->transmit(indexTx, txPower, translatedTx, polarization, frequencies[frequencies.size()-1], false, false);
					report.merge(*rp);
					delete rp;
					timerLaunch.stop();
					std::cout << "Time launch= " << timerLaunch.getTime() << std::endl;
					++launchNumber;
					indexTx++;


					/*k = 0;
					  while (k < (rx.size()-maxRx)) {
					  float3 prx = rx[k+maxRx];
					  sceneManager->updateReceiver(k, prx, sphereRadius);
					  ++k;
					  }
					  std::cout << "Launching remaining =" << (rx.size() - maxRx) << "receivers" << std::endl;
					  ResultReport* rp2 = sceneManager->transmit(indexTx, txPower, postx, polarization, frequencies[f], false, false);
					  report.merge(*rp2);
					  delete rp2;
					  */
				}
			}

			std::string rpath = outputFile + "-"  + std::to_string(frequencies[frequencies.size()-1]) + ".csv";
			std::cout << "path=" << rpath << std::endl;
			report.toCSV(rpath);
			report.clear();
		//}
		//report.toCSV(outputFile);
		report.clear();
		timer.stop();
		timerRT.stop();

		std::cout << "Time\t" << timer.getTime() << "\t" << timerRT.getTime() << std::endl;
		delete loader;

		}

void  Lille::loadScenario(bool addRandom, bool useRDN, std::string path, std::string outputFile) {
	Timer timer;
	Timer timerRT;
	ResultReport report;
	std::cout << "Load Lille with scenario " << path << "; freq=" << freq << ";tx=" << baseStation << std::endl;
	timer.start();
	std::vector<float> frequencies;
	//From 5.85 GHz to 5.93, 818 point
	for (int i = 0; i < 818; i++) {
		float f = 5.85e9 + i * 9.78e4;
		frequencies.push_back(f);

	}
	//Init context before doing anything else
	//Enable multichannel since each transmitter use different frequencies
	sceneManager->enableMultiChannel();
	sceneManager->setMinEpsilon(1e-3);
	sceneManager->setUseAntennaGain(useGain);
	sceneManager->enableGenerateRaysOnLaunch();
	//Use field or induced voltage
	//ComputeMode mode = ComputeMode::FIELD;
	ComputeMode mode = ComputeMode::VOLTAGE;
	if (useRDN) {
		RayDensityNormalizationSimulation* sim = new RayDensityNormalizationSimulation(sceneManager);
		sceneManager->setSimulation(sim);
		sim->setComputeMode(mode);
	}
	else {
		if (useDepolarization) {
			LPFlatMeshReflectionSimulation* sim = new LPFlatMeshReflectionSimulation(sceneManager);
			sceneManager->setSimulation(sim);
			sim->setComputeMode(mode);
			//Create log trace (only for one tx and one rx)
			//sim->setEnableTraceLog(false);
			//sim->setEnableSimulation(false);
			//sim->setPrintHits(true);
		}
		else {
			BasicFlatMeshReflectionSimulation* sim = new BasicFlatMeshReflectionSimulation(sceneManager);
			sceneManager->setSimulation(sim);
			sim->setComputeMode(mode);
			//sim->setPrintHits(true);
		}
	}
	//Add diffraction
	SingleDiffraction* simd = new SingleDiffraction(sceneManager);
	sceneManager->setSimulation(simd);
	simd->setComputeMode(mode);
	//Create log trace (only for one tx and one rx)
	//simd->setEnableTraceLog(false);
	//simd->setPrintHits(true);
	//Disable or enable diffraction
	//simd->setEnableSimulation(false);

	sceneManager->initContext(freq);
	//Exceptions
	//sceneManager->enableExceptions();	

	//Load files here
	ScenarioLoader* loader = new ScenarioLoader(sceneManager);
	loader->loadJSONScenario(path);
	int gainIdTx = -1;
	int gainIdRx = -1;
	if (useGain) {
		//Change antenna paths	
		std::string gp("lille/gainT");
		std::string gainPathT = gp + gsufix + ".txt";
		AntennaGain gains = sceneManager->loadGainsFromFileIndBPower(gainPathT.c_str());
		int gainIdTx = sceneManager->registerAntennaGain(gains);
		//Receiver gain
		std::string gpR("lille/gainR");
		std::string gainPathR = gpR + gsufix + ".txt";
		AntennaGain gainsRx = sceneManager->loadGainsFromFileIndBPower(gainPathR.c_str());
		int gainIdRx = sceneManager->registerAntennaGain(gainsRx);
	}

	optix::float3 polarization = make_float3(0.0f, 1.0f, 0.0f);






	std::vector<float3> rx = loadReceiversFromFile(rx_file);
	//std::vector<float3> rx = loadReceiversFromFile(rxfiles[0]);
	//std::vector<float3> rx;
	std::default_random_engine gen;
	std::uniform_real_distribution<float> df(-0.5f, 0.5f);



	int i = 0;
	int maxRx = 1500;
	for (auto p : rx) {
		//if (i == maxRx) {
		//	break;
		//}
		for (int j = 0; j < mimorx; j++) {
			//for (int j=0; j<replications; j++) {
			float3 posrx = make_float3(p.x, p.y, p.z);
			if (addRandom) {
				posrx.x = posrx.x + df(gen);
				//posrx.y=posrx.y+df(gen);
				posrx.z = posrx.z + df(gen);

			}
			posrx = make_float3(p.x + (-3.5 + j) * delta_rx, p.y, p.z);
			sceneManager->addReceiver(i, posrx, polarization, sphereRadius, sceneManager->printPower);


			if (useGain) {
				sceneManager->registerReceiverGain(i, gainIdRx);
			}
			++i;
		}
		}


		if (useRDN) {
			sceneManager->finishSceneContext();
			int rayD = 10000;
			sceneManager->setRayRange(0.0, 180.0, 0.0, 360.0, rayD, rayD);
			RayDensityNormalizationSimulation* s = dynamic_cast<RayDensityNormalizationSimulation*>(sceneManager->getSimulation(0));
			s->setInitialDensity(((float)sceneManager->getRaySphere().rayCount) / (4 * M_PIf));
			s->setFiltering(2u);

		}
		else {
			sceneManager->finishSceneContext();
			sceneManager->createRaySphere2D(0.0f, 0.1, 180.0f, 0.0f, 0.1, 360.0f);
		}

		//for (int f = 0; f < frequencies.size(); f++) {



		//std::vector<float3> rx = loadReceiversFromFile(rx_file);

		optix::float3 postx = baseStation;
		timerRT.start();
		int indexTx = i + 1;
		int launchNumber=0;
		//ResultReport report;
		for (int f = 0; f < frequencies.size(); f++) {
			//for (int f = 0; f <1; f++) {
			for (int x = 0; x < mimotx; ++x) {
				//float3 tx_aux = postx + make_float3((-3.5 + x) * delta_tx_x, 0.0, 0.0);
				//float aux_x = postx.x + (-3.5 + x) * delta_tx_x;
				float aux_x = (-3.5 + x) * delta_tx_x;
				for (int y = 0; y < mimotx; ++y) {
                                        //float aux_y=postx.y + (-3.5 + y) * delta_tx_y; 
                                        float aux_y= (-3.5 + y) * delta_tx_y; 
					//float3 tx_aux = make_float3(aux_x, aux_y, postx.z);
					float3 tx_aux = make_float3(aux_x, aux_y,0.0);
					float3 up=make_float3(0.0,1.0,0.0);
					float azimuthDegrees=45.0f; //Rotation of the mimosa panel with respect to x
					const Quaternion ar(up,azimuthDegrees);
					float3 rotatedTx=ar*tx_aux;
					float3 translatedTx=postx+rotatedTx;
					//std::cout<<"tx_aux="<<tx_aux<<"rotatedTx="<<rotatedTx<<"translatedTx="<<translatedTx<<std::endl;
					if (useGain) {
						sceneManager->registerTransmitterGain(indexTx, gainIdTx);
					}
					/*int k = 0;
					  while (k<maxRx) {
					  float3 prx = rx[k];
					  sceneManager->updateReceiver(k, prx, sphereRadius);
					  ++k;
					  }*/
					std::cout << launchNumber<<";f="<<frequencies[f]<<";tx="<<tx_aux << std::endl;

					Timer timerLaunch;
					timerLaunch.start();
					ResultReport* rp = sceneManager->transmit(indexTx, txPower, translatedTx, polarization, frequencies[f], false, false);
					report.merge(*rp);
					delete rp;
					timerLaunch.stop();
					std::cout << "Time launch= " << timerLaunch.getTime() << std::endl;
					++launchNumber;
					indexTx++;
					/*k = 0;
					  while (k < (rx.size()-maxRx)) {
					  float3 prx = rx[k+maxRx];
					  sceneManager->updateReceiver(k, prx, sphereRadius);
					  ++k;
					  }
					  std::cout << "Launching remaining =" << (rx.size() - maxRx) << "receivers" << std::endl;
					  ResultReport* rp2 = sceneManager->transmit(indexTx, txPower, postx, polarization, frequencies[f], false, false);
					  report.merge(*rp2);
					  delete rp2;
					  */
				}
			}

			std::string rpath = outputFile + "-"  + std::to_string(f) + ".csv";
			std::cout << "path=" << rpath << std::endl;
			report.toCSV(rpath);
			report.clear();
		}
		//report.toCSV(outputFile);
		report.clear();
		timer.stop();
		timerRT.stop();

		std::cout << "Time\t" << timer.getTime() << "\t" << timerRT.getTime() << std::endl;
		delete loader;

		}


		std::vector<Basestation> Lille::loadTransmittersFromFile(std::string file) {
			ScenarioLoader* sl=new ScenarioLoader(sceneManager);
			std::ifstream infile(file);
			if (!infile.good()) {
				infile.close();
				std::cout<<"Error opening "<<file<<std::endl;
				throw  opal::Exception("Lille::loadTransmittersFromFile(): error opening file");
			}
			std::cout<<"Loading transmitters from  "<<file<<std::endl;
			std::string line;
			std::vector<Basestation> tx;
			int id=0;
			while (std::getline(infile, line) ){
				std::cout<<id<<"\t"<<line<<std::endl;
				Basestation b;
				b.id=id;
				++id;
				std::string delimiters("\t");
				std::istringstream iline;
				std::string val;

				iline.str(line);
				getline(iline,val,'\t');
				b.antenna_id=std::stoul(val);
				getline(iline,val,'\t');
				b.postx.x=std::stof(val);
				getline(iline,val,'\t');
				b.postx.y=std::stof(val);
				getline(iline,val,'\t');
				b.postx.z =std::stof(val);
				getline(iline,val,'\t');
				b.azimuth=std::stof(val);	
				getline(iline,val,'\t');
				b.freq=std::stof(val);
				//Values are in MHz
				b.freq=b.freq*1e6;
				
				getline(iline,val,'\t');
				b.site=std::stoul(val);
				getline(iline,val,'\t');
				b.tech=val;
				getline(iline,val,'\t');
				b.lat=std::stof(val);	
				getline(iline,val,'\t');
				b.lon=std::stof(val);	
				b.antenna_index=0;
				tx.push_back(b);
				std::string g="5G";
				const char* found=std::strstr(b.tech.c_str(),g.c_str());
				if (found) {
					for (int i=1; i<7; ++i) {
						Basestation bs=b;
						++id;
						bs.id=id;
						bs.antenna_index=i;
						tx.push_back(bs);
					}
				}
			}
			infile.close();
			delete sl;
			return tx;
		}
		std::vector<float3> Lille::loadReceiversFromFile(std::string file) {
			ScenarioLoader* sl=new ScenarioLoader(sceneManager);
			std::ifstream infile(file);
			if (!infile.good()) {
				infile.close();
				std::cout<<"Error opening "<<file<<std::endl;
				throw  opal::Exception("loadReceiversFromFile(): error opening file");
			}
			//std::cout<<"Loading receivers from  "<<file<<std::endl;
			std::string line;
			std::vector<float3> rx;
			int i=0;
			while (std::getline(infile, line) ){
				optix::float3 v=sl->readFloat3(line);
				//std::cout<<i<<";"<<v<<std::endl;
				rx.push_back(v);
				++i;
			}
			infile.close();
			delete sl;
			return rx;
		}
		void Lille::runTests(std::string test, std::string outputFile) {
			std::vector<float> frequencies;
			//From 5.85 GHz to 5.93, 818 point
			for (int i=0; i<818; i++) {
				float f=5.85e9+i*9.78e4;
				frequencies.push_back(f);

			}

	std::vector <std::string> rxfiles;
	rxfiles.push_back("p3/frame1_last.txt");
	rxfiles.push_back("p3/frame2_last.txt");
	rxfiles.push_back("p3/frame3_last.txt");
	rxfiles.push_back("p3/frame4_last.txt");
	rxfiles.push_back("p3/frame5_last.txt");
	rxfiles.push_back("p3/frame6_last.txt");
	//rxfiles.push_back("p3/fram_laste2p.txt");
			std::vector<int> tokens=parseTestString(test);

			this->replications=1;

			freq=frequencies[tokens[0]];
			gsufix=std::to_string(tokens[0]);
			std::string path;
			path="p3/p3_new.json";
			//rx_file="p3/frame1.txt";
			rx_file=rxfiles[tokens[3]];
			txPower=.001f; //1 W


			bool useRDN=false;
			if (tokens[1]==1) {
				useRDN=true;
			}
			useGain=false;
			if (tokens[2]==1) {
				useGain=true;
			}
			//true here adds a random displacement to receivers positions
			bool addRandom=false;
			//loadScenario(addRandom,useRDN, path, outputFile);
			angleOfDeparture(addRandom,useRDN, path, outputFile);
		}
		void Lille::test() {
			this->replications=1;

			freq=5.9e9;
			std::string path("p3_new.json");;
			std::string outputFile("test.csv");;
			bool useRDN=false;
			useGain=false;
			loadScenario(false,useRDN, path, outputFile);
		}

