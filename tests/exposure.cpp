/***************************************************************/
//
//Copyright (c) 2021 Esteban Egea-Lopez http://girtel.upct.es/~simulations
//
/**************************************************************/
#include "exposure.h"
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
#include <algorithm>
#include <random>
#include <numeric>

using namespace opal;
using namespace optix;

Exposure::Exposure(OpalSceneManager*   sceneManager, float sphereRadius, bool useDepolarization) : Lille(sceneManager, sphereRadius, useDepolarization) {
	this->peakFraction=0.6f;
	this->lowFraction=0.05f;
	this->maxUsersPerCell=90;
	this->fraction=peakFraction;
	this->maxFreq=90;
}

void Exposure::configureSimulations(bool useRDN ) {
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
	//For a large scenario we need to adjust to the available memory in the GPU
// 	simd->setAdjustToMemorySize(true);
//	simd->setMaxMemoryFraction(0.5f);	
	
	//Create log trace (only for one tx and one rx)
	//simd->setEnableTraceLog(false);
	//simd->setPrintHits(true);
	//Disable or enable diffraction
	//simd->setEnableSimulation(false);

	sceneManager->initContext(freq);
	//Exceptions
	//sceneManager->enableExceptions();	

}
void Exposure::setLaunchSizeAndFinish(bool useRDN ) {
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
}
std::vector<float3> Exposure::sampleActiveReceivers(std::vector<float3>& grid, Basestation base,float distance, float fraction) {
	std::vector<float3> active;
	auto points=getPointsAtDistance(grid,distance, base.postx);
	std::cout<<base.id<<"; points at "<<distance<<" m ="<<points.size()<<std::endl;
	int au=floor(points.size()*fraction);
	if (au>maxUsersPerCell) {
		au=maxUsersPerCell;
	}
	//Generate a number of slots equal to total capacity
	std::vector<int> genSlots(points.size());
	std::iota(begin(genSlots), end(genSlots), 0);
	//Generate a random permutation of all the points, to avoid duplicate numbers
	std::mt19937 gen(4557);
	std::shuffle(begin(genSlots), end(genSlots), gen);

	//Take only the first au users
	for (int i=0; i<au; i++) {
		//std::cout<<i<<"gen="<<genSlots[i]<<";"<<points[genSlots[i]]<<std::endl;
		active.push_back(points[genSlots[i]]);

	}
	return active;

}
std::vector<int> Exposure::assignActiveReceiversToBeam(std::vector<float3>& au, Basestation base)  {
	//Group all the active users  in corresponding beams and return the number of beam in a vector.
	//Note that not all beams may be used, basestation will only transmit for the beams used
	std::vector<int> beams(7);
	for (int i=0; i< beams.size(); ++i) {
		beams[i]=0;
	}


	//To get the active users in a sector we change to a local basestation frame
	//
	float deg2rad=M_PIf/180.f;
	float azimRad=base.azimuth*deg2rad;
	//z axis go along azimuth orientation
	float3 zo=normalize(make_float3(sin(azimRad),0.0,cos(azimRad)));
	float xdeg=(base.azimuth+90)*deg2rad;
	float3 xo=normalize(make_float3(sin(xdeg),0.0,cos(xdeg)));
	//std::cout<<"dot="<<dot(zo,xo)<<std::endl;
	float3 yo=make_float3(0.0,1.0,0.0);
	Matrix4x4 t=GeometryUtils::transformToBasis(xo,yo,zo, base.postx);



	//Compute azimuths of all points with respecto to the base station orientation
	for (auto p: au) {

		float3 tp =GeometryUtils::transformPoint(p,t);
		//Now, after transforming the points to the basestation local frame, only points with positive z are
		//considered to be in the sector
		//std::cout<<p<<"tp="<<tp<<"t="<<t<<std::endl;
		if (tp.z>=0) {
			float azt=atan2(tp.x,tp.z);
			float aztd=azt/deg2rad;
			float delta=17.0/2;
			//std::cout<<p<<"; tp="<<tp<<"t="<<t<<"az="<<(azt/deg2rad)<<std::endl;
			if (aztd>-delta && aztd<=delta) {
				beams[0] +=1;
			} else if ( aztd>delta && aztd <=delta+17) {
				beams[2] +=1;
			} else if ( aztd>delta+17 && aztd <=34+delta) {
				beams[4] +=1;
			} else if ( aztd>34+delta && aztd <=90) {
				beams[6] +=1;
			} else if ( aztd>-17-delta && aztd <=-delta) {
				beams[1] +=1;
			} else if ( aztd> (-34-delta) && aztd <(-17-delta)) {
				beams[3] +=1;
			} else if ( aztd>=-90 && aztd <(-34-delta)) {
				beams[5] +=1;
			} else {
				std::cout<<"ERROR; tp="<<tp<<"t="<<t<<"az="<<(azt/deg2rad)<<std::endl;
				throw  opal::Exception("Exposure::assignActiveReceiversToBeam");
			}
		}

		//		azimuths.push_back(a);
	}
	int users=0;
	for (auto b: beams) {
		users +=b;
		//std::cout<<b<<std::endl;
	}
		std::cout<<"Total users in sector="<<users<<std::endl;

	//	std::sort(azimuths.begin(), azimuths.end());
	return beams;
}
std::vector<float3> Exposure::getPointsAtDistance(std::vector<float3>& grid,float distance, float3 tx) {
	std::vector<float3> points;
	for (auto p: grid) {
		float d= length(tx-p);
		if (d<=distance) {
			points.push_back(p);
		}
	}	
	return points;

}
float Exposure::getMeanBasestationDistace(std::vector<Basestation>& bases) {
	float d=0;
	float n=0;
	for (int i=0; i< bases.size(); ++i) {
		float aux=0;
		for (int j=0; j< bases.size(); ++j) {
			if (j!=i) {
				aux += length(bases[i].postx-bases[j].postx);
			}
		}
		d += (aux/(bases.size()-1));
		++n;
	}
	return (d/n);
}

void  Exposure::executeExposure(bool addRandom, bool useRDN, std::string scenarioPath, std::string outputFile, std::string txPath) {
	Timer timer;
	Timer timerRT;
	std::cout << "Load exposure with scenario " << scenarioPath << "; txPath=" << txPath <<  std::endl;
	timer.start();
	configureSimulations(addRandom);

	//Load files here
	std::cout<<"Loading basestations from "<<txPath<<std::endl;
	std::vector<Basestation> tx=loadBasestationsFromFile(txPath,false);
	for (auto bs : tx) {
		std::cout<<bs.toString();
	}
	int num_tx=tx.size();
	float avd=getMeanBasestationDistace(tx);
	float range=avd/2.0;
	std::cout<<"Loaded a total of "<<num_tx<<" basestations"<<"with average distance "<<avd<<std::endl;
	//Register antenna patterns
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



	std::cout<<"Loading  grid points from "<<rx_file<<std::endl;
	std::vector<float3> grid = loadReceiversFromFile(rx_file);
	std::cout<<"Loaded a total of "<<grid.size()<<" grid points"<<std::endl;
	std::cout<<"Loading  sensors points from "<<sensors_file<<std::endl;
	std::vector<float3> sensors = loadReceiversFromFile(sensors_file);
	std::cout<<"Loaded a total of "<<sensors.size()<<" sensors"<<std::endl;

		
	//Basestation bs;
	//bs.azimuth=45;
	//bs.postx=make_float3(0.0,0.0,0.0);
	//std::vector<float3> au;
	//au.push_back(make_float3(-1.0,0.0,1.0));
	//au.push_back(make_float3(cos(M_PIf/4),1.0,sin(M_PIf/4)));
	//au.push_back(make_float3(-10,0.0,1.0));
	////au.push_back(make_float3(cos(M_PIf/4),1.0,sin(M_PIf/4)));
	////au.push_back(make_float3(-sin(M_PIf/3),1.0,cos(M_PIf/3)));
	Basestation aux=tx[tx.size()-1];
	int i = aux.id+1;
	//int maxReceivers = 6790;
	int maxReceivers = 895;
	int totalReceiverPoints=grid.size()+sensors.size();
	if (maxReceivers >=totalReceiverPoints) {
		for (auto p : grid) {

			float3 posrx = make_float3(p.x, p.y, p.z);

			sceneManager->addReceiver(i, posrx, polarization, sphereRadius, sceneManager->printPower);
			//No gain for receivers. Uncomment if necessary
			//if (useGain) {
			//	sceneManager->registerReceiverGain(i, gainIdRx);
			//}
			++i;
		}

		for (auto p : sensors) {

			float3 posrx = make_float3(p.x, p.y, p.z);

			sceneManager->addReceiver(i, posrx, polarization, sphereRadius, sceneManager->printPower);
			++i;
		}
	} else {
		//Merge into grid
		std::cout<<"grid size="<<grid.size()<<std::endl;
		grid.insert(grid.end(), sensors.begin(), sensors.end());
		std::cout<<"grid size="<<grid.size()<<std::endl;
		//Just add dummy receivers, later we update them
		for (int k=0; k<maxReceivers;++k) {
			float3 posrx = grid[k];

			sceneManager->addReceiver(i, posrx, polarization, sphereRadius, sceneManager->printPower);
			++i;
		}
	}

	//Load scenario
	ScenarioLoader* loader = new ScenarioLoader(sceneManager);
	loader->loadJSONScenario(scenarioPath);

	setLaunchSizeAndFinish(useRDN);
	

	//Loop over all basestations
	for (auto bs : tx) {
		timerRT.start();
		optix::float3 postx = bs.postx;
		//sample active users for this basestation (sector)
		std::vector<float3> active_rx=sampleActiveReceivers(grid,bs, range, fraction);
		//Assign receiver to beam 
		std::vector<int> users_in_beam= assignActiveReceiversToBeam(active_rx,bs);
		//Loop over all the beams that have active users
		for (int  beam=0; beam<users_in_beam.size(); ++beam) {
			timerRT.start();
			if (users_in_beam[beam]==0) {
				continue;
			}
			if (useGain) {
				sceneManager->registerTransmitterGain(bs.id, gainIdTx[beam]);
			}
			//Loop over all receivers in case we have to split them due to memory constraints
			if (maxReceivers <=totalReceiverPoints) {
				int numberOfBatches=floor(grid.size()/maxReceivers);
				std::cout<<"Number of batches="<<numberOfBatches<<std::endl;
				ResultReport report;
				for (int indexBatch=0; indexBatch<numberOfBatches; ++indexBatch ) {
					int receiverIndex=aux.id+1;
					while (receiverIndex<maxReceivers) {
						float3 posrx = grid[indexBatch*maxReceivers + receiverIndex];
						std::cout<<"update pos="<<posrx<<"index="<<(indexBatch*maxReceivers + receiverIndex)<<std::endl;
						sceneManager->updateReceiver(receiverIndex, posrx); 
						++receiverIndex;
					}

					//Transmitter position
					optix::float3 postx = bs.postx;

					Timer launchTimer;
					launchTimer.start();
					//Orientate antenna pattern
					//No tilt, we used antenna index up to 6, only 7 had a tilt of -5 (upward)
					Matrix<4,4> t=sceneManager->orientateAntennaPattern(bs.azimuth,0.0);
					sceneManager->getActiveTransmitters()[0]->transformToPolarization=t;
					ResultReport* rp = sceneManager->transmit(bs.id, txPower, postx, polarization, bs.freq, false, true);
					//Multiple transmitters. Merge report
					report.merge(*rp);
					launchTimer.stop();
					std::cout <<"Launch time="<< launchTimer.getTime() << std::endl;
					//delete rp; //TODO: if we delete we cannot save later the report, but it is leaked now, should collect and delete after writing
					//Start for new receivers
				}
				//Save one report per basestation and index
				std::string rpath=outputFile+"-"+std::to_string(bs.id)+"-"+std::to_string(beam)+".csv";
				report.toCSV(rpath);

			} else {


				////Orientate antenna pattern
				////No tilt, we used antenna index up to 6, only 7 had a tilt of -5 (upward)
				////Diagramas are already shifted  from the zero azimuth, so we only need to orientate them to the basestation azimuth
				Matrix<4,4> t=sceneManager->orientateAntennaPattern(bs.azimuth,0.0);
				sceneManager->getActiveTransmitters()[0]->transformToPolarization=t;
				ResultReport* rp = sceneManager->transmit(bs.id, txPower, postx, polarization, bs.freq, false, true);

				//Save one report per basestation and index
				std::string rpath=outputFile+"-"+std::to_string(bs.id)+"-"+std::to_string(beam)+".csv";
				rp->toCSV(rpath);

				//report.merge(*rp);
				delete rp;
			}
			timerRT.stop();
			std::cout << "Basestation time=" << timerRT.getTime() << std::endl;
		}
		//Save additional data per basestation
		std::string dpath=outputFile+"-"+std::to_string(bs.id)+"-data.csv";
		std::string header="active_users,";
		std::string data=std::to_string(active_rx.size())+",";

		for (int  beam=0; beam<users_in_beam.size()-1; ++beam) {
			header += "beam_"+std::to_string(beam)+",";
			data += std::to_string(users_in_beam[beam])+",";

		}
		header += "beam_"+std::to_string(users_in_beam.size()-1);
		data += std::to_string(users_in_beam[users_in_beam.size()-1]);
		std::cout<<header<<std::endl;
		std::cout<<data<<std::endl;
		saveData(dpath,header,data);
		//Save one report per basestation
		//std::string rpath=outputFile+"-"+std::to_string(bs.id)+".csv";
		//report.toCSV(rpath);
		//report.clear();
		//Save additional data
	}


	//Multiple transmitters. save full report csv
//	std::string rpath=outputFile+"-last.csv";
//	report.toCSV(rpath);
	timer.stop();

	std::cout << "Total time\t" << timer.getTime() << "\t" << timerRT.getTime() << std::endl;
	delete loader;
}
void Exposure::saveData(std::string outputFile, std::string header, std::string data) {
	std::ofstream myfile;
	myfile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
	try {
		myfile.open(outputFile);
		myfile<<header<<std::endl;
		myfile<<data<<std::endl;
		myfile.close();
	} catch (std::ifstream::failure e) {
		std::cerr << "Exception opening/reading/closing data file\n"; 

	}
}
void Exposure::run5GBeamSelection(std::string test, std::string outputFile) {
	std::vector<int> tokens = parseTestString(test);

	this->replications = 1;

	std::string path;
	path = "emf_5G/waz_scene2.json";
	this->rx_file = "emf_5G/rx-waz.txt";
	this->sensors_file = "emf_5G/waz_50sens.txt";
	//std::string txPath="emf_5G/wazemmes_tx_info_unity.csv";
	std::string txPath="emf_5G/tx_waz_5g.txt";
	if (tokens[0]==1) {
		txPath="emf_5G/euratech_tx_info_unity.csv";
	}
	this->fraction=peakFraction;
	if (tokens[1]==1) {
		this->fraction=lowFraction;
	}
	txPower = 1.0f; //1 W
	this->useGain = true;


	bool useRDN = false;
	bool addRandom=false;
	executeExposure(addRandom, useRDN, path, outputFile, txPath);
}
void Exposure::run5GOneOperator(std::string test, std::string outputFile) {
	std::vector<int> tokens = parseTestString(test);

	this->replications = 1;

	std::string path;
	path = "emf_5G/waz_scene2.json";
	this->rx_file = "emf_5G/rx-waz.txt";
	this->sensors_file = "emf_5G/waz_50sens.txt";
	//std::string txPath="emf_5G/orange1.txt";
	std::string txPath="emf_5G/orange.txt";
	if (tokens[0]==1) {
		txPath="emf_5G/orange2.txt";
	}
	this->fraction=peakFraction;
	if (tokens[1]==1) {
		this->fraction=lowFraction;
	}
	txPower = 1.0f; //1 W
	this->useGain = true;


	bool useRDN = false;
	bool addRandom=false;
	executeExposureOneOperator(addRandom, useRDN, path, outputFile, txPath);
}
void  Exposure::executeExposureOneOperator(bool addRandom, bool useRDN, std::string scenarioPath, std::string outputFile, std::string txPath) {
	Timer timer;
	Timer timerRT;
	std::cout << "Load exposure one operator with scenario " << scenarioPath << "; txPath=" << txPath <<  std::endl;
	timer.start();
	configureSimulations(addRandom);

	//Load files here
	std::cout<<"Loading basestations from "<<txPath<<std::endl;
	std::vector<Basestation> tx=loadBasestationsFromFile(txPath,false);
	for (auto bs : tx) {
		std::cout<<bs.toString();
	}
	int num_tx=tx.size();
	float avd=getMeanBasestationDistace(tx);
	float range=avd/2.0;
	std::cout<<"Loaded a total of "<<num_tx<<" basestations"<<"with average distance "<<avd<<std::endl;
	//Register antenna patterns
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



	std::cout<<"Loading  grid points from "<<rx_file<<std::endl;
	std::vector<float3> grid = loadReceiversFromFile(rx_file);
	std::cout<<"Loaded a total of "<<grid.size()<<" grid points"<<std::endl;
	std::cout<<"Loading  sensors points from "<<sensors_file<<std::endl;
	std::vector<float3> sensors = loadReceiversFromFile(sensors_file);
	std::cout<<"Loaded a total of "<<sensors.size()<<" sensors"<<std::endl;

		
	//Basestation bs;
	//bs.azimuth=45;
	//bs.postx=make_float3(0.0,0.0,0.0);
	//std::vector<float3> au;
	//au.push_back(make_float3(-1.0,0.0,1.0));
	//au.push_back(make_float3(cos(M_PIf/4),1.0,sin(M_PIf/4)));
	//au.push_back(make_float3(-10,0.0,1.0));
	////au.push_back(make_float3(cos(M_PIf/4),1.0,sin(M_PIf/4)));
	////au.push_back(make_float3(-sin(M_PIf/3),1.0,cos(M_PIf/3)));
	Basestation aux=tx[tx.size()-1];
	int i = aux.id+1;
	int maxReceivers = 6790;
	//int maxReceivers = 895;
	int totalReceiverPoints=grid.size()+sensors.size();
	if (maxReceivers >=totalReceiverPoints) {
		for (auto p : grid) {

			float3 posrx = make_float3(p.x, p.y, p.z);

			sceneManager->addReceiver(i, posrx, polarization, sphereRadius, sceneManager->printPower);
			//No gain for receivers. Uncomment if necessary
			//if (useGain) {
			//	sceneManager->registerReceiverGain(i, gainIdRx);
			//}
			++i;
		}

		for (auto p : sensors) {

			float3 posrx = make_float3(p.x, p.y, p.z);

			sceneManager->addReceiver(i, posrx, polarization, sphereRadius, sceneManager->printPower);
			++i;
		}
	} else {
		//Merge into grid
		std::cout<<"grid size="<<grid.size()<<std::endl;
		grid.insert(grid.end(), sensors.begin(), sensors.end());
		std::cout<<"grid size="<<grid.size()<<std::endl;
		//Just add dummy receivers, later we update them
		for (int k=0; k<maxReceivers;++k) {
			float3 posrx = grid[k];

			sceneManager->addReceiver(i, posrx, polarization, sphereRadius, sceneManager->printPower);
			++i;
		}
	}

	//Load scenario
	ScenarioLoader* loader = new ScenarioLoader(sceneManager);
	loader->loadJSONScenario(scenarioPath);

	setLaunchSizeAndFinish(useRDN);
	

	//Loop over all basestations
	for (auto bs : tx) {
		timerRT.start();
		optix::float3 postx = bs.postx;
		//sample active users for this basestation (sector)
		//std::vector<float3> active_rx=sampleActiveReceivers(grid,bs, range, fraction);
		//Assign receiver to beam 
		//std::vector<int> users_in_beam= assignActiveReceiversToBeam(active_rx,bs);
		//Loop over all beams
		std::cout << "Executing basestation: " << bs.toString() << std::endl;
		for (int  beam=0; beam<7; ++beam) {
			ResultReport report;
			//Loop over all frequencies
			for (int f = 0; f<maxFreq; ++f) {
				float freq=bs.freq  + (f*100e6);
				timerRT.start();
				if (useGain) {
					sceneManager->registerTransmitterGain(bs.id, gainIdTx[beam]);
				}
				//Loop over all receivers in case we have to split them due to memory constraints
				if (maxReceivers <=totalReceiverPoints) {
					int numberOfBatches=floor(grid.size()/maxReceivers);
					std::cout<<"Number of batches="<<numberOfBatches<<std::endl;
					for (int indexBatch=0; indexBatch<numberOfBatches; ++indexBatch ) {
						int receiverIndex=aux.id+1;
						while (receiverIndex<maxReceivers) {
							float3 posrx = grid[indexBatch*maxReceivers + receiverIndex];
							std::cout<<"update pos="<<posrx<<"index="<<(indexBatch*maxReceivers + receiverIndex)<<std::endl;
							sceneManager->updateReceiver(receiverIndex, posrx); 
							++receiverIndex;
						}

						//Transmitter position
						optix::float3 postx = bs.postx;

						Timer launchTimer;
						launchTimer.start();
						//Orientate antenna pattern
						//No tilt, we used antenna index up to 6, only 7 had a tilt of -5 (upward)
						Matrix<4,4> t=sceneManager->orientateAntennaPattern(bs.azimuth,0.0);
						sceneManager->getActiveTransmitters()[0]->transformToPolarization=t;
						ResultReport* rp = sceneManager->transmit(bs.id, txPower, postx, polarization, freq, false, true);
						//Multiple transmitters. Merge report
						report.merge(*rp);
						launchTimer.stop();
						std::cout <<"Launch time="<< launchTimer.getTime() << std::endl;
						//delete rp; //TODO: if we delete we cannot save later the report, but it is leaked now, should collect and delete after writing
						//Start for new receivers
					}

				} else {


					////Orientate antenna pattern
					////No tilt, we used antenna index up to 6, only 7 had a tilt of -5 (upward)
					////Diagramas are already shifted  from the zero azimuth, so we only need to orientate them to the basestation azimuth
					Matrix<4,4> t=sceneManager->orientateAntennaPattern(bs.azimuth,0.0);
					sceneManager->getActiveTransmitters()[0]->transformToPolarization=t;
					Timer launchTimer;
					launchTimer.start();
					ResultReport* rp = sceneManager->transmit(bs.id, txPower, postx, polarization, freq, false, true);

					//Save one report per basestation and index
					//rp->toCSV(rpath);

					report.merge(*rp);
					launchTimer.stop();
					std::cout <<"Launch time="<< launchTimer.getTime() << std::endl;
					//report.merge(*rp);
					delete rp;
				}
			}
			//Save one report per basestation and index
			std::string rpath=outputFile+"-"+std::to_string(bs.id)+"-"+std::to_string(beam)+".csv";
			report.toCSV(rpath);
			timerRT.stop();
			std::cout << "Basestation time=" << timerRT.getTime() << std::endl;
		}
		//Save additional data per basestation
		//std::string dpath=outputFile+"-"+std::to_string(bs.id)+"-data.csv";
		//std::string header="active_users,";
		//std::string data=std::to_string(active_rx.size())+",";

		//for (int  beam=0; beam<users_in_beam.size()-1; ++beam) {
		//	header += "beam_"+std::to_string(beam)+",";
		//	data += std::to_string(users_in_beam[beam])+",";

		//}
		//header += "beam_"+std::to_string(users_in_beam.size()-1);
		//data += std::to_string(users_in_beam[users_in_beam.size()-1]);
		//std::cout<<header<<std::endl;
		//std::cout<<data<<std::endl;
		//saveData(dpath,header,data);
		////Save one report per basestation
		//std::string rpath=outputFile+"-"+std::to_string(bs.id)+".csv";
		//report.toCSV(rpath);
		//report.clear();
		//Save additional data
	}


	//Multiple transmitters. save full report csv
//	std::string rpath=outputFile+"-last.csv";
//	report.toCSV(rpath);
	timer.stop();

	std::cout << "Total time\t" << timer.getTime() << "\t" << timerRT.getTime() << std::endl;
	delete loader;
}

std::vector<Basestation> Exposure::loadBasestationsFromFile(std::string file, bool replicatePerDiagram) {
	ScenarioLoader* sl=new ScenarioLoader(sceneManager);
	std::ifstream infile(file);
	if (!infile.good()) {
		infile.close();
		std::cout<<"Error opening "<<file<<std::endl;
		throw  opal::Exception("Exposure::loadBasestationsFromFile error opening file");
	}
	std::cout<<"Loading basestations from  "<<file<<std::endl;
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
		//std::cout<<";val="<<val<<std::endl;
		b.site=std::stoul(val);

		getline(iline,val,'\t');
		//std::cout<<";val="<<val<<std::endl;
		b.postx.x=std::stof(val);

		getline(iline,val,'\t');
		//std::cout<<";val="<<val<<std::endl;
		b.postx.y=std::stof(val);
		getline(iline,val,'\t');
		//std::cout<<";val="<<val<<std::endl;
		b.postx.z =std::stof(val);
		getline(iline,val,'\t');
		//std::cout<<";val="<<val<<std::endl;
		b.azimuth=std::stof(val);	
		getline(iline,val,'\t');
		//std::cout<<";val="<<val<<std::endl;
		b.freq=std::stof(val);

		//Values are in MHz
		b.freq=b.freq*1e6;

		getline(iline,val,'\t');
		//std::cout<<";val="<<val<<std::endl;
		b.freq_end=std::stof(val);
		b.freq_end=b.freq_end*1e6;
		getline(iline,val,'\t');
		//std::cout<<";val="<<val<<std::endl;
		b.antenna_id=std::stoul(val);

		getline(iline,val,'\t');
		//std::cout<<";val="<<val<<std::endl;
		b.tech=val;

		getline(iline,val,'\t');
		//std::cout<<";val="<<val<<std::endl;
		b.lat=std::stof(val);	

		getline(iline,val,'\t');
		//std::cout<<";val="<<val<<std::endl;
		b.lon=std::stof(val);	

		b.antenna_index=0;
		tx.push_back(b);
		if (replicatePerDiagram) {
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
	}
	infile.close();
	delete sl;
	return tx;
}

