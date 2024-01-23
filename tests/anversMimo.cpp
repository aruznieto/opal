/***************************************************************/
//
//Copyright (c) 2019 Esteban Egea-Lopez http://girtel.upct.es/~simulations
//
/**************************************************************/
#include "anversMimo.h"
#include "tests.h"
#include "../timer.h"
#include <memory>
#include <iostream>
#include <fstream>
#include "../flatSimulation.h"
#include "../curvedMeshSimulation.h"
#include "../curvedFlatMeshSimulation.h"
#include "../rayDensityNormalizationSimulation.h"
#include "../singleDiffraction.h"
#include "../util.h"
#include "../results.h"
using namespace opal;
using namespace optix;
AnversMimo::AnversMimo(OpalSceneManager* sceneManager, float sphereRadius, bool useDepolarization) : AnversTests(sceneManager,sphereRadius,useDepolarization) {
	this->delta_rx=0.075f;
	this->delta_tx=0.025f;
	this->mimotx=8;
	this->mimotx_x=4;
	this->mimotx_y=8;
	this->mimorx=8;
	std::cout <<" Creating AnversMimo" <<std::endl;
}
void AnversMimo::runLille(std::string test, std::string outputPath) {
	this->useAntennaGain=false;
	//this->useRDN=false;
	std::vector<float3> postx(1);
	std::vector<float3> posrx(1);
	std::vector<float3> pol(4);
	std::vector<float> freq;

	//Origin seem to be at NW corner of the tunnel cross section
	//postx[0]=make_float3(0.50f,2.2f-6.0f, 0.0f);
	//postx[0]=make_float3(0.5f+delta_tx*3.5,(6.0 -(1+delta_tx*3.5))-6.0f, 0.0f);
	postx[0]=make_float3(0.5f,-1.0f, 0.0f); //Border is at 1 m from the top wall
	tx=postx[0];
	//posrx[0]=make_float3(3.2f,2.7f-6.0f, 20.0f);
	//posrx[0]=make_float3(2.55,2.7f-6.0f, 50.0f);
	posrx[0]=make_float3(2.55,1.6f-6.0f, 50.0f);
	rx=posrx[0];
	const float3 VR=normalize(make_float3(1.0,1.0,0.0));
	const float3 VL=normalize(make_float3(-1.0,1.0,0.0));
	pol[0]=VR;
	pol[1]=VL;
	pol[2]=V;
	pol[3]=H;
	float f=5.9e9;
	//float f=1.35e9;
	//for (int i=0; i<818;++i) {
	//	freq.push_back(f);
	//	f +=80e6;
	//}
	//float deltaf=6.666666666666031e+06;
	float deltaf=3.225806451612473e+06;
	for (int i=0; i<32;++i) {
		freq.push_back(f);
		f +=deltaf;
	}
	this->delta_rx=0.075f;
	this->delta_tx=0.025f;
	this->mimotx_x=4;
	this->mimotx_y=8;
	this->mimorx=4;
	std::cout<<"Total number of frequencies "<<freq.size()<<std::endl;
	//Parse tests
	std::vector<int> tokens=parseTestString(test);
	frequency=freq[tokens[0]];
	//Polarization
	polarizationTx=pol[tokens[1]];
	polarizationRx=pol[tokens[2]]; //New token...
	bool field=true;
	bool sectorized=false;
	bool forward=false;
	bool useReflection=true;
	bool useDiffraction=false;
	emptyTunnel=false;
	bool multitransmitter=false;
	this->useRDN=true;
	this->increaseRadius=true;
	//std::string filePath=outputPath+"trucks-"+std::to_string(tokens[0])+"-"+std::to_string(tokens[1])+"-"+std::to_string(tokens[2])+".csv";
	if (tokens[3]==0) {
		std::string filePath=outputPath+"trucks-"+std::to_string(tokens[1])+"-"+std::to_string(tokens[2])+".csv";
		runFlatTrucksLille(true,field,useReflection,useDiffraction, sectorized, forward,multitransmitter, filePath);
	} else {
		posrx[0]=make_float3(2.55,1.6f-6.0f, 0.0f);
		rx=posrx[0];
		std::string filePath=outputPath+"cars-"+std::to_string(tokens[1])+"-"+std::to_string(tokens[2])+".csv";
		runFlatCarTrucksLille(true,field,useReflection,useDiffraction, sectorized, forward,multitransmitter, filePath);
	}
}
void AnversMimo::runTests(std::string test, std::string outputPath) {
	this->useAntennaGain=false;
	//this->useRDN=false;
	std::vector<float3> postx(1);
	std::vector<float3> posrx(1);
	std::vector<float3> pol(4);
	std::vector<float> freq;

	//Origin seem to be at NW corner of the tunnel cross section
	postx[0]=make_float3(0.50f,2.2f-6.0f, 0.0f);
	//postx[0]=make_float3(0.5f+delta_tx*3.5,(6.0 -(1+delta_tx*3.5))-6.0f, 0.0f);
	tx=postx[0];
	//posrx[0]=make_float3(3.2f,2.7f-6.0f, 20.0f);
	posrx[0]=make_float3(2.55,2.7f-6.0f, 50.0f);
	rx=posrx[0];
	const float3 VR=normalize(make_float3(1.0,1.0,0.0));
	const float3 VL=normalize(make_float3(-1.0,1.0,0.0));
	pol[0]=VR;
	pol[1]=VL;
	pol[2]=V;
	pol[3]=H;
	//float f=5.85e9;
	float f=5.9e9;
	//for (int i=0; i<818;++i) {
	//	freq.push_back(f);
	//	f +=80e6;
	//}
	//float deltaf=6.666666666666031e+06;
	float deltaf=3.225806451612473e+06;
	for (int i=0; i<32;++i) {
		freq.push_back(f);
		f +=deltaf;
	}
	this->delta_rx=0.075f;
	this->delta_tx=0.025f;
	this->mimotx=8;
	this->mimorx=8;
	std::cout<<"Total number of frequencies "<<freq.size()<<std::endl;
	//Parse tests
	std::vector<int> tokens=parseTestString(test);
	frequency=freq[tokens[0]];
	//Polarization
	polarizationTx=pol[tokens[1]];
	polarizationRx=pol[tokens[2]]; //New token...
	bool field=true;
	bool sectorized=false;
	bool forward=false;
	bool useReflection=true;
	bool useDiffraction=false;
	emptyTunnel=false;
	bool multitransmitter=false;
	std::string filePath=outputPath+"trucks-"+std::to_string(tokens[0])+"-"+std::to_string(tokens[1])+"-"+std::to_string(tokens[2])+".csv";
	runFlatTrucks(true,field,useReflection,useDiffraction, sectorized, forward,multitransmitter, filePath);
}
void AnversMimo::runFlatEmpty(bool half, bool computeField, bool useReflection, bool useDiffraction, bool sectorized, bool forward, bool multitransmitter) {
	OpalSimulation* sim;
	if (useRDN) {
		sim = new RayDensityNormalizationSimulation(sceneManager);
		//RayDensityNormalizationSimulation* sim = new RayDensityNormalizationSimulation(sceneManager);
	} else {	
		sim= new LPFlatMeshReflectionSimulation(sceneManager);
		//LPFlatMeshReflectionSimulation* sim= new LPFlatMeshReflectionSimulation(sceneManager);
	}

	//sim->setExecutionMethod(RDNExecutionMode::HITINFO);
	//sim->setEnableTraceLog(true);
	//sim->setPrintHits(true);
	sceneManager->setSimulation(sim);
	ComputeMode mode=ComputeMode::VOLTAGE;
	if (computeField) {
		std::cout<<"Computing FIELD" <<std::endl;
		mode=ComputeMode::FIELD;
	}
	sim->setComputeMode(mode);
	sim->setEnableSimulation(useReflection);
	sceneManager->enableGenerateRaysOnLaunch();	
	//Add diffraction simulation
	SingleDiffraction* simd= new SingleDiffraction(sceneManager);
	sceneManager->setSimulation(simd);
	simd->setComputeMode(mode);
	//simd->setEnableTraceLog(true);

	if (this->useAntennaGain) {
		sceneManager->setUseAntennaGain(true);
	}
	sceneManager->initContext(frequency);
	simd->setEnableSimulation(useDiffraction);
	if (emptyTunnel) {
		simd->setEnableSimulation(false);
	}
	//simd->setPrintHits(true);
	//sceneManager->getSimulation(0)->setPrintHits(true);


	Timer timer;
	std::cout<<"\tf="<<(frequency/1e6)<< " MHz; tx=" <<tx<<"; polarization="<<polarizationTx<<std::endl;
	std::cout<<"\trx="<<rx<<"; radius="<<sphereRadius;
	if (increaseRadius) {
		std::cout<<"+ "<<sphereDelta;
	}
	std::cout<<"; polarization="<<polarizationRx<<std::endl;
	sceneManager->setMinEpsilon(1e-4f);
	//std::string path("meshes/anvers/trucks");
	//loadAnversScenario(path);



	//Flat with trucks
	//std::string path("straight-trucks-metal-one.json");
	//std::string path("straight-trucks-metal-4.json");
	//loadStraightAnversTunnel();	
	if (!emptyTunnel) {
		std::string path("straight-trucks-metal.json");
		//std::string path("straight-single-truck.json");
		loadAnversJsonScenario(path);
	} else {

		float width=10.2f;
		float height=6.0f;
		float length=1070.0f;
		Matrix4x4 tm;
		tm.setRow(0, make_float4(width, 0, 0, width/2));
		tm.setRow(1, make_float4(0, height, 0, -height/2));
		tm.setRow(2, make_float4(0, 0, length, length/2.0f));
		tm.setRow(3, make_float4(0, 0, 0, 1));
		MaterialEMProperties emProp1;
		emProp1.dielectricConstant = make_float2(5.0f, -60.0f*sceneManager->getChannelParameters().waveLength*0.01f);
		emProp1.tattenuation = make_float2(0.1f,-75.f );
		loadTransformedSquareTunnel(emProp1, tm);
	}
	int gainId;	
	if (this->useAntennaGain && !emptyTunnel) {
		if (forward) {
			AntennaGain gains=sceneManager->loadGainsFromFileIndBPower("forward.txt", true);
			gainId=sceneManager->registerAntennaGain(gains);
		} else {
			AntennaGain gains=sceneManager->loadGainsFromFileIndBPower("backward.txt", true);
			gainId=sceneManager->registerAntennaGain(gains);
		}
		//sceneManager->registerTransmitterGain(0,gainId);
	}


	uint nrx=250;	
	uint expand=4;
	if (useRDN) {
		nrx=1024;

		//nrx=8000;
		expand=1024/(nrx/mimorx); //Total 1024 positions to make it divisible by mimorx
		sphereDelta=2e-3;
	} else {
		sphereDelta=2e-4;
	}	
	for (int i=0;i<nrx;++i) {
		//sceneManager->addReceiver(i,make_float3(3.72f,2.7f-6.0f, 20.0f),polarization, sphereRadius, sceneManager->printPower);
		sceneManager->addReceiver(i,rx,polarizationRx, sphereRadius, sceneManager->printPower);
		if (this->useAntennaGain && !emptyTunnel) {
			sceneManager->registerReceiverGain(i,gainId);
		}
	}



	sceneManager->finishSceneContext();

	int rayD=10000;	
	if (useRDN) {
		if (!sectorized) {
			if (half) { 
				std::cout <<"**** Anvers Tunnel with  Half Sphere RDN ***"<<std::endl;	
				sceneManager->setRayRange(0.0,180.0,-90.0,90.0,rayD,rayD);
				dynamic_cast<RayDensityNormalizationSimulation*>(sim)->setInitialDensity(((float)sceneManager->getRaySphere().rayCount)/(2*M_PIf));
			} else {
				std::cout <<"**** Anvers Tunnel with Isotropic RDN ***"<<std::endl;	
				sceneManager->setRayRange(0.0,180.0,0.0,360.0,rayD,rayD);
				dynamic_cast<RayDensityNormalizationSimulation*>(sim)->setInitialDensity(((float)sceneManager->getRaySphere().rayCount)/(4*M_PIf));
			}

			dynamic_cast<RayDensityNormalizationSimulation*>(sim)->setFiltering(filtering);
		}
	} else {
		sceneManager->createRaySphere2D(0,0.1,180,-90,0.1,90);
	}



	//Multitransmitter test

	//Enable multitransmitter
	//if (multitransmitter) {
	//	sceneManager->enableMultitransmitter();
	//}
	//	sceneManager->getTransmitterManager()->registerTransmitter(3,tx, polarizationTx, 1.0f);
	//	sceneManager->getTransmitterManager()->addTransmitterToGroup(3,1.0f, tx, polarizationTx);
	//	sceneManager->getTransmitterManager()->registerTransmitter(4,make_float3(1,-3.8, 0), polarizationTx, 1.0f);
	//	sceneManager->getTransmitterManager()->addTransmitterToGroup(4,1.0f,make_float3(1,-3.8, 0), polarizationTx);
	////Add transmitters
	////for (int i=0; i<mimotx; ++i) {
	////	float3 tx_pos=tx+make_float3((-3.5+i)*delta_tx,0.0f,0.0f);
	////	sceneManager->getTransmitterManager()->registerTransmitter(i+nrx,tx_pos, polarizationTx, 1.0f);
	////	sceneManager->getTransmitterManager()->addTransmitterToGroup(i+nrx,1.0f, tx_pos, polarizationTx);
	////}

	//	sceneManager->addReceiver(1,make_float3(3.4,(2.7-6.0), 150),polarizationRx,sphereRadius, sceneManager->printPower);
	//	sceneManager->addReceiver(2,make_float3(3.2,(2.7-6.0), 150),polarizationRx,sphereRadius, sceneManager->printPower);
	//	//	if (this->useAntennaGain) {
	//	//		sceneManager->registerReceiverGain(1,gainId);
	//	//	}
	//	//sceneManager->transmit(0, 1.0f, tx, polarizationTx, false);
	//	//sceneManager->transmit(0, 1.0f, make_float3(1,-3.8, 0), polarizationTx, false);
	//	ResultReport* report=sceneManager->groupTransmit(false);
	//	report->sortByRx();
	//	std::cout<<"Printing report" <<std::endl;
	//	std::cout<<report->toString()<<std::endl;
	//	delete report;
	//return;	
	//Free lane
	//uint nrx=380;	

	timer.start();

	float zStep=0.951f;
	uint launches=0;
	float tl=0.0f;
	ResultReport report;
	if (multitransmitter) {
		//Enable multitransmitter
		sceneManager->enableMultitransmitter();

		//Add transmitters
		for (int i=0; i<mimotx; ++i) {
				float3 tx_pos=tx+make_float3((-3.5+i)*delta_tx,0.0f,0.0f);
				sceneManager->getTransmitterManager()->registerTransmitter(i+nrx,tx_pos, polarizationTx, 1.0f);
				sceneManager->getTransmitterManager()->addTransmitterToGroup(i+nrx,1.0f, tx_pos, polarizationTx);
			
		}
		//float3 tx_pos=tx+make_float3((-3.5+tt)*delta_tx,0.0f,0.0f);
		for (int j=0;j<mimorx;++j) {
			float zinit=50.0f;
			float radius=sphereRadius;
			Timer tlaunch;
			tlaunch.start();
			//for (int i=0;i<nrx;++i) {
			for (int i=0;i<nrx;++i) {
				float3 posrx;
				//Cars on Free lane
				posrx=make_float3(2.55f+(-3.5+j)*delta_rx,(2.7f-6.0f),zinit );
				if (increaseRadius) {
					sceneManager->updateReceiver(i, posrx,radius);
				} else {
					sceneManager->updateReceiver(i, posrx);
				}
				zinit += zStep;
				radius += sphereDelta;
			}
			if (sectorized) {
				SphereScanConfiguration c;
				//throw  opal::Exception("AnversTests::runAnversFlatStaticBetweenTrucksRDNIsotropic(): not implemented sectorized for LPFlatMeshReflectionSimulation yet");
				tl +=runSectorizedLaunch(c, 0, 1.0, tx, polarizationTx, sim);
			} else {
				//First launch
				//sceneManager->transmit(0, 1.0f, postx, polarizationTx, false);
				//sceneManager->transmit(0, 1.0f, tx, polarizationTx, false);
				ResultReport* rp=sceneManager->groupTransmit(false);
				std::cout<<"Merging reports"<<std::endl;
				report.merge(*rp);
				delete rp;
				//sceneManager->transmit(tt+rxid, 1.0f, tx_pos, polarizationTx, false);
				for (int i=0; i<mimotx; ++i) {
					float3 tx_pos=tx+make_float3((-3.5+i)*delta_tx,0.0f,0.0f);
					sceneManager->getTransmitterManager()->addTransmitterToGroup(i+nrx,1.0f, tx_pos, polarizationTx);
				}
			}

			++launches;
			tlaunch.stop();
			std::cout<<launches<<";z="<<zinit<<";Time launch="<<tlaunch.getTime()<<std::endl;
		}
		/*float zinit=50.0f;
		  for (int k=0;k<expand;++k) {
		  int rxid=0;
		  Timer tlaunch;
		  tlaunch.start();
		//for (int i=0;i<nrx;++i) {
		for (int i=0;i<floor(nrx/mimorx);++i) {
		float3 posrx;
		for (int j=0;j<mimorx;++j) {
		//Cars on Free lane
		posrx=make_float3(2.55f+(-3.5+j)*delta_rx,(2.7f-6.0f),zinit );
		if (increaseRadius) {
		sceneManager->updateReceiver(rxid, posrx,sphereRadius);
		} else {
		sceneManager->updateReceiver(rxid, posrx);
		}
		++rxid;
		}
		zinit += zStep;
		sphereRadius += sphereDelta;
		}
		if (sectorized) {
		SphereScanConfiguration c;
		//throw  opal::Exception("AnversTests::runAnversFlatStaticBetweenTrucksRDNIsotropic(): not implemented sectorized for LPFlatMeshReflectionSimulation yet");
		tl +=runSectorizedLaunch(c, 0, 1.0, tx, polarizationTx, sim);
		} else {
		//First launch
		//sceneManager->transmit(0, 1.0f, postx, polarizationTx, false);
		//sceneManager->transmit(0, 1.0f, tx, polarizationTx, false);
		ResultReport* rp=sceneManager->groupTransmit(false);
		std::cout<<"Merging reports"<<std::endl;
		report.merge(*rp);
		delete rp;
		//sceneManager->transmit(tt+rxid, 1.0f, tx_pos, polarizationTx, false);
		for (int i=0; i<mimotx; ++i) {
		float3 tx_pos=tx+make_float3((-3.5+i)*delta_tx,0.0f,0.0f);
		sceneManager->getTransmitterManager()->addTransmitterToGroup(i+nrx,1.0f, tx_pos, polarizationTx);
		}
		}

		++launches;
		tlaunch.stop();
		std::cout<<launches<<";z="<<zinit<<";Time launch="<<tlaunch.getTime()<<std::endl;
		}*/
		report.sort();
		} else {
			for (int tt=0; tt<mimotx; ++tt) {
				float3 tx_pos=tx+make_float3((-3.5+tt)*delta_tx,0.0f,0.0f);
				Timer tlaunch;
				for (int j=0;j<mimorx;++j) {
					float zinit=50.0f;
					float radius=sphereRadius;
					tlaunch.start();
					for (int i=0;i<nrx;++i) {
						float3 posrx;
						//Cars on Free lane
						posrx=make_float3(2.55f+(-3.5+j)*delta_rx,(2.7f-6.0f),zinit );
						if (increaseRadius) {
							sceneManager->updateReceiver(i, posrx,radius);
							radius += sphereDelta;
						} else {
							sceneManager->updateReceiver(i, posrx);
						}
						zinit += zStep;
					}
					if (sectorized) {
						SphereScanConfiguration c;
						//throw  opal::Exception("AnversTests::runAnversFlatStaticBetweenTrucksRDNIsotropic(): not implemented sectorized for LPFlatMeshReflectionSimulation yet");
						tl +=runSectorizedLaunch(c, 0, 1.0, tx, polarizationTx, sim);
					} else {
						//First launch
						//sceneManager->transmit(0, 1.0f, postx, polarizationTx, false);
						//sceneManager->transmit(0, 1.0f, tx, polarizationTx, false);
						//	sceneManager->groupTransmit(false);
						ResultReport* rp=sceneManager->transmit(tt+nrx, 1.0f, tx_pos, polarizationTx, false);

						std::cout<<"Merging reports"<<std::endl;
						report.merge(*rp);
						delete rp;
					}
				}

				++launches;
				tlaunch.stop();
				std::cout<<launches<<";tx="<<tt<<";Time launch="<<tlaunch.getTime()<<std::endl;
			}

		} 
		//        byRxZ sf;
		//	report.template sortBy<byRxZ>(sf);
		if (sectorized) {
			std::cout<<"Total time="<<tl<<". Time/launch="<<(tl/launches)<<std::endl;
		} else {	
			timer.stop();
			std::cout<<"Time="<<timer.getTime()<<". Time/launch="<<(timer.getTime()/launches)<<std::endl;
		}
		std::ofstream file("mimos6.txt");
		file<<report.toString();
		file.close();
		}

void AnversMimo::runFlatTrucks(bool half, bool computeField, bool useReflection, bool useDiffraction, bool sectorized, bool forward, bool multitransmitter, std::string filePath) {
	Timer totalTime;
	totalTime.start();
	OpalSimulation* sim;
	if (useRDN) {
		sim = new RayDensityNormalizationSimulation(sceneManager);
		//RayDensityNormalizationSimulation* sim = new RayDensityNormalizationSimulation(sceneManager);
	} else {	
		sim= new LPFlatMeshReflectionSimulation(sceneManager);
		//LPFlatMeshReflectionSimulation* sim= new LPFlatMeshReflectionSimulation(sceneManager);
	}

	//sim->setExecutionMethod(RDNExecutionMode::HITINFO);
	//sim->setEnableTraceLog(true);
	//sim->setPrintHits(true);
	sceneManager->setSimulation(sim);
	ComputeMode mode=ComputeMode::VOLTAGE;
	if (computeField) {
		std::cout<<"Computing FIELD" <<std::endl;
		mode=ComputeMode::FIELD;
	}
	sim->setComputeMode(mode);
	sim->setEnableSimulation(useReflection);
	sceneManager->enableGenerateRaysOnLaunch();	
	//Add diffraction simulation
	SingleDiffraction* simd= new SingleDiffraction(sceneManager);
	sceneManager->setSimulation(simd);
	simd->setComputeMode(mode);
	//simd->setEnableTraceLog(true);

	if (this->useAntennaGain) {
		sceneManager->setUseAntennaGain(true);
	}
	sceneManager->initContext(frequency);
	simd->setEnableSimulation(useDiffraction);
	if (emptyTunnel) {
		simd->setEnableSimulation(false);
	}
	//simd->setPrintHits(true);
	//sceneManager->getSimulation(0)->setPrintHits(true);


	Timer timer;
	std::cout<<"\tf="<<(frequency/1e6)<< " MHz; tx=" <<tx<<"; polarization="<<polarizationTx<<std::endl;
	std::cout<<"\trx="<<rx<<"; radius="<<sphereRadius;
	if (increaseRadius) {
		std::cout<<"+ "<<sphereDelta;
	}
	std::cout<<"; polarization="<<polarizationRx<<std::endl;
	sceneManager->setMinEpsilon(1e-4f);
	//std::string path("meshes/anvers/trucks");
	//loadAnversScenario(path);


	Timer scene;
	scene.start();
	//Flat with trucks
	//std::string path("straight-trucks-metal-one.json");
	//std::string path("straight-trucks-metal-4.json");
	//loadStraightAnversTunnel();	
	if (!emptyTunnel) {
		std::string path("anvers/straight-trucks-metal.json");
		//std::string path("straight-single-truck.json");
		loadAnversJsonScenario(path);
	}
	std::cout<<"Loading tunnel walls"<<std::endl;
	float width=10.2f;
	float height=6.0f;
	float length=1070.0f;
	Matrix4x4 tm;
	tm.setRow(0, make_float4(width, 0, 0, width/2));
	tm.setRow(1, make_float4(0, height, 0, -height/2));
	tm.setRow(2, make_float4(0, 0, length, length/2.0f));
	tm.setRow(3, make_float4(0, 0, 0, 1));
	MaterialEMProperties emProp1;
	emProp1.dielectricConstant = make_float2(5.0f, -60.0f*sceneManager->getChannelParameters().waveLength*0.01f);
	emProp1.tattenuation = make_float2(0.1f,-75.f );
	loadTransformedSquareTunnel(emProp1, tm);
	scene.stop();
	std::cout<<"Creating scene time="<<scene.getTime()<<std::endl;

	int gainId;	
	if (this->useAntennaGain && !emptyTunnel) {
		if (forward) {
			AntennaGain gains=sceneManager->loadGainsFromFileIndBPower("forward.txt", true);
			gainId=sceneManager->registerAntennaGain(gains);
		} else {
			AntennaGain gains=sceneManager->loadGainsFromFileIndBPower("backward.txt", true);
			gainId=sceneManager->registerAntennaGain(gains);
		}
		//sceneManager->registerTransmitterGain(0,gainId);
	}


	uint nrx=20;	
	uint expand=4;
	float zStep=0.951f;
	float initZ=50;
	float endZ=1000;
	//int positions=(endZ-initZ)/zStep;
	int positions=1000;
	if (useRDN) {
		nrx=positions*mimorx;
		std::cout<<"positions="<<positions<<"mimorx="<<mimorx<<"nrx="<<nrx<<std::endl;
		//nrx=8000;
		//expand=1024/(nrx/mimorx); //Total 1024 positions to make it divisible by mimorx
		sphereDelta=2e-3;
	} else {
		sphereDelta=2e-4;
	}
 	Timer rec;
	rec.start();
	//Reuse polarization matrix for faster loading of large number of receivers
	optix::Matrix<4,4> polMatrix=sceneManager->computeMatrixFromWorldToPolarization(polarizationRx); 
		
	for (int i=0;i<nrx;++i) {
		//sceneManager->addReceiver(i,make_float3(3.72f,2.7f-6.0f, 20.0f),polarization, sphereRadius, sceneManager->printPower);
		sceneManager->addReceiver(i,rx,polarizationRx, sphereRadius, sceneManager->printPower, &polMatrix);
		if (this->useAntennaGain && !emptyTunnel) {
			sceneManager->registerReceiverGain(i,gainId);
		}
	}
	rec.stop();
	std::cout<<"Adding receivers time="<<rec.getTime()<<std::endl;



	sceneManager->finishSceneContext();

	int rayD=10000;	
	if (useRDN) {
		if (!sectorized) {
			if (half) { 
				std::cout <<"**** Anvers Tunnel with  Half Sphere RDN ***"<<std::endl;	
				sceneManager->setRayRange(0.0,180.0,-90.0,90.0,rayD,rayD);
				dynamic_cast<RayDensityNormalizationSimulation*>(sim)->setInitialDensity(((float)sceneManager->getRaySphere().rayCount)/(2*M_PIf));
			} else {
				std::cout <<"**** Anvers Tunnel with Isotropic RDN ***"<<std::endl;	
				sceneManager->setRayRange(0.0,180.0,0.0,360.0,rayD,rayD);
				dynamic_cast<RayDensityNormalizationSimulation*>(sim)->setInitialDensity(((float)sceneManager->getRaySphere().rayCount)/(4*M_PIf));
			}

			dynamic_cast<RayDensityNormalizationSimulation*>(sim)->setFiltering(filtering);
		}
	} else {
		sceneManager->createRaySphere2D(0,0.1,180,-90,0.1,90);
	}




	timer.start();

	uint launches=0;
	float tl=0.0f;
	ResultReport report;
	if (multitransmitter) {
		//Enable multitransmitter
		sceneManager->enableMultitransmitter();
		Timer tlaunch;

		//Register all transmitters
		int idtx=0;
		for (int i=0; i<mimotx; ++i) {
			for (int k=0; k<mimotx; ++k) {
				float3 tx_pos=tx+make_float3((-3.5+i)*delta_tx,(-3.5+k)*delta_tx,0.0f);
				sceneManager->getTransmitterManager()->registerTransmitter(idtx+nrx,tx_pos, polarizationTx, 1.0f);
				//sceneManager->getTransmitterManager()->addTransmitterToGroup(idtx+nrx,1.0f, tx_pos, polarizationTx);
				++idtx;
			}
		}
		
		//given a buffer(mimorx, positions) for a (x,y) the linear index i= y * width + x
		for (int x=0;x<mimorx;++x) {
			float zinit=0.0f;
			float radius=sphereRadius;
			for (int y=0;y<positions;++y) {
				float3 posrx;
				//Cars on Free lane
				posrx=rx+make_float3((-3.5+x)*delta_rx,0.0f,zinit );
				int idrx=y*mimorx+x;
				if (increaseRadius) {
					sceneManager->updateReceiver(idrx, posrx,radius);
					radius += sphereDelta;
				} else {
					sceneManager->updateReceiver(idrx, posrx);
				}
				zinit += zStep;
			}
		}
		if (sectorized) {
			SphereScanConfiguration c;
			//throw  opal::Exception("AnversTests::runAnversFlatStaticBetweenTrucksRDNIsotropic(): not implemented sectorized for LPFlatMeshReflectionSimulation yet");
			tl +=runSectorizedLaunch(c, 0, 1.0, tx, polarizationTx, sim);
		} else {
			//First launch
			//sceneManager->transmit(0, 1.0f, postx, polarizationTx, false);
			//sceneManager->transmit(0, 1.0f, tx, polarizationTx, false);
			//sceneManager->transmit(tt+rxid, 1.0f, tx_pos, polarizationTx, false);
			int idtx=0;
			for (int i=0; i<mimotx; ++i) {
				tlaunch.start();
				for (int k=0; k<mimotx; ++k) {
					float3 tx_pos=tx+make_float3((-3.5+i)*delta_tx,(-3.5+k)*delta_tx,0.0f);
					sceneManager->getTransmitterManager()->addTransmitterToGroup(idtx+nrx,1.0f, tx_pos, polarizationTx);
					++idtx;
				}
				ResultReport* rp=sceneManager->groupTransmit(false);
				std::cout<<"Merging reports"<<std::endl;
				report.merge(*rp);
				delete rp;
				++launches;
				tlaunch.stop();
				std::cout<<launches<<";Time launch="<<tlaunch.getTime()<<std::endl;
			}
		}


		report.sort();
	} else {
		int idtx=0;
		for (int tt=0; tt<mimotx; ++tt) {
			for (int ff=0; ff<mimotx; ++ff) {
				float3 tx_pos=tx+make_float3((-3.5+tt)*delta_tx,(-3.5+ff)*delta_tx,0.0f);
				Timer tlaunch;
				//given a buffer(mimorx, positions) for a (x,y) the linear index i= y * width + x
				for (int x=0;x<mimorx;++x) {
					float zinit=0.0f;
					float radius=sphereRadius;
					tlaunch.start();
					for (int y=0;y<positions;++y) {
						float3 posrx;
						//Cars on Free lane
						posrx=rx+make_float3((-3.5+x)*delta_rx,0.0f,zinit );
						int idrx=y*mimorx+x;
						if (increaseRadius) {
							sceneManager->updateReceiver(idrx, posrx,radius);
							radius += sphereDelta;
						} else {
							sceneManager->updateReceiver(idrx, posrx);
						}
						zinit += zStep;
					}
				}
				if (sectorized) {
					SphereScanConfiguration c;
					//throw  opal::Exception("AnversTests::runAnversFlatStaticBetweenTrucksRDNIsotropic(): not implemented sectorized for LPFlatMeshReflectionSimulation yet");
					tl +=runSectorizedLaunch(c, 0, 1.0, tx, polarizationTx, sim);
				} else {
					//First launch
					//sceneManager->transmit(0, 1.0f, postx, polarizationTx, false);
					//sceneManager->transmit(0, 1.0f, tx, polarizationTx, false);
					//	sceneManager->groupTransmit(false);
					ResultReport* rp=sceneManager->transmit(idtx+nrx, 1.0f, tx_pos, polarizationTx, false);
					++idtx;
					std::cout<<"Merging reports"<<std::endl;
					report.merge(*rp);
					delete rp;
				}
				++launches;
				tlaunch.stop();
				std::cout<<launches<<";tx="<<idtx<<";Time launch="<<tlaunch.getTime()<<std::endl;
			}

		}
	} //else 
	//        byRxZ sf;
	//	report.template sortBy<byRxZ>(sf);
	if (sectorized) {
		std::cout<<"Total time="<<tl<<". Time/launch="<<(tl/launches)<<std::endl;
	} else {	
		timer.stop();
		std::cout<<"Time="<<timer.getTime()<<". Time/launch="<<(timer.getTime()/launches)<<std::endl;
	}
	//std::ofstream file("mimos6.txt");
	//file<<report.toString();
	//file.close();
	std::cout<<"Saving results to "<<filePath<<std::endl;
	report.toCSV(filePath);
	totalTime.stop();
	std::cout<<"TotalTime="<<totalTime.getTime()<<"; setup time="<<(totalTime.getTime()-timer.getTime())<<std::endl;
}
void AnversMimo::runFlatTrucksLille(bool half, bool computeField, bool useReflection, bool useDiffraction, bool sectorized, bool forward, bool multitransmitter, std::string filePath) {
	Timer totalTime;
	totalTime.start();
	OpalSimulation* sim;
	if (useRDN) {
		sim = new RayDensityNormalizationSimulation(sceneManager);
		//RayDensityNormalizationSimulation* sim = new RayDensityNormalizationSimulation(sceneManager);
	} else {	
		sim= new LPFlatMeshReflectionSimulation(sceneManager);
		//LPFlatMeshReflectionSimulation* sim= new LPFlatMeshReflectionSimulation(sceneManager);
	}

	//sim->setExecutionMethod(RDNExecutionMode::HITINFO);
	//sim->setEnableTraceLog(true);
	//sim->setPrintHits(true);
	sceneManager->setSimulation(sim);
	ComputeMode mode=ComputeMode::VOLTAGE;
	if (computeField) {
		std::cout<<"Computing FIELD" <<std::endl;
		mode=ComputeMode::FIELD;
	}
	sim->setComputeMode(mode);
	sim->setEnableSimulation(useReflection);
	sceneManager->enableGenerateRaysOnLaunch();	
	//Add diffraction simulation
	SingleDiffraction* simd= new SingleDiffraction(sceneManager);
	sceneManager->setSimulation(simd);
	simd->setComputeMode(mode);
	simd->setEnableSimulation(useDiffraction);
	//simd->setEnableTraceLog(true);

	if (this->useAntennaGain) {
		sceneManager->setUseAntennaGain(true);
	}
	sceneManager->initContext(frequency);
	if (emptyTunnel) {
		simd->setEnableSimulation(false);
	}
	//simd->setPrintHits(true);
	//sceneManager->getSimulation(0)->setPrintHits(true);


	Timer timer;
	sceneManager->setMinEpsilon(1e-4f);
	//std::string path("meshes/anvers/trucks");
	//loadAnversScenario(path);


	Timer scene;
	scene.start();
	//Flat with trucks
	//std::string path("straight-trucks-metal-one.json");
	//std::string path("straight-trucks-metal-4.json");
	//loadStraightAnversTunnel();	
	if (!emptyTunnel) {
		std::string path("anvers/straight-trucks-metal.json");
		//std::string path("anvers/trucks-random.json");
		//std::string path("straight-single-truck.json");
		loadAnversJsonScenario(path);
	}
	std::cout<<"Loading tunnel walls"<<std::endl;
	float width=10.2f;
	float height=6.0f;
	float length=1070.0f;
	Matrix4x4 tm;
	tm.setRow(0, make_float4(width, 0, 0, width/2));
	tm.setRow(1, make_float4(0, height, 0, -height/2));
	tm.setRow(2, make_float4(0, 0, length, length/2.0f));
	tm.setRow(3, make_float4(0, 0, 0, 1));
	MaterialEMProperties emProp1;
	emProp1.dielectricConstant = make_float2(5.0f, -60.0f*sceneManager->getChannelParameters().waveLength*0.01f);
	emProp1.tattenuation = make_float2(0.1f,-75.f );
	loadTransformedSquareTunnel(emProp1, tm);
	scene.stop();
	std::cout<<"Creating scene time="<<scene.getTime()<<std::endl;

	int gainId;	
	if (this->useAntennaGain && !emptyTunnel) {
		if (forward) {
			AntennaGain gains=sceneManager->loadGainsFromFileIndBPower("anvers/forward.txt", true);
			gainId=sceneManager->registerAntennaGain(gains);
		} else {
			AntennaGain gains=sceneManager->loadGainsFromFileIndBPower("anvers/backward.txt", true);
			gainId=sceneManager->registerAntennaGain(gains);
		}
		//sceneManager->registerTransmitterGain(0,gainId);
	}


	uint nrx=20;	
	uint expand=4;
	float zStep=0.15f; //15 cm
	//float zStep=0.07f; //15 cm
	//float zStep=0.951f; //15 cm
	float initZ=50;
	float endZ=1000;
	//int positions=floor((endZ-initZ)/zStep);
	int positions=1000; //We will do 6 per 1000 receivers, we lose the last 50 meters
	if (useRDN) {
		nrx=positions; //just add the 6333 receivers
		std::cout<<"positions="<<positions<<"mimorx="<<mimorx<<"nrx="<<nrx<<std::endl;
		//nrx=8000;
		//expand=1024/(nrx/mimorx); //Total 1024 positions to make it divisible by mimorx
		sphereDelta=3e-4;
		//sphereDelta=2e-3;
	} else {
		sphereDelta=2e-4;
	}
 	Timer rec;
	rec.start();
	//Reuse polarization matrix for faster loading of large number of receivers
	//optix::Matrix<4,4> polMatrix=sceneManager->computeMatrixFromWorldToPolarization(polarizationRx); 
		
	for (int i=0;i<nrx;++i) {
		sceneManager->addReceiver(i,make_float3(3.72f,2.7f-6.0f, 20.0f),polarizationRx, sphereRadius, sceneManager->printPower);
		//sceneManager->addReceiver(i,rx,polarizationRx, sphereRadius, sceneManager->printPower, &polMatrix);
		if (this->useAntennaGain && !emptyTunnel) {
			sceneManager->registerReceiverGain(i,gainId);
		}
	}
	rec.stop();
	std::cout<<"Adding receivers time="<<rec.getTime()<<std::endl;
	std::cout<<"\tf="<<(frequency/1e6)<< " MHz; tx=" <<tx<<"; polarization="<<polarizationTx<<std::endl;
	std::cout<<"\trx="<<rx<<"; radius="<<sphereRadius;
	if (increaseRadius) {
		std::cout<<"+ "<<sphereDelta;
	}
	std::cout<<"; polarization="<<polarizationRx<<std::endl;



	sceneManager->finishSceneContext();

	int rayD=10000;	
	if (useRDN) {
		if (!sectorized) {
			if (half) { 
				std::cout <<"**** Anvers Tunnel with  Half Sphere RDN ***"<<std::endl;	
				sceneManager->setRayRange(0.0,180.0,-90.0,90.0,rayD,rayD);
				dynamic_cast<RayDensityNormalizationSimulation*>(sim)->setInitialDensity(((float)sceneManager->getRaySphere().rayCount)/(2*M_PIf));
			} else {
				std::cout <<"**** Anvers Tunnel with Isotropic RDN ***"<<std::endl;	
				sceneManager->setRayRange(0.0,180.0,0.0,360.0,rayD,rayD);
				dynamic_cast<RayDensityNormalizationSimulation*>(sim)->setInitialDensity(((float)sceneManager->getRaySphere().rayCount)/(4*M_PIf));
			}

			dynamic_cast<RayDensityNormalizationSimulation*>(sim)->setFiltering(filtering);
		}
	} else {
		sceneManager->createRaySphere2D(0,0.1,180,-90,0.1,90);
	}




	timer.start();

	uint launches=0;
	float tl=0.0f;
	ResultReport report;
	int idtx=0;
//	for (int tt=0; tt<1; ++tt) {
	for (int tt=0; tt<mimotx_x; ++tt) {
		//for (int ff=0; ff<1; ++ff) {
		for (int ff=0; ff<mimotx_y; ++ff) {
			float3 tx_pos=tx+make_float3((-1.5+tt)*delta_tx,(-3.5+ff)*delta_tx,0.0f);
			//given a buffer(mimorx, positions) for a (x,y) the linear index i= y * width + x
			//for (int x=0;x<1;++x) {
			for (int x=0;x<mimorx;++x) {
				Timer batchTimer;
				batchTimer.start();
				float zinit=0.0f;
				float radius=sphereRadius;
				int idrx=0;
				int icr=0;
				for (int batch=0; batch<6; ++batch) {
					Timer tlaunch;
					for (int y=0;y<nrx;++y) {
						float3 posrx;
						//Cars on Free lane
						posrx=rx+make_float3((-1.5+x)*delta_rx,0.0f,zinit );
						//posrx=make_float3(2.28750000000000,	-3.30000000000000, 50+zinit);
						//posrx=make_float3(3.2,	-3.30000000000000, 50+zinit);

						//int idrx=y*mimorx+x;
						if (increaseRadius) {
							sceneManager->updateReceiver(y, posrx,radius);
							radius += sphereDelta;
						} else {
							sceneManager->updateReceiver(y, posrx);
						}
						zinit += zStep;
						++idrx;
					}
					if (sectorized) {
						SphereScanConfiguration c;
						//throw  opal::Exception("AnversTests::runAnversFlatStaticBetweenTrucksRDNIsotropic(): not implemented sectorized for LPFlatMeshReflectionSimulation yet");
						tl +=runSectorizedLaunch(c, 0, 1.0, tx, polarizationTx, sim);
					} else {
						//First launch
						//sceneManager->transmit(0, 1.0f, postx, polarizationTx, false);
						//sceneManager->transmit(0, 1.0f, tx, polarizationTx, false);
						//	sceneManager->groupTransmit(false);
						tlaunch.start();
						//tx_pos=make_float3(0.412500000000000,-3.88750000000000,	0);
						//tx_pos=make_float3(0.5,-3.8,	0);
						std::cout<<"Launching with transmitter at"<<tx_pos<<std::endl;
						ResultReport* rp=sceneManager->transmit(idtx+nrx, 1.0f, tx_pos, polarizationTx, false);
						std::cout<<"Merging reports"<<std::endl;
						report.merge(*rp);
						delete rp;
						++idtx;
						++launches;
						tlaunch.stop();
						std::cout<<launches<<";tx="<<idtx<<";Time launch="<<tlaunch.getTime()<<std::endl;
					}
				}//batch
				 //Completed all 6000 receivers
				batchTimer.stop();
				std::string reportPath=filePath+"-"+std::to_string(x)+"-"+std::to_string(ff)+"-"+std::to_string(tt)+".csv";
				std::cout<<"Saving results to "<<reportPath<<", time="<<batchTimer.getTime()<< std::endl;
				report.toCSV(reportPath);
				report.clear();
			} //mimorx

		} //mimotx_y

		}//mimotx_x
		 //        byRxZ sf;
		 //	report.template sortBy<byRxZ>(sf);
		if (sectorized) {
			std::cout<<"Total time="<<tl<<". Time/launch="<<(tl/launches)<<std::endl;
		} else {	
			timer.stop();
			std::cout<<"Time="<<timer.getTime()<<". Time/launch="<<(timer.getTime()/launches)<<std::endl;
		}
		//std::ofstream file("mimos6.txt");
		//file<<report.toString();
		//file.close();
		totalTime.stop();
	std::cout<<"TotalTime="<<totalTime.getTime()<<"; setup time="<<(totalTime.getTime()-timer.getTime())<<std::endl;
}

void AnversMimo::runFlatCarTrucksLille(bool half, bool computeField, bool useReflection, bool useDiffraction, bool sectorized, bool forward, bool multitransmitter, std::string filePath) {
	Timer totalTime;
	totalTime.start();
	OpalSimulation* sim;
	if (useRDN) {
		sim = new RayDensityNormalizationSimulation(sceneManager);
		//RayDensityNormalizationSimulation* sim = new RayDensityNormalizationSimulation(sceneManager);
	} else {	
		sim= new LPFlatMeshReflectionSimulation(sceneManager);
		//LPFlatMeshReflectionSimulation* sim= new LPFlatMeshReflectionSimulation(sceneManager);
	}

	//sim->setExecutionMethod(RDNExecutionMode::HITINFO);
	//sim->setEnableTraceLog(true);
	//sim->setPrintHits(true);
	sceneManager->setSimulation(sim);
	ComputeMode mode=ComputeMode::VOLTAGE;
	if (computeField) {
		std::cout<<"Computing FIELD" <<std::endl;
		mode=ComputeMode::FIELD;
	}
	sim->setComputeMode(mode);
	sim->setEnableSimulation(useReflection);
	sceneManager->enableGenerateRaysOnLaunch();	
	//Add diffraction simulation
	SingleDiffraction* simd= new SingleDiffraction(sceneManager);
	sceneManager->setSimulation(simd);
	simd->setComputeMode(mode);
	simd->setEnableSimulation(useDiffraction);
	//simd->setEnableTraceLog(true);

	if (this->useAntennaGain) {
		sceneManager->setUseAntennaGain(true);
	}
	sceneManager->initContext(frequency);
	if (emptyTunnel) {
		simd->setEnableSimulation(false);
	}
	//simd->setPrintHits(true);
	//sceneManager->getSimulation(0)->setPrintHits(true);


	Timer timer;
	sceneManager->setMinEpsilon(1e-4f);
	//std::string path("meshes/anvers/trucks");
	//loadAnversScenario(path);


	Timer scene;
	scene.start();
	//Flat with trucks
	//std::string path("straight-trucks-metal-one.json");
	//std::string path("straight-trucks-metal-4.json");
	//loadStraightAnversTunnel();	
	if (!emptyTunnel) {
		std::string path("anvers/trucks-cars.json");
		//std::string path("straight-single-truck.json");
		loadAnversJsonScenario(path);
	}
	std::cout<<"Loading tunnel walls"<<std::endl;
	float width=10.2f;
	float height=6.0f;
	float length=1070.0f;
	Matrix4x4 tm;
	tm.setRow(0, make_float4(width, 0, 0, width/2));
	tm.setRow(1, make_float4(0, height, 0, -height/2));
	tm.setRow(2, make_float4(0, 0, length, length/2.0f));
	tm.setRow(3, make_float4(0, 0, 0, 1));
	MaterialEMProperties emProp1;
	emProp1.dielectricConstant = make_float2(5.0f, -60.0f*sceneManager->getChannelParameters().waveLength*0.01f);
	emProp1.tattenuation = make_float2(0.1f,-75.f );
	loadTransformedSquareTunnel(emProp1, tm);
	scene.stop();
	std::cout<<"Creating scene time="<<scene.getTime()<<std::endl;

	int gainId;	
	if (this->useAntennaGain && !emptyTunnel) {
		if (forward) {
			AntennaGain gains=sceneManager->loadGainsFromFileIndBPower("anvers/forward.txt", true);
			gainId=sceneManager->registerAntennaGain(gains);
		} else {
			AntennaGain gains=sceneManager->loadGainsFromFileIndBPower("anvers/backward.txt", true);
			gainId=sceneManager->registerAntennaGain(gains);
		}
		//sceneManager->registerTransmitterGain(0,gainId);
	}


	uint nrx=20;	
	uint expand=4;
	float zStep=0.15f; //15 cm
	float initZ=50;
	float endZ=1000;

	std::vector<float> zcenter;
	float aux=initZ;
	std::vector<float> zv;
	for (int c=0; c<15; ++c) {
		aux += 30+4.5;
		zcenter.push_back(aux);
		//std::cout<<aux<<std::endl;
		for (int k=-50; k<=50; ++k) {
			zv.push_back(aux+(k*0.15));
			//std::cout<<(aux+(k*0.15))<<std::endl;
		}
		aux+=30;

	}
	//int positions=floor((endZ-initZ)/zStep);
	int positions=zv.size(); //1515 
	if (useRDN) {
		nrx=positions; //just add the 6333 receivers
		std::cout<<"positions="<<positions<<"mimorx="<<mimorx<<"nrx="<<nrx<<std::endl;
		//nrx=8000;
		//expand=1024/(nrx/mimorx); //Total 1024 positions to make it divisible by mimorx
		sphereDelta=2e-3;
	} else {
		sphereDelta=2e-4;
	}
 	Timer rec;
	rec.start();
	std::cout<<"\tf="<<(frequency/1e6)<< " MHz; tx=" <<tx<<"; polarization="<<polarizationTx<<std::endl;
	std::cout<<"\trx="<<rx<<"; radius="<<sphereRadius;
	if (increaseRadius) {
		std::cout<<"+ "<<sphereDelta;
	}
	std::cout<<"; polarization="<<polarizationRx<<std::endl;
	//Reuse polarization matrix for faster loading of large number of receivers
	optix::Matrix<4,4> polMatrix=sceneManager->computeMatrixFromWorldToPolarization(polarizationRx); 
		
	for (int i=0;i<nrx;++i) {
		//sceneManager->addReceiver(i,make_float3(3.72f,2.7f-6.0f, 20.0f),polarization, sphereRadius, sceneManager->printPower);
		sceneManager->addReceiver(i,rx,polarizationRx, sphereRadius, sceneManager->printPower, &polMatrix);
		if (this->useAntennaGain && !emptyTunnel) {
			sceneManager->registerReceiverGain(i,gainId);
		}
	}
	rec.stop();
	std::cout<<"Adding receivers time="<<rec.getTime()<<std::endl;



	sceneManager->finishSceneContext();

	int rayD=10000;	
	if (useRDN) {
		if (!sectorized) {
			if (half) { 
				std::cout <<"**** Anvers Tunnel with  Half Sphere RDN ***"<<std::endl;	
				sceneManager->setRayRange(0.0,180.0,-90.0,90.0,rayD,rayD);
				dynamic_cast<RayDensityNormalizationSimulation*>(sim)->setInitialDensity(((float)sceneManager->getRaySphere().rayCount)/(2*M_PIf));
			} else {
				std::cout <<"**** Anvers Tunnel with Isotropic RDN ***"<<std::endl;	
				sceneManager->setRayRange(0.0,180.0,0.0,360.0,rayD,rayD);
				dynamic_cast<RayDensityNormalizationSimulation*>(sim)->setInitialDensity(((float)sceneManager->getRaySphere().rayCount)/(4*M_PIf));
			}

			dynamic_cast<RayDensityNormalizationSimulation*>(sim)->setFiltering(filtering);
		}
	} else {
		sceneManager->createRaySphere2D(0,0.1,180,-90,0.1,90);
	}




	timer.start();

	uint launches=0;
	float tl=0.0f;
	ResultReport report;
	if (multitransmitter) {
		//Enable multitransmitter
		sceneManager->enableMultitransmitter();
		Timer tlaunch;

		//Register all transmitters
		int idtx=0;
		for (int i=0; i<mimotx; ++i) {
			for (int k=0; k<mimotx; ++k) {
				float3 tx_pos=tx+make_float3((-3.5+i)*delta_tx,(-3.5+k)*delta_tx,0.0f);
				sceneManager->getTransmitterManager()->registerTransmitter(idtx+nrx,tx_pos, polarizationTx, 1.0f);
				//sceneManager->getTransmitterManager()->addTransmitterToGroup(idtx+nrx,1.0f, tx_pos, polarizationTx);
				++idtx;
			}
		}
		
		//given a buffer(mimorx, positions) for a (x,y) the linear index i= y * width + x
		for (int x=0;x<mimorx;++x) {
			float zinit=0.0f;
			float radius=sphereRadius;
			for (int y=0;y<positions;++y) {
				float3 posrx;
				//Cars on Free lane
				posrx=rx+make_float3((-3.5+x)*delta_rx,0.0f,zinit );
				int idrx=y*mimorx+x;
				if (increaseRadius) {
					sceneManager->updateReceiver(idrx, posrx,radius);
					radius += sphereDelta;
				} else {
					sceneManager->updateReceiver(idrx, posrx);
				}
				zinit += zStep;
			}
		}
		if (sectorized) {
			SphereScanConfiguration c;
			//throw  opal::Exception("AnversTests::runAnversFlatStaticBetweenTrucksRDNIsotropic(): not implemented sectorized for LPFlatMeshReflectionSimulation yet");
			tl +=runSectorizedLaunch(c, 0, 1.0, tx, polarizationTx, sim);
		} else {
			//First launch
			//sceneManager->transmit(0, 1.0f, postx, polarizationTx, false);
			//sceneManager->transmit(0, 1.0f, tx, polarizationTx, false);
			//sceneManager->transmit(tt+rxid, 1.0f, tx_pos, polarizationTx, false);
			int idtx=0;
			for (int i=0; i<mimotx; ++i) {
				tlaunch.start();
				for (int k=0; k<mimotx; ++k) {
					float3 tx_pos=tx+make_float3((-3.5+i)*delta_tx,(-3.5+k)*delta_tx,0.0f);
					sceneManager->getTransmitterManager()->addTransmitterToGroup(idtx+nrx,1.0f, tx_pos, polarizationTx);
					++idtx;
				}
				ResultReport* rp=sceneManager->groupTransmit(false);
				std::cout<<"Merging reports"<<std::endl;
				report.merge(*rp);
				delete rp;
				++launches;
				tlaunch.stop();
				std::cout<<launches<<";Time launch="<<tlaunch.getTime()<<std::endl;
			}
		}


		report.sort();
	} else {
		int idtx=0;
		for (int tt=0; tt<mimotx_x; ++tt) {
			for (int ff=0; ff<mimotx_y; ++ff) {
				float3 tx_pos=tx+make_float3((-1.5+tt)*delta_tx,(-3.5+ff)*delta_tx,0.0f);
				//given a buffer(mimorx, positions) for a (x,y) the linear index i= y * width + x
				for (int x=0;x<mimorx;++x) {
					float zinit=0.0f;
					float radius=sphereRadius;
					int idrx=0;
					Timer tlaunch;
					for (int y=0;y<nrx;++y) {
						float3 posrx;
						//Cars on Free lane
						posrx=rx+make_float3((-1.5+x)*delta_rx,0.0f, zv[y] );
						//std::cout<<"zv="<<zv[y]<<"posrx="<<posrx<<std::endl;
						//int idrx=y*mimorx+x;
						if (increaseRadius) {
							sceneManager->updateReceiver(y, posrx,radius);
							radius += sphereDelta;
						} else {
							sceneManager->updateReceiver(y, posrx);
						}
						zinit += zStep;
						++idrx;
					}
					if (sectorized) {
						SphereScanConfiguration c;
						//throw  opal::Exception("AnversTests::runAnversFlatStaticBetweenTrucksRDNIsotropic(): not implemented sectorized for LPFlatMeshReflectionSimulation yet");
						tl +=runSectorizedLaunch(c, 0, 1.0, tx, polarizationTx, sim);
					} else {
						//First launch
						//sceneManager->transmit(0, 1.0f, postx, polarizationTx, false);
						//sceneManager->transmit(0, 1.0f, tx, polarizationTx, false);
						//	sceneManager->groupTransmit(false);
						tlaunch.start();
						std::cout<<"Launching with transmitter "<<tx_pos<<std::endl;
						ResultReport* rp=sceneManager->transmit(idtx+nrx, 1.0f, tx_pos, polarizationTx, false);
						std::cout<<"Merging reports"<<std::endl;
						report.merge(*rp);
						delete rp;
						++idtx;
						++launches;
						tlaunch.stop();
						std::cout<<launches<<";tx="<<idtx<<";Time launch="<<tlaunch.getTime()<<std::endl;
					}
					std::string reportPath=filePath+"-"+std::to_string(x)+"-"+std::to_string(ff)+"-"+std::to_string(tt)+".csv";
					std::cout<<"Saving results to "<<reportPath<< std::endl;
					report.toCSV(reportPath);
					report.clear();
				}
			}

		}
	} //else 
	//        byRxZ sf;
	//	report.template sortBy<byRxZ>(sf);
	if (sectorized) {
		std::cout<<"Total time="<<tl<<". Time/launch="<<(tl/launches)<<std::endl;
	} else {	
		timer.stop();
		std::cout<<"Time="<<timer.getTime()<<". Time/launch="<<(timer.getTime()/launches)<<std::endl;
	}
	//std::ofstream file("mimos6.txt");
	//file<<report.toString();
	//file.close();
	totalTime.stop();
	std::cout<<"TotalTime="<<totalTime.getTime()<<"; setup time="<<(totalTime.getTime()-timer.getTime())<<std::endl;
}

