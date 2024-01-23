/***************************************************************/
//
//Copyright (c) 2019 Esteban Egea-Lopez http://girtel.upct.es/~simulations
//
/**************************************************************/
#include "rouxAP.h"
#include "../timer.h"
#include <memory>
#include <fstream>
#include "../flatSimulation.h"
#include "../curvedMeshSimulation.h"
#include "../curvedFlatMeshSimulation.h"
#include "../rayDensityNormalizationSimulation.h"
#include "../util.h"
using namespace opal;
using namespace optix;
RouxAPTests::RouxAPTests(OpalSceneManager* sceneManager, float sphereRadius, bool useDepolarization) : TunnelsBase(sceneManager,sphereRadius,useDepolarization) {
	this->rx=make_float3(0.0f,0.0f,0.0f);
	this->tx=make_float3(0.0f,0.0f,0.0f);
	this->polarizationTx = V;
	this->polarizationRx = V;
	this->discriminateAngle=2.5;
	this->increaseRadius =true;
	this->sphereDelta = 1e-4;
	this->filtering=2; //Not divide
	this->frequency=900e6;
	this->sphereDistanceError = 1e-7;
}
void RouxAPTests::loadCubeEquivalentTunnel() {
	MaterialEMProperties emProp1;
	emProp1.dielectricConstant = make_float2(5.0f, -60.0f*sceneManager->getChannelParameters().waveLength*0.01f);
	emProp1.tattenuation = make_float2(0.1f,-75.f );
	
	loadSquareTunnel(8.0f, 5.6, 1500, emProp1);
}
void RouxAPTests::loadRouxTunnel() {
	//Roux tunnel is a  partial cylinder: we put an open cylinder and an plane cutting it at some distance from origin
	// So, the tunnel runs on Z axis (front, in Unity). It is centered on the origin, and the entrance to the tunnel is at [0, 0, 0]
	//Load a quad as plane
	
	//Get concrete
	//MaterialEMProperties emProp1 = sceneManager->ITUparametersToMaterial(5.31,0,0.0326,0.8905);
	//else
	MaterialEMProperties emProp1;
	emProp1.dielectricConstant = make_float2(5.0f, -60.0f*sceneManager->getChannelParameters().waveLength*0.01f);
	emProp1.tattenuation = make_float2(0.1f,-75.f );
	
	loadCircularTunnel(4.3f, 1500, emProp1);
	
	//Add floor
	//Horizontal plane (ZX) as quad at origin 
	int quadind[6] = { 0,1,2,1,0,3 };
	optix::float3 quadh[4] = { make_float3(-0.5f,0.0f,-0.5f),make_float3(0.5f,0.f,0.5f) ,make_float3(0.5f,0.f,-0.5f) ,make_float3(-0.5f,0.0f,0.5f) };

//Translate to put center of quad at mid length and displace 1.8 meters below origin	
//Scale 10xlength: even though the plane inside the tunnel is just 8 meters wide
	Matrix4x4 tm;
	tm.setRow(0, make_float4(10, 0, 0, 0.f));
	tm.setRow(1, make_float4(0, 1, 0, -2.0f)); //1.8 below origin. In AP it is 2 m.
	//tm.setRow(1, make_float4(0, 1, 0, -2.f));
	tm.setRow(2, make_float4(0, 0, 1500, 750.f));
	tm.setRow(3, make_float4(0, 0, 0,  1)); 
	sceneManager->addStaticMesh(4, quadh, 6, quadind, tm, emProp1 );
	
	//std::vector<float3> v(std::begin(quadh), std::end(quadh));
	//std::vector<int> ind(std::begin(quadind), std::end(quadind));
	//std::cout<<"Writing quad" <<std::endl;
	//sceneManager->writeMeshToPLYFile("quad2.ply", v,ind, tm);
}
void RouxAPTests::runRDNCube() {
	//float freq = 1.31e9f;
	RayDensityNormalizationSimulation* sim = new RayDensityNormalizationSimulation(sceneManager);
	sceneManager->setSimulation(sim);	
	sceneManager->enableGenerateRaysOnLaunch();	
	sceneManager->initContext(frequency);
	//sceneManager->getSimulation()->setPrintHits(true);


	loadCubeEquivalentTunnel();

	Timer timer;
	//Sphere scanning
	float initElevation=0.0f;
	float initAzimuth=-90.0f;
	float endElevation=180.0f;
	float endAzimuth=90.0f;
	float deltaEl=10.0f;
	float deltaAz=10.0f;
	float asEl=0.01f;
	float asAz=0.01f;
	float overlap=0.5f;
	//float overlap=0.0f;
	std::cout<<"**** Simulating Roux AP Cube tunnel with sectorized RDN ***** "<<std::endl;
	//std::cout<<"Transmitting at "<<(frequency/1e6)<< " MHz from" <<tx<<" with polarization="<<polarizationTx<<std::endl;
	//std::cout<<"Receiving with radius="<<sphereRadius<<" with polarization="<<polarizationRx<<std::endl;
	//std::cout<<"Scanning the sphere with ASElevation="<<asEl<< " and ASAzimuth="<<asAz<<std::endl;
	std::cout<<"\tf="<<(frequency/1e6)<< " MHz; tx=" <<tx<<"; polarization="<<polarizationTx<<std::endl;
	std::cout<<"\trx="<<rx<<"; radius="<<sphereRadius;
	if (increaseRadius) {
		std::cout<<"+ "<<sphereDelta;
	}
	std::cout<<"; polarization="<<polarizationRx<<std::endl;
	std::cout<<"\t Scanning the sphere with ASElevation="<<asEl<< " and ASAzimuth="<<asAz<<std::endl;
	sceneManager->setMinEpsilon(1e-4f);
	
	float currentElevation=initElevation;
	float currentAzimuth=initAzimuth;
	int rayD=10000;
//	OpalRaySphereGenerator* gen=sceneManager->getRaySphereGenerator();
//	gen->generateRandomUniformOnDevice(currentElevation,currentElevation+deltaEl,currentAzimuth,currentAzimuth+deltaAz,rayD*rayD);
//	sceneManager->createRaySphereFromExternalBuffer(rayD,rayD,gen->getDevicePointer());
//	RayDensityNormalizationSimulation* sim=dynamic_cast<RayDensityNormalizationSimulation*>(sceneManager->getSimulation());
//	sim->setInitialDensity(sceneManager->getRaySphere().rayCount,currentAzimuth,currentAzimuth+deltaAz,currentElevation, currentElevation+deltaEl);
//	sim->setFiltering(filtering);	
	
	//sceneManager->setPrintEnabled(1024*1024*1024, make_uint3(976,552,0));
	sceneManager->finishSceneContext();
	timer.start();
	

	sceneManager->setRayRange(currentElevation,currentElevation+deltaEl,currentAzimuth,currentAzimuth+deltaAz,rayD,rayD);
	//RayDensityNormalizationSimulation* sim=dynamic_cast<RayDensityNormalizationSimulation*>(sceneManager->getSimulation());
	sim->setInitialDensity(sceneManager->getRaySphere().rayCount,currentAzimuth,currentAzimuth+deltaAz,currentElevation, currentElevation+deltaEl);
	sim->setFiltering(filtering);	
	
	//sceneManager->addReceiver(1,make_float3(3.72,-22.4871,569.75),polarization, sphereRadius, sceneManager->printPower);
	//sceneManager->addReceiver(1,make_float3(3.72,-22.4788,570.5),polarizationRx, sphereRadius, sceneManager->printPower);
	
	uint nrx=50;	
	for (int i=1;i<=nrx;++i) {
		//sceneManager->addReceiver(i,make_float3(3.72f,2.7f-6.0f, 20.0f),polarization, sphereRadius, sceneManager->printPower);
		sceneManager->addReceiver(i,rx,polarizationRx, sphereRadius, sceneManager->printPower);
	}
	
	float zinit=50.0f;
	uint launches=0;
	tx.y=-0.8f;
	rx.y=-0.8f;
	
	for (int i=0;i<floor(500/nrx);++i) {
		//Update z
		for (int j=1;j<=nrx;++j) {
			rx.z=zinit;

			
			if (increaseRadius) {
				sceneManager->updateReceiver(j, rx,sphereRadius);
				sphereRadius += sphereDelta;
			} else {
				sceneManager->updateReceiver(j, rx);
			}
			zinit += 2.0f;
		}
		//First launch
		sceneManager->transmit(0, 1.0f, tx, polarizationTx, true);
		//Now loop to fill the solid angle
		currentAzimuth += deltaAz;
		//Trace all elevations
		while (currentElevation<endElevation) {

			//Trace all azimuth	
			while(currentAzimuth<endAzimuth) {

			//	std::cout<<"Tracing angle (el/az)="<<(currentElevation-overlap)<<","<<(currentElevation+deltaEl+overlap)<<"/"<<(currentAzimuth-overlap)<<","<<(currentAzimuth+deltaAz+overlap)<<std::endl;
				sceneManager->setRayRange(currentElevation,currentElevation+deltaEl,currentAzimuth,currentAzimuth+deltaAz,rayD,rayD);
				//gen->generateRandomUniformOnDevice(currentElevation,currentElevation+deltaEl,currentAzimuth,currentAzimuth+deltaAz,rayD*rayD);
				//sceneManager->createRaySphereFromExternalBuffer(rayD,rayD,gen->getDevicePointer());
				sim->setInitialDensity(sceneManager->getRaySphere().rayCount,currentAzimuth,currentAzimuth+deltaAz,currentElevation, currentElevation+deltaEl);
				sceneManager->transmit(0, 1.0f, tx, polarizationTx, true);
				currentAzimuth += deltaAz;
			}
			currentAzimuth=initAzimuth;
			currentElevation += deltaEl;
		}
		sceneManager->endPartialLaunch(1u);

		currentElevation=initElevation;
		currentAzimuth=initAzimuth;
		++launches;

	} 
	
	timer.stop();
	std::cout<<"Time="<<timer.getTime()<<". Time/launch="<<(timer.getTime()/launches)<<std::endl;
	
}
void RouxAPTests::runRDNIsotropic() {
	//float freq = 1.31e9f;
	RayDensityNormalizationSimulation* sim = new RayDensityNormalizationSimulation(sceneManager);
	sceneManager->setSimulation(sim);	
	sceneManager->enableGenerateRaysOnLaunch();	
	

	sceneManager->initContext(frequency);
	//sceneManager->getSimulation()->setPrintHits(true);


	loadRouxTunnel();

	Timer timer;
	std::cout<<"**** Simulating Roux AP tunnel with Isotropic RDN ***** "<<std::endl;
	//std::cout<<"Transmitting at "<<(frequency/1e6)<< " MHz from" <<tx<<" with polarization="<<polarizationTx<<std::endl;
	//std::cout<<"Receiving with radius="<<sphereRadius<<" with polarization="<<polarizationRx<<std::endl;
	//std::cout<<"Scanning the sphere with ASElevation="<<asEl<< " and ASAzimuth="<<asAz<<std::endl;
	std::cout<<"\tf="<<(frequency/1e6)<< " MHz; tx=" <<tx<<"; polarization="<<polarizationTx<<std::endl;
	std::cout<<"\trx="<<rx<<"; radius="<<sphereRadius;
	if (increaseRadius) {
		std::cout<<"+ "<<sphereDelta;
	}
	std::cout<<"; polarization="<<polarizationRx<<std::endl;
	std::cout<<"******* ***** "<<std::endl;
	//std::cout<<"\t Scanning the sphere with ASElevation="<<asEl<< " and ASAzimuth="<<asAz<<std::endl;
	sceneManager->setMinEpsilon(1e-4f);
	
	//int rayD=31622;
	int rayD=10000;
	//int rayD=3;
  
        //Not necessary if ray generation on launch is enabled
	//OpalRaySphereGenerator* gen=sceneManager->getRaySphereGenerator();
	//gen->generateRandomUniformOnDevice(0.0,180.0,-90.0,90.0,rayD*rayD); 
	//gen->generateRandomUniformOnDevice(0.0,90.0,0.0,360.0,rayD*rayD); 
	//sceneManager->createRaySphereFromExternalBuffer(rayD,rayD,gen->getDevicePointer());
	//sim->setInitialDensity(sceneManager->getRaySphere().rayCount/(4*M_PIf));
	//RayDensityNormalizationSimulation* sim=dynamic_cast<RayDensityNormalizationSimulation*>(sceneManager->getSimulation());
	//float dens=sim->setInitialDensity(sceneManager->getRaySphere().rayCount,-90,90,0,180);
	//sim->setFiltering(filtering);	
	
	//sceneManager->setPrintEnabled(1024*1024*1024, make_uint3(976,552,0));
	sceneManager->finishSceneContext();
	timer.start();
	
	sceneManager->setRayRange(0.0,180.0,-90.0,90.0,rayD,rayD);
	//RayDensityNormalizationSimulation* sim=dynamic_cast<RayDensityNormalizationSimulation*>(sceneManager->getSimulation());
	float dens=sim->setInitialDensity(sceneManager->getRaySphere().rayCount,-90,90,0,180);
	sim->setFiltering(filtering);	

	
	//sceneManager->addReceiver(1,make_float3(3.72,-22.4871,569.75),polarization, sphereRadius, sceneManager->printPower);
	//sceneManager->addReceiver(1,make_float3(3.72,-22.4788,570.5),polarizationRx, sphereRadius, sceneManager->printPower);
	
	//uint nrx=1;	
	uint nrx=100;	
	for (int i=1;i<=nrx;++i) {
		//sceneManager->addReceiver(i,make_float3(3.72f,2.7f-6.0f, 20.0f),polarization, sphereRadius, sceneManager->printPower);
		sceneManager->addReceiver(i,rx,polarizationRx, sphereRadius, sceneManager->printPower);
	}
	
	float zinit=50.0f;
	uint launches=0;
	
	for (int i=0;i<floor(700/nrx);++i) {
		//Update z
		for (int j=1;j<=nrx;++j) {
			rx.z=zinit;

			
			if (increaseRadius) {
				//float rad=sim->sphereRadiusForExpectedDirectHits(100,dens,length(tx-rx));
				//sceneManager->updateReceiver(j, rx,rad);
				//sphereRadius = sqrt(length(tx-rx)*sphereDistanceError);
				sceneManager->updateReceiver(j, rx,sphereRadius);
				if (sphereRadius <= 4.3) {
					sphereRadius += sphereDelta;
				}
			} else {
				sceneManager->updateReceiver(j, rx);
			}
			zinit += 2.0f;
		}
		//First launch
		sceneManager->transmit(0, 5.0f, tx, polarizationTx, false);
		++launches;

	} 
	
	timer.stop();
	std::cout<<"Time="<<timer.getTime()<<". Time/launch="<<(timer.getTime()/launches)<<std::endl;
	
}
void RouxAPTests::runRDNFixedRay() {
	//float freq = 1.31e9f;
	RayDensityNormalizationSimulation* sim = new RayDensityNormalizationSimulation(sceneManager);
	sceneManager->setSimulation(sim);	
	sceneManager->enableGenerateRaysOnLaunch(); //No memory used to store the rays: more ray density
	sceneManager->initContext(frequency);
	//sceneManager->getSimulation()->setPrintHits(true);
	sim->setExecutionMethod(RDNExecutionMode::NOATOMIC);


	loadRouxTunnel();
	Timer timer;
	//Sphere scanning
	float initElevation=0.0f;
	float initAzimuth=-90.0f;
	float endElevation=180.0f;
	float endAzimuth=90.0f;
	float deltaEl=10.0f;
	//float deltaAz=10.0f;
	float asEl=0.01f;
	float asAz=0.01f;
	float overlap=0.5f;
	//float overlap=0.0f;
	std::cout<<"**** Simulating Roux AP tunnel with sectorized RDN and Fixed Ray***** "<<std::endl;
	//std::cout<<"Transmitting at "<<(frequency/1e6)<< " MHz from" <<tx<<" with polarization="<<polarizationTx<<std::endl;
	//std::cout<<"Receiving with radius="<<sphereRadius<<" with polarization="<<polarizationRx<<std::endl;
	//std::cout<<"Scanning the sphere with ASElevation="<<asEl<< " and ASAzimuth="<<asAz<<std::endl;
	std::cout<<"\tf="<<(frequency/1e6)<< " MHz; tx=" <<tx<<"; polarization="<<polarizationTx<<std::endl;
	std::cout<<"\trx="<<rx<<"; radius="<<sphereRadius;
	if (increaseRadius) {
		std::cout<<"+ "<<sphereDelta<<" distanceError "<<sphereDistanceError;
	}
	std::cout<<"; polarization="<<polarizationRx<<std::endl;
	sceneManager->setMinEpsilon(1e-4f);
	uint nrx=750;	
	int rayD=2000;
	sim->setFixedRayBuffer(nrx, rayD);	
	sceneManager->finishSceneContext();
	
	float currentElevation=initElevation;
	float currentAzimuth=initAzimuth;
	//int rayD=10000;
	//make constant density
	float deg2rad=M_PIf/180.0f;
	float rayGoal=1e9; //In rays/sr (stereo-radian)
	//int rayD=floor(sqrt(rayGoal*(cosf(deg2rad*currentElevation)-cosf(deg2rad*(currentElevation+deltaEl))))); //Missing here the azimuth part of the density...
	
 	float	deltaAz=rayD*rayD/(deg2rad*rayGoal*(cosf(deg2rad*currentElevation)-cosf(deg2rad*(currentElevation+deltaEl)))); 

	std::cout<<"\t Scanning the sphere with ASElevation="<<asEl<< " and ASAzimuth="<<asAz<<"rays="<<rayD*rayD<<std::endl;
	
	//OpalRaySphereGenerator* gen=sceneManager->getRaySphereGenerator();
	//gen->generateRandomUniformOnDevice(currentElevation,currentElevation+deltaEl,currentAzimuth,currentAzimuth+deltaAz,rayD*rayD);
	//sceneManager->createRaySphereFromExternalBuffer(rayD,rayD,gen->getDevicePointer());
	//RayDensityNormalizationSimulation* sim=dynamic_cast<RayDensityNormalizationSimulation*>(sceneManager->getSimulation());
	//sim->setInitialDensity(sceneManager->getRaySphere().rayCount,currentAzimuth,currentAzimuth+deltaAz,currentElevation, currentElevation+deltaEl);
	//sim->setFiltering(filtering);	
	
	//sceneManager->setPrintEnabled(1024*1024*1024, make_uint3(976,552,0));
	sceneManager->setRayRange(currentElevation,currentElevation+deltaEl,currentAzimuth,currentAzimuth+deltaAz,rayD,rayD);
	sim->setInitialDensity(sceneManager->getRaySphere().rayCount,currentAzimuth,currentAzimuth+deltaAz,currentElevation, currentElevation+deltaEl);
	sim->setFiltering(filtering);	
	
	for (int i=1;i<=nrx;++i) {
		//sceneManager->addReceiver(i,make_float3(3.72f,2.7f-6.0f, 20.0f),polarization, sphereRadius, sceneManager->printPower);
		sceneManager->addReceiver(i,rx,polarizationRx, sphereRadius, sceneManager->printPower);
	}
	timer.start();
	

	
	//sceneManager->addReceiver(1,make_float3(3.72,-22.4871,569.75),polarization, sphereRadius, sceneManager->printPower);
	//sceneManager->addReceiver(1,make_float3(3.72,-22.4788,570.5),polarizationRx, sphereRadius, sceneManager->printPower);
	
	float zinit=50.0f;
	uint launches=0;
	
	for (int i=0;i<floor(750/nrx);++i) {
                Timer t2;
		t2.start();
		//Update z
		for (int j=1;j<=nrx;++j) {
			rx.z=zinit;

			
			if (increaseRadius) {
				//sphereRadius = sqrt(length(tx-rx)*sphereDistanceError);
				sceneManager->updateReceiver(j, rx,sphereRadius);
				sphereRadius += sphereDelta;
			} else {
				sceneManager->updateReceiver(j, rx);
			}
			zinit += 2.0f;
		}
		//First launch
		sceneManager->transmit(0, 5.0f, tx, polarizationTx, true);
		//Now loop to fill the solid angle
		currentAzimuth += deltaAz;
		//Trace all elevations
		while (currentElevation<endElevation) {

			//Trace all azimuth	
			while(currentAzimuth<endAzimuth) {

				//std::cout<<"Tracing angle (el/az)="<<(currentElevation-overlap)<<","<<(currentElevation+deltaEl+overlap)<<"/"<<(currentAzimuth-overlap)<<","<<(currentAzimuth+deltaAz+overlap)<<std::endl;
				//Constant density
			//	rayD=floor(1e9*(cosf(deg2rad*currentElevation)-cosf(deg2rad*(currentElevation+deltaEl))));
				//rayD=floor(sqrt(rayGoal*(cosf(deg2rad*currentElevation)-cosf(deg2rad*(currentElevation+deltaEl)))));
				deltaAz=rayD*rayD/(deg2rad*rayGoal*(cosf(deg2rad*currentElevation)-cosf(deg2rad*(currentElevation+deltaEl)))); //Missing here the azimuth part of the density...
				if ((currentAzimuth+deltaAz)<endAzimuth) {
					sceneManager->setRayRange(currentElevation,currentElevation+deltaEl,currentAzimuth,currentAzimuth+deltaAz,rayD,rayD);
				} else {
					//Make sure that we do not start over
					sceneManager->setRayRange(currentElevation,currentElevation+deltaEl,currentAzimuth,endAzimuth,rayD,rayD);
				}
			//	gen->generateRandomUniformOnDevice(currentElevation,currentElevation+deltaEl,currentAzimuth,currentAzimuth+deltaAz,rayD*rayD);
			//	sceneManager->createRaySphereFromExternalBuffer(rayD,rayD,gen->getDevicePointer());
				sim->setInitialDensity(sceneManager->getRaySphere().rayCount,currentAzimuth,currentAzimuth+deltaAz,currentElevation, currentElevation+deltaEl);
				sceneManager->transmit(0, 5.0f, tx, polarizationTx, true);
				currentAzimuth += deltaAz;
		                ++launches;
			}
			currentAzimuth=initAzimuth;
			currentElevation += deltaEl;
		}
		sceneManager->endPartialLaunch(1u);
		t2.stop();
		std::cout<<"Batch Time="<<t2.getTime()<<std::endl;
		currentElevation=initElevation;
		currentAzimuth=initAzimuth;

	} 
	
	timer.stop();
	std::cout<<"Time="<<timer.getTime()<<". Time/launch="<<(timer.getTime()/launches)<<std::endl;
	
}
void RouxAPTests::runRDN() {
	//float freq = 1.31e9f;
	RayDensityNormalizationSimulation* sim = new RayDensityNormalizationSimulation(sceneManager);
	sceneManager->setSimulation(sim);	
	sceneManager->enableGenerateRaysOnLaunch(); //No memory used to store the rays: more ray density
	
	sceneManager->initContext(frequency);
	//sceneManager->getSimulation()->setPrintHits(true);

	//RayDensityNormalizationSimulation* sim=dynamic_cast<RayDensityNormalizationSimulation*>(sceneManager->getSimulation());
	//sim->setExecutionMethod(RDNExecutionMode::HITINFO);

	loadRouxTunnel();
	Timer timer;
	//Sphere scanning
	float initElevation=0.0f;
	float initAzimuth=-90.0f;
	float endElevation=180.0f;
	float endAzimuth=90.0f;
	float deltaEl=10.0f;
	float deltaAz=10.0f;
	float asEl=0.01f;
	float asAz=0.01f;
	float overlap=0.5f;
	//float overlap=0.0f;
	std::cout<<"**** Simulating Roux AP tunnel with sectorized RDN ***** "<<std::endl;
	//std::cout<<"Transmitting at "<<(frequency/1e6)<< " MHz from" <<tx<<" with polarization="<<polarizationTx<<std::endl;
	//std::cout<<"Receiving with radius="<<sphereRadius<<" with polarization="<<polarizationRx<<std::endl;
	//std::cout<<"Scanning the sphere with ASElevation="<<asEl<< " and ASAzimuth="<<asAz<<std::endl;
	std::cout<<"\tf="<<(frequency/1e6)<< " MHz; tx=" <<tx<<"; polarization="<<polarizationTx<<std::endl;
	std::cout<<"\trx="<<rx<<"; radius="<<sphereRadius;
	if (increaseRadius) {
		std::cout<<"+ "<<sphereDelta<<" distanceError "<<sphereDistanceError;
	}
	std::cout<<"; polarization="<<polarizationRx<<std::endl;
	sceneManager->setMinEpsilon(1e-4f);
	
	float currentElevation=initElevation;
	float currentAzimuth=initAzimuth;
	//int rayD=10000;
	//make constant density
	float deg2rad=M_PIf/180.0f;
	float rayGoal=10e9; //In rays/sr (stereo-radian)
	float solidAngle=deg2rad*deltaAz*(cosf(deg2rad*currentElevation)-cosf(deg2rad*(currentElevation+deltaEl)));
	int rayD=floor(sqrt(rayGoal*solidAngle)); 
	std::cout<<"\t Scanning the sphere with ASElevation="<<asEl<< " and ASAzimuth="<<asAz<<"rays="<<rayD*rayD<<std::endl;
	
	//OpalRaySphereGenerator* gen=sceneManager->getRaySphereGenerator();
	//gen->generateRandomUniformOnDevice(currentElevation,currentElevation+deltaEl,currentAzimuth,currentAzimuth+deltaAz,rayD*rayD);
	//sceneManager->createRaySphereFromExternalBuffer(rayD,rayD,gen->getDevicePointer());
	//RayDensityNormalizationSimulation* sim=dynamic_cast<RayDensityNormalizationSimulation*>(sceneManager->getSimulation());
	//sim->setInitialDensity(sceneManager->getRaySphere().rayCount,currentAzimuth,currentAzimuth+deltaAz,currentElevation, currentElevation+deltaEl);
	//sim->setFiltering(filtering);	
	
	//sceneManager->setPrintEnabled(1024*1024*1024, make_uint3(976,552,0));
	sceneManager->finishSceneContext();
	timer.start();
	
	sceneManager->setRayRange(currentElevation,currentElevation+deltaEl,currentAzimuth,currentAzimuth+deltaAz,rayD,rayD);
	//RayDensityNormalizationSimulation* sim=dynamic_cast<RayDensityNormalizationSimulation*>(sceneManager->getSimulation());
	sim->setInitialDensity(sceneManager->getRaySphere().rayCount,currentAzimuth,currentAzimuth+deltaAz,currentElevation, currentElevation+deltaEl);
	sim->setFiltering(filtering);	

	
	//sceneManager->addReceiver(1,make_float3(3.72,-22.4871,569.75),polarization, sphereRadius, sceneManager->printPower);
	//sceneManager->addReceiver(1,make_float3(3.72,-22.4788,570.5),polarizationRx, sphereRadius, sceneManager->printPower);
	
	uint nrx=700;	
	for (int i=1;i<=nrx;++i) {
		//sceneManager->addReceiver(i,make_float3(3.72f,2.7f-6.0f, 20.0f),polarization, sphereRadius, sceneManager->printPower);
		sceneManager->addReceiver(i,rx,polarizationRx, sphereRadius, sceneManager->printPower);
	}
	
	float zinit=50.0f;
	uint launches=0;
	
	for (int i=0;i<floor(700/nrx);++i) {
                Timer t2;
		t2.start();
		//Update z
		for (int j=1;j<=nrx;++j) {
			rx.z=zinit;

			
			if (increaseRadius) {
				//sphereRadius = sqrt(length(tx-rx)*sphereDistanceError);
				sceneManager->updateReceiver(j, rx,sphereRadius);
				sphereRadius += sphereDelta;
				//if (sphereRadius>2.5) {
				//	sphereRadius=2.5;
				//}
			} else {
				sceneManager->updateReceiver(j, rx);
			}
			zinit += 2.0f;
		}
		//First launch
		sceneManager->transmit(0, 5.0f, tx, polarizationTx, true);
		//Now loop to fill the solid angle
		currentAzimuth += deltaAz;
		//Trace all elevations
		while (currentElevation<endElevation) {

			//Trace all azimuth	
			while(currentAzimuth<endAzimuth) {

				//std::cout<<"Tracing angle (el/az)="<<(currentElevation-overlap)<<","<<(currentElevation+deltaEl+overlap)<<"/"<<(currentAzimuth-overlap)<<","<<(currentAzimuth+deltaAz+overlap)<<std::endl;
				//Constant density
			//	rayD=floor(1e9*(cosf(deg2rad*currentElevation)-cosf(deg2rad*(currentElevation+deltaEl))));
				solidAngle=deg2rad*deltaAz*(cosf(deg2rad*currentElevation)-cosf(deg2rad*(currentElevation+deltaEl)));
				rayD=floor(sqrt(rayGoal*solidAngle)); //Missing here the azimuth part of the density...
				sceneManager->setRayRange(currentElevation,currentElevation+deltaEl,currentAzimuth,currentAzimuth+deltaAz,rayD,rayD);
			//	gen->generateRandomUniformOnDevice(currentElevation,currentElevation+deltaEl,currentAzimuth,currentAzimuth+deltaAz,rayD*rayD);
			//	sceneManager->createRaySphereFromExternalBuffer(rayD,rayD,gen->getDevicePointer());
				sim->setInitialDensity(sceneManager->getRaySphere().rayCount,currentAzimuth,currentAzimuth+deltaAz,currentElevation, currentElevation+deltaEl);
				sceneManager->transmit(0, 5.0f, tx, polarizationTx, true);
				currentAzimuth += deltaAz;
		                ++launches;
			}
			currentAzimuth=initAzimuth;
			currentElevation += deltaEl;
		}
		sceneManager->endPartialLaunch(1u);
		t2.stop();
		std::cout<<"Batch Time="<<t2.getTime()<<std::endl;
		currentElevation=initElevation;
		currentAzimuth=initAzimuth;

	} 
	
	timer.stop();
	std::cout<<"Time="<<timer.getTime()<<". Time/launch="<<(timer.getTime()/launches)<<std::endl;
	
}
void RouxAPTests::runADSingleRay() {
	LPCurvedFlatMeshReflectionSimulation* sim = new LPCurvedFlatMeshReflectionSimulation(sceneManager);
	sceneManager->setSimulation(sim);
	//sceneManager->enableGenerateRaysOnLaunch();	
	
	sceneManager->initContext(frequency);
	sceneManager->getSimulation()->setPrintHits(true);


	loadRouxTunnel();
	//loadCubeEquivalentTunnel();

	Timer timer;
	//Sphere scanning
	float initElevation=0.0f;
	float initAzimuth=-90.0f;
	float endElevation=180.0f;
	float endAzimuth=90.0f;
	float deltaEl=10.0f;
	float deltaAz=10.0f;
	float asEl=0.01f;
	float asAz=0.01f;
	float overlap=0.5f;
	//float overlap=0.0f;
	std::cout<<"**** Simulating Single Ray Roux tunnel with CURVEDFLATWALLS ***** "<<std::endl;
	std::cout<<"Transmitting at "<<(frequency/1e6)<< " MHz from" <<tx<<" with polarization="<<polarizationTx<<std::endl;
	std::cout<<"Receiving with radius="<<sphereRadius<<" with polarization="<<polarizationRx<<std::endl;
	std::cout<<"Scanning the sphere with ASElevation="<<asEl<< " and ASAzimuth="<<asAz<<std::endl;
	sceneManager->setMinEpsilon(1e-4f);
	
	float currentElevation=initElevation;
	float currentAzimuth=initAzimuth;
	
	
	//***Single ray transmit****
	//float3 mRay=make_float3(0.8986926675, -0.4382142425, 0.01788338087);
	//float3 mRay=make_float3(-0.312150329351425f, 0.0343762412667274f, 0.949410557746887f);
	//float3 mRay=normalize(make_float3(1.0f, 1.0f, 1.0f)); 
	//float3 mRay=normalize(make_float3(0.0f, -1.0f, 1.0f)); 
	//float3 mRay=normalize(make_float3(1.0f, 0.0f, 1.0f)); 
	float3 mRay=normalize(make_float3(-0.891520f, 0.573570f, 1.0f)); 
	sceneManager->createRaySphere2D(1,1,&mRay);
	
	//sceneManager->setPrintEnabled(1024*1024*1024, make_uint3(976,552,0));
	sceneManager->finishSceneContext();
	timer.start();
	
	//Angle for separating duplicate rays

	//LPCurvedMeshReflectionSimulation* sim=dynamic_cast<LPCurvedMeshReflectionSimulation*>(sceneManager->getSimulation());
	sim->setMaxAngleForDuplicateRays(discriminateAngle*M_PIf/180.f);
	
	
	
	rx.z=52.0f;
	//rx.x=0.0024;
	//rx.y=2.5112;
	//rx.z=10.9444;
	uint launches=1;
	sceneManager->addReceiver(1,rx,polarizationRx, sphereRadius, sceneManager->printPower);
	sceneManager->transmit(0, 1.0f, tx, polarizationTx, false);
	
	timer.stop();
	std::cout<<"Time="<<timer.getTime()<<". Time/launch="<<(timer.getTime()/launches)<<std::endl;
	

}
void RouxAPTests::runADCube() {
	//LPCurvedFlatMeshReflectionSimulation* sim = new LPCurvedFlatMeshReflectionSimulation(sceneManager);
	//sceneManager->setSimulation(sim);
	LPFlatMeshReflectionSimulation* sim = new LPFlatMeshReflectionSimulation(sceneManager);
	sceneManager->setSimulation(sim);
	
	sceneManager->initContext(frequency);
	//sceneManager->getSimulation()->setPrintHits(true);

	//Cube tunnel is centered, so we have to displace tx and rx to get the same configuration
	loadCubeEquivalentTunnel();

	Timer timer;
	//Receiver polarization
	//optix::float3 polarization = H; 	
	//optix::float3 polarization = pol; 	
	
	//Transmitter
	//optix::float3 postx = make_float3(0.50f,2.52f-6.0f, 0.0f);
	//float3 polarizationTx = H; 
	//float3 polarizationTx = polarization; 
	//Sphere scanning
	float initElevation=0.0f;
	float initAzimuth=-90.0f;
	float endElevation=180.0f;
	float endAzimuth=90.0f;
	float deltaEl=10.0f;
	float deltaAz=10.0f;
	float asEl=0.01f;
	float asAz=0.01f;
	float overlap=0.5f;
	//float overlap=0.0f;
	std::cout<<"**** Simulating Roux Square tunnel with CURVEDFLATWALLS ***** "<<std::endl;
	//std::cout<<"Transmitting at "<<(frequency/1e6)<< " MHz from" <<tx<<" with polarization="<<polarizationTx<<std::endl;
	//std::cout<<"Receiving with radius="<<sphereRadius<<" with polarization="<<polarizationRx<<std::endl;
	//std::cout<<"Scanning the sphere with ASElevation="<<asEl<< " and ASAzimuth="<<asAz<<std::endl;
	std::cout<<"\tf="<<(frequency/1e6)<< " MHz; tx=" <<tx<<"; polarization="<<polarizationTx<<std::endl;
	std::cout<<"\trx="<<rx<<"; radius="<<sphereRadius;
	if (increaseRadius) {
		std::cout<<" + "<<sphereDelta;
	}
	std::cout<<"; polarization="<<polarizationRx<<std::endl;
	std::cout<<"\tScanning the sphere with ASElevation="<<asEl<< " and ASAzimuth="<<asAz<<std::endl;
	std::cout<<"******************** "<<std::endl;
	
	float currentElevation=initElevation;
	float currentAzimuth=initAzimuth;
	//std::cout<<"Tracing angle (el/az)="<<(currentElevation-overlap)<<","<<(currentElevation+deltaEl+overlap)<<"/"<<(currentAzimuth-overlap)<<","<<(currentAzimuth+deltaAz+overlap)<<std::endl;
	
	sceneManager->setMinEpsilon(1e-4f);
	sceneManager->createRaySphere2D(currentElevation-overlap,asEl,currentElevation+deltaEl+overlap,currentAzimuth-overlap,asAz,currentAzimuth+deltaAz+overlap);
	
	//sceneManager->setPrintEnabled(1024*1024*1024, make_uint3(976,552,0));
	sceneManager->finishSceneContext();
	timer.start();
	
	//Angle for separating duplicate rays

	//LPCurvedMeshReflectionSimulation* sim=dynamic_cast<LPCurvedMeshReflectionSimulation*>(sceneManager->getSimulation());
	//sim->setMaxAngleForDuplicateRays(discriminateAngle*M_PIf/180.f);
	
	
	
	
	uint nrx=100;	
	for (int i=1;i<=nrx;++i) {
		//sceneManager->addReceiver(i,make_float3(3.72f,2.7f-6.0f, 20.0f),polarization, sphereRadius, sceneManager->printPower);
		sceneManager->addReceiver(i,rx,polarizationRx, sphereRadius, sceneManager->printPower);
	}
	
	float zinit=50.0f;
	uint launches=0;
	tx.y=-0.8f;
	rx.y=-0.8f;
	for (int i=0;i<floor(500/nrx);++i) {
		//Update z
		for (int j=1;j<=nrx;++j) {

			rx.z=zinit;
			if (increaseRadius) {
				sceneManager->updateReceiver(j, rx,sphereRadius);
				sphereRadius += sphereDelta;
			} else {
				sceneManager->updateReceiver(j, rx);
			}
			zinit += 2.0f;
		}
		//First launch
		sceneManager->transmit(0, 1.0f, tx, polarizationTx, true);

		//Now loop to fill the solid angle
		currentAzimuth += deltaAz;
		//Trace all elevations
		while (currentElevation<endElevation) {

			//Trace all azimuth	
			while(currentAzimuth<endAzimuth) {
				//std::cout<<"Tracing angle (el/az)="<<(currentElevation-overlap)<<","<<(currentElevation+deltaEl+overlap)<<"/"<<(currentAzimuth-overlap)<<","<<(currentAzimuth+deltaAz+overlap)<<std::endl;
				sceneManager->createRaySphere2D(currentElevation-overlap,asEl,currentElevation+deltaEl+overlap,currentAzimuth-overlap,asAz,currentAzimuth+deltaAz+overlap);
				sceneManager->transmit(0, 1.0f, tx, polarizationTx, true);
				currentAzimuth += deltaAz;
			}
			currentAzimuth=initAzimuth;
			currentElevation += deltaEl;
		}
		sceneManager->endPartialLaunch(1u);

		currentElevation=initElevation;
		currentAzimuth=initAzimuth;
		++launches;

	} 
	
	timer.stop();
	std::cout<<"Time="<<timer.getTime()<<". Time/launch="<<(timer.getTime()/launches)<<std::endl;
	
}
void RouxAPTests::runAD() {
	LPCurvedFlatMeshReflectionSimulation* sim = new LPCurvedFlatMeshReflectionSimulation(sceneManager);
	sceneManager->setSimulation(sim);
	sceneManager->enableGenerateRaysOnLaunch();	
	
	sceneManager->initContext(frequency);
	//sceneManager->getSimulation()->setPrintHits(true);


	loadRouxTunnel();
	Timer timer;
	//Sphere scanning
	float initElevation=0.0f;
	float initAzimuth=-90.0f;
	float endElevation=180.0f;
	float endAzimuth=90.0f;
	float deltaEl=10.0f;
	float deltaAz=10.0f;
	float asEl=0.01f;
	float asAz=0.01f;
	float overlap=0.5f;
	//float overlap=0.0f;
	std::cout<<"**** Simulating Roux tunnel with CURVEDFLATWALLS ***** "<<std::endl;
	//std::cout<<"Transmitting at "<<(frequency/1e6)<< " MHz from" <<tx<<" with polarization="<<polarizationTx<<std::endl;
	//std::cout<<"Receiving with radius="<<sphereRadius<<" with polarization="<<polarizationRx<<std::endl;
	//std::cout<<"Scanning the sphere with ASElevation="<<asEl<< " and ASAzimuth="<<asAz<<std::endl;
	std::cout<<"\tf="<<(frequency/1e6)<< " MHz; tx=" <<tx<<"; polarization="<<polarizationTx<<std::endl;
	std::cout<<"\trx="<<rx<<"; radius="<<sphereRadius;
	if (increaseRadius) {
		std::cout<<" + "<<sphereDelta;
	}
	std::cout<<"; polarization="<<polarizationRx<<std::endl;
	std::cout<<"\tScanning the sphere with ASElevation="<<asEl<< " and ASAzimuth="<<asAz<<std::endl;
	std::cout<<"\tAD="<<discriminateAngle<<std::endl;
	std::cout<<"******************** "<<std::endl;
	sceneManager->setMinEpsilon(1e-4f);
	
	float currentElevation=initElevation;
	float currentAzimuth=initAzimuth;
	//std::cout<<"Tracing angle (el/az)="<<(currentElevation-overlap)<<","<<(currentElevation+deltaEl+overlap)<<"/"<<(currentAzimuth-overlap)<<","<<(currentAzimuth+deltaAz+overlap)<<std::endl;
	
	sceneManager->createRaySphere2D(currentElevation-overlap,asEl,currentElevation+deltaEl+overlap,currentAzimuth-overlap,asAz,currentAzimuth+deltaAz+overlap);
	
	//sceneManager->setPrintEnabled(1024*1024*1024, make_uint3(976,552,0));
	sceneManager->finishSceneContext();
	timer.start();
	
	//Angle for separating duplicate rays

	//LPCurvedMeshReflectionSimulation* sim=dynamic_cast<LPCurvedMeshReflectionSimulation*>(sceneManager->getSimulation());
	sim->setMaxAngleForDuplicateRays(discriminateAngle*M_PIf/180.f);
	
	
	
	
	uint nrx=10;	
	for (int i=1;i<=nrx;++i) {
		sceneManager->addReceiver(i,rx,polarizationRx, sphereRadius, sceneManager->printPower);
	}
	
	float zinit=50.0f;
	uint launches=0;

	for (int i=0;i<floor(500/nrx);++i) {
		//Update z
		for (int j=1;j<=nrx;++j) {
			rx.z=zinit;
			if (increaseRadius) {
				sceneManager->updateReceiver(j, rx,sphereRadius);
				sphereRadius += sphereDelta;
			} else {
				sceneManager->updateReceiver(j, rx);
			}
			zinit += 2.0f;
		}
		//First launch
		sceneManager->transmit(0, 1.0f, tx, polarizationTx, true);

		//Now loop to fill the solid angle
		currentAzimuth += deltaAz;
		//Trace all elevations
		while (currentElevation<endElevation) {

			//Trace all azimuth	
			while(currentAzimuth<endAzimuth) {
				//std::cout<<"Tracing angle (el/az)="<<(currentElevation-overlap)<<","<<(currentElevation+deltaEl+overlap)<<"/"<<(currentAzimuth-overlap)<<","<<(currentAzimuth+deltaAz+overlap)<<std::endl;
				sceneManager->createRaySphere2D(currentElevation-overlap,asEl,currentElevation+deltaEl+overlap,currentAzimuth-overlap,asAz,currentAzimuth+deltaAz+overlap);
				sceneManager->transmit(0, 1.0f, tx, polarizationTx, true);
				currentAzimuth += deltaAz;
			}
			currentAzimuth=initAzimuth;
			currentElevation += deltaEl;
		}
		sceneManager->endPartialLaunch(1u);

		currentElevation=initElevation;
		currentAzimuth=initAzimuth;
		++launches;

	} 
	
	timer.stop();
	std::cout<<"Time="<<timer.getTime()<<". Time/launch="<<(timer.getTime()/launches)<<std::endl;
	
}
void RouxAPTests::runTests(std::string test) {
	std::vector<float3> postx(4);
	std::vector<float3> posrx(4);
	std::vector<float3> pol(2);
	std::vector<float> freq(5);
	std::vector<float> disc(3);
	std::vector<float> rad(3);

	posrx[0]=postx[0]=make_float3(0.0f,0.0f, 0.0f); //C
	posrx[1]=postx[1]=make_float3(2.0f,0.0f, 0.0f); //NC
	posrx[2]=postx[2]=make_float3(3.0,3.0, 0.0f); //CTW. Is at 45 degrees 
	posrx[3]=postx[3]=make_float3(3.1,0.0, 0.0f); //VNC. 

	pol[0]=V;
	pol[1]=H;

	freq[0]=450e6;
	freq[1]=510e6;
	freq[2]=900e6;
	freq[3]=2e9;
	freq[4]=5e9;

	disc[0]=2.5f;
	disc[1]=15.0f;
	disc[2]=60.0f;
	
	rad[0]=1e-5;
	rad[1]=1e-7;
	rad[2]=1e-8;
	

	sphereDelta=sphereRadius/100;
	//Parse tests
	std::vector<int> tokens=parseTestString(test);

	frequency=freq[tokens[0]];
	polarizationTx=pol[tokens[1]];
	polarizationRx=pol[tokens[1]];
	tx=postx[tokens[2]];
	rx=posrx[tokens[3]];
	discriminateAngle=disc[tokens[4]];
	sphereDistanceError = rad[tokens[6]];
	//sphereDistanceError = 1e-8;
	//token[4] use also for filtering in RDN... change it
	filtering = static_cast<unsigned int>(tokens[4]);
	if (tokens[5]==0) {
		//runRDNCube();
		runRDN();
		//runRDNFixedRay();
	} else if (tokens[5]==2) {
		runRDNIsotropic();
		//runRDNCube();
	} else {
		//runAD();
		//runADCube();
		runADSingleRay();
	}

}
