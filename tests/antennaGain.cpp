/***************************************************************/
//
//Copyright (c) 2021 Esteban Egea-Lopez http://girtel.upct.es/~simulations
//
/**************************************************************/

#include "antennaGain.h"
#include <limits>
#include "../singleDiffraction.h"
#include "../flatSimulation.h"
#include "../curvedFlatMeshSimulation.h"
#include "../basicSimulation.h"
#include "../rayDensityNormalizationSimulation.h"
#include "../util.h"
#include "../timer.h"

using namespace optix;
using namespace opal;
AntennaGainTests::AntennaGainTests(OpalSceneManager*   sceneManager, float sphereRadius) {
	this->sceneManager=sceneManager;
	this->sphereRadius=sphereRadius;
}
void AntennaGainTests::freeSpace(bool useDepolarization) {
	Timer timer;
	float freq = 5.9e9f;
	std::cout<<"Running free space test"<<std::endl;
        //sceneManager->enableGenerateRaysOnLaunch();	
	
	//Init context before doing anything else
	if (useDepolarization) {
		//LPCurvedFlatMeshReflectionSimulation* sim = new LPCurvedFlatMeshReflectionSimulation(sceneManager);
		LPFlatMeshReflectionSimulation* sim = new LPFlatMeshReflectionSimulation(sceneManager);
		sceneManager->setSimulation(sim);
	} else {
		BasicFlatMeshReflectionSimulation* sim = new BasicFlatMeshReflectionSimulation(sceneManager);
		sceneManager->setSimulation(sim);
	}
	
	//To see clearly the effects, it is necessary to uncomment the rtPrintf and/or printf functions in the optix/polarization files that
	//show the electric field vectors in the ray (ver_v, hor_v) and in the receivers, as well as the antenna gain

	// Consider our assumptions
	//- the orientation of the electric field vector on the ray when it is launched depends on the polarization of the antenna, since
	//  we are considering only linear polarizations and assume that the polarization is given by the antenna spatial orientation, which we actually call polarization,
	//  that is, we assume that a linear antenna (a dipole, for instance) oriented in Y axis (0,1,0) has a linear polarization=(0,1,0)
	//  If this antenna is rotated to (1,0,0), we assume the polarization is now (1,0,0) too.
	//- The antenna pattern file is assumed to be given with respect to a standard azimuth (Z is 0 degrees) and elevation (relative to Y). It may be independent of the polarization.
	// Imagine a directive antenna, let us say, with a beam pointing in the Z direction (0 azimuth). It may have linear polarization (0,1,0). If we orientate the antenna 
	//  90 degrees azimuth the beam is now pointing to X, but the linear polarization is still (0,1,0) it has not changed. We can tilt the beam -45 degrees to point to "ground" 
	//  and the polarization may be now (0,1,1) if we have physically rotated the antenna, or  still (0,1,0) if we have electrically rotated the beam 
	//  A dipole has a zero in the radiation pattern along the antenna orientation. If we rotate the antenna we rotate both the  polarization and the radiation pattern
	//  In Opal, we assume that if we change the polarization, we have physically rotated the antenna, and so the radiation pattern is also rotated accordingily,
	//  so when we use updateReceiver or transmit with a different polarization, the diagram is automatically rotated by Opal
	//  If one just want to rotate the radiation pattern without changing the polarization, one has to use 
	//


	sceneManager->setUseAntennaGain(true);
	sceneManager->enableMultiChannel();
	
	sceneManager->initContext(freq);
	sceneManager->setPrintRecords(true);
	timer.start();
	optix::float3 postx = make_float3(0.0f, 10.0f, 0.0f);
	optix::float3 polarizationTx = make_float3(0.0f, 1.0f, 0.0f); //Perpendicular to the floor. Assuming as in Unity that forward is z-axis and up is y-axis
	//optix::float3 polarizationTx = normalize(make_float3(1.0f, 1.0f, 0.0f)); 
	//optix::float3 polarization = normalize(make_float3(1.0f, 1.0f, 0.0f)); 
	optix::float3 polarization = normalize(make_float3(0.0f, 1.0f, 0.0f)); 
	optix::float3 posrx = make_float3(0.0f, 10.0f, 10.0f);
	sceneManager->addReceiver(1, posrx, polarization, sphereRadius, sceneManager->printPower);
	
	
	AntennaGain gains=sceneManager->loadGainsFromFileIndBPower("dipole.txt");
	int gainId=sceneManager->registerAntennaGain(gains);
	sceneManager->registerReceiverGain(1,gainId);
	sceneManager->registerTransmitterGain(0,gainId);
	//sceneManager->createRaySphere2DSubstep(1, 1); //0.1 degree delta step
	

	//***Single ray transmit****
	float3 mRay=normalize(make_float3(0.0,0,1));
	
	
	//float3 mRay=normalize(posrx-postx);
	sceneManager->createRaySphere2D(1,1,&mRay);
	
	sceneManager->finishSceneContext();
	//sceneManager->transmit(0, 1.0f, postx, polarizationTx, false);

	//Now test orientating the diagram but without changing the polarization
	Matrix<4,4> t=sceneManager->orientateAntennaPattern(90.0,-45.0);
	sceneManager->transmit(0, 1.0f, postx, polarizationTx, freq, false,t);

	timer.stop();
	std::cout<<"Time="<<timer.getTime()<<std::endl;

}
void AntennaGainTests::testAntennaOrientation() {
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
	//std::string gp("emf_5G/diagram.txt");
	std::string gp("dipole.txt");
	//std::string gainPathT = gp + gsufix + ".txt";
	AntennaGain gains = sceneManager->loadGainsFromFileIndBPower(gp.c_str());
	int gainIdTx = sceneManager->registerAntennaGain(gains);
	optix::float3 posrx = make_float3(10.0f, 10.0f, 0.0f);
	sceneManager->addReceiver(1, posrx, polarization, sphereRadius, sceneManager->printPower);
	sceneManager->registerTransmitterGain(0,gainIdTx);
	//***Single ray transmit****
	float3 mRay=normalize(make_float3(1,0,0));
	sceneManager->createRaySphere2D(1,1,&mRay);

	optix::float3 postx = make_float3(0.0f, 10.0f, 0.0f);
	
	sceneManager->finishSceneContext();
	//Matrix<4,4> t=sceneManager->orientateAntennaPattern(90.0,0.0);

	//Test passing point
	float3 mobile=make_float3(0,0,0);
	//float3 mobile=make_float3(-1,sqrt(2)*tan(M_PIf/4),-1);
	
	//The point position has be relative to the antenna source
	Matrix<4,4> t=sceneManager->pointAntennaPatternTo(postx, mobile);
	//Test passing the axis
	//float3 azimuthZero=make_float3(1,0,0);
	//float3 azimuthRotation = make_float3(0,1,0);
	//Matrix<4,4> t=sceneManager->orientateAntennaPattern(90.0,45.0,azimuthRotation, azimuthZero);

	sceneManager->getActiveTransmitters()[0]->transformToPolarization=t;
	sceneManager->transmit(0, 1.0f, postx, polarization, false,true);
	return;
}
void AntennaGainTests::freeSpaceRDN() {
	Timer timer;
	float freq = 5.9e9f;
	std::cout<<"Running free space test"<<std::endl;
	RayDensityNormalizationSimulation* sim= new RayDensityNormalizationSimulation(sceneManager);
	
	sim->setComputeMode(ComputeMode::VOLTAGE);
        sceneManager->enableGenerateRaysOnLaunch();	
	sceneManager->setSimulation(sim);
	sceneManager->setUseAntennaGain(true);
	
	sceneManager->initContext(freq);
	
	//sceneManager->getSimulation()->setPrintHits(true);	
	sceneManager->finishSceneContext();

	uint rayD=1000u;
	sceneManager->setRayRange(0.0,180.0,0.0,360.0,rayD,rayD);
	sim->setInitialDensity(((float)sceneManager->getRaySphere().rayCount)/(4*M_PIf));
	sim->setFiltering(2);

	timer.start();
	optix::float3 postx = make_float3(0.0f, 10.0f, 0.0f);
	optix::float3 polarizationTx = make_float3(0.0f, 1.0f, 0.0f); //Perpendicular to the floor. Assuming as in Unity that forward is z-axis and up is y-axis
	optix::float3 polarization = make_float3(0.0f, 1.0f, 0.0f); //Perpendicular to the floor. Assuming as in Unity that forward is z-axis and up is y-axis
	optix::float3 posrx = make_float3(0.0f, 10.0f, 10.0f);
	sceneManager->addReceiver(1, posrx, polarization, sphereRadius, sceneManager->printPower);
	AntennaGain gains=sceneManager->loadGainsFromFileIndBPower("gain17514.txt");
	int gainId=sceneManager->registerAntennaGain(gains);
	sceneManager->registerReceiverGain(1,gainId);
	sceneManager->registerTransmitterGain(0,gainId);
	
	sceneManager->transmit(0, 1.0f, postx, polarizationTx, false);
	timer.stop();
	std::cout<<"Time="<<timer.getTime()<<std::endl;
	//		postx = make_float3(-18.0f, 10.0f, 50.0f);
	//		sceneManager->transmit(0, 1.0f, postx, polarization);

}
