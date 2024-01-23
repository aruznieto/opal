/***************************************************************/
//
//Copyright (c) 2021 Esteban Egea-Lopez http://girtel.upct.es/~simulations
//
/**************************************************************/
#include "anversDoa.h"
#include "../timer.h"
#include "tests.h"
#include "../flatSimulation.h"
#include "../ext/flatSimulationDoa.h"
#include "../ext/curvedFlatMeshSimulationDoa.h"
#include "../curvedMeshSimulation.h"
#include "../curvedFlatMeshSimulation.h"
#include "../rayDensityNormalizationSimulation.h"
#include "../singleDiffraction.h"
#include "../util.h"

AnversDoaTests::AnversDoaTests(OpalSceneManager* sceneManager, float sphereRadius, bool useDepolarization) : AnversTests(sceneManager,sphereRadius,useDepolarization) {
	this->emptyTunnel=false;
}

void AnversDoaTests::runTests(std::string test) {
	std::vector<float3> postx(4);
	std::vector<float3> posrx(2);
	std::vector<float3> pol(2);
	std::vector<float> freq(10);
	std::vector<float> disc(3);

	//Origin seem to be at NW corner of the tunnel cross section
	postx[0]=make_float3(0.50f,2.2f-6.0f, 0.0f);
	//postx[0]=make_float3(0.50f,2.52f-6.0f, 0.0f);
	postx[1]=make_float3(0.50f,2.18f-6.0f, 0.0f);
	postx[2]=make_float3(0.50f,1.84f-6.0f, 0.0f);
	postx[3]=make_float3(0.50f,1.50f-6.0f, 0.0f);


	posrx[0]=make_float3(3.2f,2.5f-6.0f, 200.0f);
	posrx[1]=make_float3(3.2f,2.5f-6.0f, 800.0f);
	
	pol[0]=V;
	pol[1]=H;
	float f=1.31e9;
	for (int i=0; i<10;++i) {
		freq[i]=f;
		f +=8e6;
	}
	disc[0]=2.5f;
	disc[1]=15.0f;
	disc[2]=30.0f;
	



	//Parse tests
	std::vector<int> tokens=parseTestString(test);
      
	frequency=freq[tokens[0]];
//Polarization
	polarizationTx=pol[tokens[1]];
	polarizationRx=pol[tokens[1]];
	if (tokens[1]==2) {
		polarizationTx=V;
		polarizationRx=H;
	}
	if (tokens[1]==3) {
		polarizationTx=H;
		polarizationRx=V;
	}
	tx=postx[tokens[2]];
	//Uset token[3] for empty tunnel
	if (tokens[3]==0) {
		emptyTunnel=true;
	} else {
		emptyTunnel=false;
	}
	rx=posrx[0];
	discriminateAngle=disc[tokens[4]];
	//token[4] use also for filtering in RDN... change it
	filtering = tokens[4];
	if (tokens[6]==0) {
		tunnelAsFlat=false;
	} else {
		tunnelAsFlat=true;
	}
	//std::cout<<"Running Anvers July Empty Curved Anvers tunnel with f="<<frequency<<"; tx ="<<tx<<"; rx="<<rx<<"; polTx="<<polarizationTx<<"; polRx="<<polarizationRx<<std::endl;
	if (tokens[5]==0) {
		runRDN();
	} else if (tokens[5]==1){
		runRDNIsotropic(true, false);
	} else if (tokens[5]==2) {
		runIsotropicTunnel(true);
	} else {
		runIsotropicTunnel(false);
	}
	
}
void AnversDoaTests::runDoaFlatTrucks(bool half, bool computeField, bool forward,  bool sectorized, bool useReflection, bool useDiffraction) { 
	
	OpalSimulation* sim;
	if (useRDN) {
		sim = new RayDensityNormalizationSimulation(sceneManager);
		dynamic_cast<RayDensityNormalizationSimulation*>(sim)->setExecutionMethod(RDNExecutionMode::HITINFO);	
	} else {	
		//sim= new LPFlatMeshReflectionSimulation(sceneManager);
		sim= new DoALPFlatMeshReflectionSimulation(sceneManager);
	}
        
	//sim->setExecutionMethod(RDNExecutionMode::HITINFO);
	//sim->setEnableTraceLog(true);
	sim->setPrintHits(true);
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
	simd->setPrintHits(true);
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
	//std::cout<<"**** Simulating Flat Anvers tunnel with trucks with isotropic RDN and diffraction ***** "<<std::endl;
	//std::cout<<"**** Simulating Flat Anvers tunnel with single middle truck with isotropic RDN and diffraction ***** "<<std::endl;
	//std::cout<<"**** Simulating Flat Anvers empty tunnel with same rx positions and  single iddle truck with isotropic RDN  ***** "<<std::endl;
	//std::cout<<"**** Simulating Flat Anvers tunnel with trucks and cars on free-lane with isotropic RDN and diffraction ***** "<<std::endl;
	if (emptyTunnel ) {
		std::cout<<"**** Simulating Flat Anvers tunnel empty and cars on free-lane ***** "<<std::endl;
	} else {
		std::cout<<"**** Simulating Flat Anvers tunnel with trucks and cars on free-lane  ***** "<<std::endl;
	}
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
	if (!emptyTunnel) {
        	std::string path("straight-trucks-metal.json");
        //std::string path("straight-single-truck.json");
		loadAnversJsonScenario(path);
	}
	//loadStraightAnversTunnel();	

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
    	int rayD=10000;	
	//With ray generation on launch
	sceneManager->finishSceneContext();

	//With LPFlatMeshReflectionSimulation
	if (useRDN) {
		if (!sectorized) {
			if (half) { 
				std::cout <<"**** Anvers Tunnel with  Half Sphere RDN ***"<<std::endl;	
				sceneManager->setRayRange(0.0,180.0,-90.0,90.0,rayD,rayD);
			} else {
				std::cout <<"**** Anvers Tunnel with Isotropic RDN ***"<<std::endl;	
				sceneManager->setRayRange(0.0,180.0,0.0,360.0,rayD,rayD);
			}

			//sceneManager->createRaySphereFromExternalBuffer(rayD,rayD,gen->getDevicePointer());
			//RayDensityNormalizationSimulation* sim=dynamic_cast<RayDensityNormalizationSimulation*>(sceneManager->getSimulation());
			if (half) { 
				dynamic_cast<RayDensityNormalizationSimulation*>(sim)->setInitialDensity(((float)sceneManager->getRaySphere().rayCount)/(2*M_PIf));
			} else {
				dynamic_cast<RayDensityNormalizationSimulation*>(sim)->setInitialDensity(((float)sceneManager->getRaySphere().rayCount)/(4*M_PIf));
			}
			dynamic_cast<RayDensityNormalizationSimulation*>(sim)->setFiltering(filtering);
		}
	} else {
		sceneManager->createRaySphere2D(0,0.1,180,-90,0.1,90);
	}

	
	timer.start();
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
	
	
//		sceneManager->addReceiver(1,make_float3(3.2,(2.7-6.0), 150),polarizationRx,0.274, sceneManager->printPower);
//		if (this->useAntennaGain) {
//			sceneManager->registerReceiverGain(1,gainId);
//		}
//		sceneManager->transmit(0, 1.0f, tx, polarizationTx, false);
//	return;	
	//uint nrx=380;	
	
	//Free lane
	float zpos[] = {300.0f};	
	float sphe[] = {0.574f};	
	//float zpos[] = {300.0f, 498.0f, 524.0f, 550.0f, 576.0f};	
	//float sphe[] = {0.574f, 0.97f, 1.02f, 1.074f, 1.126f};	
	for (int i=0;i<1;++i) {
	//for (int i=0;i<5;++i) {
		//sceneManager->addReceiver(i,make_float3(3.72f,2.7f-6.0f, 20.0f),polarization, sphereRadius, sceneManager->printPower);
		//Cars on Free lane
		float3 posrx=make_float3(3.2f,(2.7f-6.0f),zpos[i] );
		sceneManager->addReceiver(i,posrx,polarizationRx, sphe[i], sceneManager->printPower);
		if (this->useAntennaGain && !emptyTunnel) {
			sceneManager->registerReceiverGain(i,gainId);
		}
	}
	float tl=0;
	int launches=1;	
	if (sectorized) {
		SphereScanConfiguration c;
		//throw  opal::Exception("AnversTests::runAnversFlatStaticBetweenTrucksRDNIsotropic(): not implemented sectorized for LPFlatMeshReflectionSimulation yet");
		tl +=runSectorizedLaunch(c, 6, 1.0, tx, polarizationTx, sim);
	} else {
		//First launch
		//sceneManager->transmit(0, 1.0f, postx, polarizationTx, false);
		sceneManager->transmit(6, 1.0f, tx, polarizationTx, false);

	}

	if (sectorized) {
		std::cout<<"Total time="<<tl<<". Time/launch="<<(tl/launches)<<std::endl;
	} else {	
		timer.stop();
		std::cout<<"Time="<<timer.getTime()<<". Time/launch="<<(timer.getTime()/launches)<<std::endl;
	}
	
}
void AnversDoaTests::runTests6GHz(std::string test, bool useGain) {
	this->useAntennaGain=useGain;
	this->useRDN=false;
	std::vector<float3> postx(4);
	std::vector<float3> posrx(1);
	std::vector<float3> pol(4);
	std::vector<float> freq(10);

	//Origin seem to be at NW corner of the tunnel cross section
	postx[0]=make_float3(0.50f,2.2f-6.0f, 0.0f);
	//postx[0]=make_float3(0.50f,2.52f-6.0f, 0.0f);
	postx[1]=make_float3(0.50f,2.18f-6.0f, 0.0f);
	postx[2]=make_float3(0.50f,1.84f-6.0f, 0.0f);
	postx[3]=make_float3(0.50f,1.50f-6.0f, 0.0f);


	posrx[0]=make_float3(3.2f,2.5f-6.0f, 20.0f);
	
	const float3 VR=normalize(make_float3(1.0,1.0,0.0));
	const float3 VL=normalize(make_float3(-1.0,1.0,0.0));
	pol[0]=VR;
	pol[1]=VL;
	pol[2]=V;
	pol[3]=H;
	float f=5.9e9;
	for (int i=0; i<10;++i) {
		freq[i]=f;
		f +=8e6;
	}
	//Parse tests
	std::vector<int> tokens=parseTestString(test);
      
	frequency=freq[tokens[0]];
//Polarization
	polarizationTx=pol[tokens[1]];
	polarizationRx=pol[tokens[7]]; //New token...
	tx=postx[tokens[2]];
	rx=posrx[tokens[3]];
	//token[4] use also for filtering in RDN... change it
	//We usually keep the same filtering so I reuse it for empty

	tunnelAsFlat=false;
	bool field=false;
	if (tokens[6]==0) {
		field=false;
	} else {
		field=true;
	}
	bool forward=false;
//Forward means that the antenna radiation pattern is only going to take rays coming in the forward semisphere, that is, azimuth [-90,90]
//Since Z axis points along the tunnel and the transmitter is at the beginning, forward means that rays come back reflected from the end 
//of the tunnel. Backward is the opposite
//It would have probably been better calling 'forward' to the rays coming from the transmitter ('forward rays') and 'backward' those coming from the end
//But we did it the other way


	if (tokens[8]==1) {
		forward=true;
	}
	bool sectorized=false;
	if (tokens.size()>9 && tokens[9]==1) {
		sectorized=true;
		sphereDelta=2e-3;
	}
	bool useReflection=true;
	bool useDiffraction=true;

	if (tokens.size()>10 ) {
		if (tokens[10]==1) {
			//Only reflection
			useDiffraction=false;
		} else if (tokens[10]==2) {
			//Only diffraction
			useReflection=false;
		}
	}
		//sphereDelta=2e-4;
	if (tokens[4]==0) {
		emptyTunnel=false;
	} else {
		emptyTunnel=true;
	}
	if (tokens[5]==0) {
		//Flat tunnel with trucks or empty
		std::cout<<"Running Flat Anvers at 5.9 GHz /runDoaFlatTrucks()/ with f="<<frequency<<"; tx ="<<tx<<"; rx="<<rx<<"; polTx="<<polarizationTx<<"; polRx="<<polarizationRx<<"field="<<field<<"forward="<<forward<<"sectorized="<<sectorized<<"useReflection="<<useReflection<<"useDiffraction="<<useDiffraction<<"emptyTunnel="<<emptyTunnel<<std::endl;
		runDoaFlatTrucks(true,field, forward, sectorized, useReflection, useDiffraction);
	}  	
}
void AnversDoaTests::runIsotropicTunnel(bool curvedSim) {
	if (curvedSim) {
		DoALPCurvedFlatMeshReflectionSimulation* sim= new DoALPCurvedFlatMeshReflectionSimulation(sceneManager);
		/******** With angle discrimination ********/
		sim->setMaxAngleForDuplicateRays(discriminateAngle*M_PIf/180.f);
		sceneManager->setSimulation(sim);	
		std::cout<<"**** Simulating Anvers tunnel with CURVEDFLATWALLS ***** "<<std::endl;
	} else {
		DoALPFlatMeshReflectionSimulation* sim= new DoALPFlatMeshReflectionSimulation(sceneManager);
		//sim->setEnableTraceLog(true);
		sceneManager->setSimulation(sim);	
		std::cout<<"**** Simulating Anvers tunnel with FLAT ***** "<<std::endl;
	}
	sceneManager->enableGenerateRaysOnLaunch();	
	sceneManager->getSimulation()->setPrintHits(true);

	if (emptyTunnel==false) {
		//Add diffraction simulation
		SingleDiffraction* simd= new SingleDiffraction(sceneManager);
		sceneManager->setSimulation(simd);
		simd->setEnableSimulation(true);
		simd->setPrintHits(true);
		//simd->setEnableTraceLog(true);
	}

	sceneManager->initContext(frequency);

	if (emptyTunnel) {
		loadAnversTunnel();
	} else {
		
		if (curvedSim) {
			std::string path("trucks-metal.json");
			loadAnversJsonScenario(path);
		} else {
			std::string path("trucks-metal-flat.json");
			loadAnversJsonScenario(path);
		}
	}	


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
	float overlap=0.0f;
	//float overlap=0.0f;
	std::cout<<"Transmitting at "<<tx<<" with polarization="<<polarizationTx<<std::endl;
	std::cout<<"Receiving with radius="<<sphereRadius<<" with polarization="<<polarizationRx<<std::endl;
	std::cout<<"Scanning the sphere with ASElevation="<<asEl<< " and ASAzimuth="<<asAz<<std::endl;
	sceneManager->setMinEpsilon(1e-4f);
	
	float currentElevation=initElevation;
	float currentAzimuth=initAzimuth;
	//std::cout<<"Tracing angle (el/az)="<<(currentElevation-overlap)<<","<<(currentElevation+deltaEl+overlap)<<"/"<<(currentAzimuth-overlap)<<","<<(currentAzimuth+deltaAz+overlap)<<std::endl;
	sceneManager->finishSceneContext();
	
	sceneManager->createRaySphere2D(currentElevation,asEl,endElevation,currentAzimuth,asAz,endAzimuth);
	RaySphere sp=sceneManager->getRaySphere();
	std::cout<<"Ray density = "<<(sp.rayCount/(2.0*M_PIf))<<std::endl;	
	timer.start();
	
	//Angle for separating duplicate rays

	
	//float zinit=rx.z;
	//
	//float3 posrx;
	//
	//uint launches=0;
	////To get the correct z
	//float R=2900.0f;
	//float l=112.0f; //2*R*sin a/2=2*R*0.04f;
	//float h=2.32092875f; //h=R-0.5*sqrt(4*R*R-4*l*l)

	//if (zinit<=426) {
	//	posrx=make_float3(rx.x,(2.7f-6.0f)-(0.04f*zinit),zinit );
	//} else if ((zinit>426) && (zinit<=658)) {
	//	float zp=zinit-538; //z-426+112
	//	float yp=sqrt((R*R)-(zp*zp))-R+h; 
	//	float fl=-(23.04+yp); //Floor Y coordinate on curved section
	//	posrx=make_float3(rx.x, fl+2.7f,zinit );
	//} else {

	//	posrx=make_float3(rx.x,(2.7f-23.04f)+(0.04f*(zinit-658)),zinit );
	//}
	//sceneManager->addReceiver(1,posrx,polarizationRx, sphereRadius, sceneManager->printPower);
	

//Receivers, z positions: end of truck2, mid, beginning truck3
	//float zpos[] = {178.0f, 198.0f, 218.0f, 790.0f, 810.0f, 830.0f};	
	float zpos[] = {166.0f, 186.0f, 206.0f, 790.0f, 810.0f, 830.0f};	
//	float zinit=178.0f;  
	
	float3 posrx;
	
	uint launches=0;
	//To get the correct z
	float R=2900.0f;
	float l=112.0f; //2*R*sin a/2=2*R*0.04f;
	float h=2.32092875f; //h=R-0.5*sqrt(4*R*R-4*l*l)
	for (int i=0; i<6; i++) { 
		float zinit=zpos[i];
		if (zinit<=426) {
			posrx=make_float3(rx.x,(2.7f-6.0f)-(0.04f*zinit),zinit );
			sphereRadius=0.05;
		} else if ((zinit>426) && (zinit<=658)) {
			float zp=zinit-538; //z-426+112
			float yp=sqrt((R*R)-(zp*zp))-R+h; 
			float fl=-(23.04+yp); //Floor Y coordinate on curved section
			posrx=make_float3(rx.x, fl+2.7f,zinit );
		} else {

			sphereRadius=0.1;
			posrx=make_float3(rx.x,(2.7f-23.04f)+(0.04f*(zinit-658)),zinit );
		}
		sceneManager->addReceiver(i+1,posrx,polarizationRx, sphereRadius, sceneManager->printPower);
	}
		//First launch
	sceneManager->transmit(0, 1.0f, tx, polarizationTx, false);

	++launches;

	
	
	timer.stop();
	std::cout<<"Time="<<timer.getTime()<<". Time/launch="<<(timer.getTime()/launches)<<std::endl;
	
}
void AnversDoaTests::runTunnel() {
	LPCurvedFlatMeshReflectionSimulation* sim= new LPCurvedFlatMeshReflectionSimulation(sceneManager);
	sceneManager->setSimulation(sim);	
	sceneManager->enableGenerateRaysOnLaunch();	
	sceneManager->initContext(frequency);
	//sceneManager->getSimulation()->setPrintHits(false);


	loadAnversTunnel();

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
	float overlap=0.0f;
	//float overlap=0.0f;
	std::cout<<"**** Simulating Anvers tunnel with CURVEDFLATWALLS ***** "<<std::endl;
	std::cout<<"Transmitting at "<<tx<<" with polarization="<<polarizationTx<<std::endl;
	std::cout<<"Receiving with radius="<<sphereRadius<<" with polarization="<<polarizationRx<<std::endl;
	std::cout<<"Scanning the sphere with ASElevation="<<asEl<< " and ASAzimuth="<<asAz<<std::endl;
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
	
	
	
	/******** With angle discrimination ********/
	int runs=5;	
        int nrx =floor(1080/runs);
	for (int i=1;i<=nrx;++i) {
		//sceneManager->addReceiver(i,make_float3(3.72f,2.7f-6.0f, 20.0f),polarization, sphereRadius, sceneManager->printPower);
		sceneManager->addReceiver(i,rx,polarizationRx, sphereRadius, sceneManager->printPower);
	}
	
	float zinit=20.0f;
	uint launches=0;
	//To get the correct z
	float R=2900.0f;
	float l=112.0f; //2*R*sin a/2=2*R*0.04f;
	float h=2.32092875f; //h=R-0.5*sqrt(4*R*R-4*l*l)

	for (int i=0;i<runs;++i) {
		float3 posrx;
		for (int j=1;j<=nrx;++j) {

			if (zinit<=426) {
				posrx=make_float3(rx.x,(2.7f-6.0f)-(0.04f*zinit),zinit );
			} else if ((zinit>426) && (zinit<=658)) {
				float zp=zinit-538; //z-426+112
				float yp=sqrt((R*R)-(zp*zp))-R+h; 
				float fl=-(23.04+yp); //Floor Y coordinate on curved section
				posrx=make_float3(rx.x, fl+2.7f,zinit );
			} else {

				posrx=make_float3(rx.x,(2.7f-23.04f)+(0.04f*(zinit-658)),zinit );
			}
			if (increaseRadius) {
				sceneManager->updateReceiver(j, posrx,sphereRadius);
				sphereRadius += sphereDelta;
			} else {
				sceneManager->updateReceiver(j, posrx);
			}
			zinit += zStep;
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
void AnversDoaTests::runRDNIsotropic(bool half, bool computeField) {
	//float freq = 1.31e9f;
	RayDensityNormalizationSimulation* sim = new RayDensityNormalizationSimulation(sceneManager);
	sim->setExecutionMethod(RDNExecutionMode::HITINFO);	
	sceneManager->setSimulation(sim);
	ComputeMode mode=ComputeMode::VOLTAGE;
	if (computeField) {
		mode=ComputeMode::FIELD;
	}
	sim->setComputeMode(mode);
	sim->setPrintHits(true);
	sceneManager->enableGenerateRaysOnLaunch();	
	if (emptyTunnel==false) {
		
		//Add diffraction simulation
		SingleDiffraction* simd= new SingleDiffraction(sceneManager);
		sceneManager->setSimulation(simd);
		simd->setComputeMode(mode);
		simd->setEnableSimulation(true);
		simd->setPrintHits(true);
	}

	sceneManager->initContext(frequency);

	if (emptyTunnel) {
		loadAnversTunnel();
	} else {
        	std::string path("trucks-metal.json");
		loadAnversJsonScenario(path);
	}
	//loadStraightAnversTunnel(sceneManager);
	Timer timer;
	std::cout<<"**** Simulating Anvers DOA tunnel with isotropic RDN ***** "<<std::endl;
	std::cout<<"\tf="<<(frequency/1e6)<< " MHz; tx=" <<tx<<"; polarization="<<polarizationTx<<std::endl;
	std::cout<<"\trx="<<rx<<"; radius="<<sphereRadius;
	if (increaseRadius) {
		std::cout<<"+ "<<sphereDelta;
	}
	std::cout<<"; polarization="<<polarizationRx<<std::endl;
	sceneManager->setMinEpsilon(1e-4f);
	
	int rayD=10000;
	
	//With ray generation on launch
	sceneManager->finishSceneContext();
	if (half) { 
		std::cout <<"**** Anvers Tunnel with  Half Sphere RDN ***"<<std::endl;	
		sceneManager->setRayRange(0.0,180.0,-90.0,90.0,rayD,rayD);
	} else {
		std::cout <<"**** Anvers Tunnel with Isotropic RDN ***"<<std::endl;	
		sceneManager->setRayRange(0.0,180.0,0.0,360.0,rayD,rayD);
	}
	//sceneManager->createRaySphereFromExternalBuffer(rayD,rayD,gen->getDevicePointer());
	//RayDensityNormalizationSimulation* sim=dynamic_cast<RayDensityNormalizationSimulation*>(sceneManager->getSimulation());
	if (half) { 
		sim->setInitialDensity(((float)sceneManager->getRaySphere().rayCount)/(2*M_PIf));
	} else {
		sim->setInitialDensity(((float)sceneManager->getRaySphere().rayCount)/(4*M_PIf));
	}
	sim->setFiltering(filtering);
	
	//sceneManager->setPrintEnabled(1024*1024*1024, make_uint3(976,552,0));
	//sceneManager->finishSceneContext();
	timer.start();
	

	//Receivers, z positions: end of truck2, mid, beginning truck3
	float zpos[] = {166.0f, 186.0f, 206.0f, 790.0f, 810.0f, 830.0f};	
	//float zpos[] = {178.0f, 198.0f, 218.0f, 790.0f, 810.0f, 830.0f};	
//	float zinit=178.0f;  
	
	float3 posrx;
	
	uint launches=0;
	//To get the correct z
	float R=2900.0f;
	float l=112.0f; //2*R*sin a/2=2*R*0.04f;
	float h=2.32092875f; //h=R-0.5*sqrt(4*R*R-4*l*l)
	for (int i=0; i<6; i++) { 
		float zinit=zpos[i];
		if (zinit<=426) {
			sphereRadius=0.86f;
			posrx=make_float3(rx.x,(2.7f-6.0f)-(0.04f*zinit),zinit );
		} else if ((zinit>426) && (zinit<=658)) {
			float zp=zinit-538; //z-426+112
			float yp=sqrt((R*R)-(zp*zp))-R+h; 
			float fl=-(23.04+yp); //Floor Y coordinate on curved section
			posrx=make_float3(rx.x, fl+2.7f,zinit );
		} else {

			//Fixed here 
			sphereRadius=2.06f;
			posrx=make_float3(rx.x,(2.7f-23.04f)+(0.04f*(zinit-658)),zinit );
		}
		sceneManager->addReceiver(i+1,posrx,polarizationRx, sphereRadius, sceneManager->printPower);
	}
	sceneManager->transmit(0, 1.0f, tx, polarizationTx, false);

	++launches;

	
	
	timer.stop();
	std::cout<<"Time="<<timer.getTime()<<". Time/launch="<<(timer.getTime()/launches)<<std::endl;
	
}

