/***************************************************************/
//
//Copyright (c) 2019 Esteban Egea-Lopez http://girtel.upct.es/~simulations
//
/**************************************************************/
#include "anvers.h"
#include "tests.h"
#include "../timer.h"
#include <memory>
#include <fstream>
#include "../flatSimulation.h"
#include "../curvedMeshSimulation.h"
#include "../curvedFlatMeshSimulation.h"
#include "../rayDensityNormalizationSimulation.h"
#include "../singleDiffraction.h"
#include "../util.h"
using namespace opal;
using namespace optix;
AnversTests::AnversTests(OpalSceneManager* sceneManager, float sphereRadius, bool useDepolarization) : TunnelsBase(sceneManager,sphereRadius,useDepolarization) {
	this->rx=make_float3(0.0f,0.0f,0.0f);
	this->tx=make_float3(0.0f,0.0f,0.0f);
	this->polarizationTx = V;
	this->polarizationRx = V;
	this->discriminateAngle=2.5;
	this->useCurvedSection=true;
	this->increaseRadius =true;
	this->sphereDelta = 2e-3;
	this->tunnelAsFlat = false;
	this->useAntennaGain=false;
	this->filtering=2u;
	//this->zStep= 0.708f;
	this->zStep=1.0f;
	this->emptyTunnel =false;
	this->useRDN=true;
}
void AnversTests::loadAnversTunnel() {
	//Anvers is a slightly curved tunnel on the z direction with a rectangular cross section of 10.2 x 6 m (width,height) a length of around 1081 m and a maximum slope in Z of 4%
		
	MaterialEMProperties emProp1;
	emProp1.dielectricConstant = make_float2(5.0f, -60.0f*sceneManager->getChannelParameters().waveLength*0.01f);
	emProp1.tattenuation = make_float2(0.1f,-75.f );
	
	//All the meshes are already scaled, so the transform matrices should be identities
	//Origin is at top left corner, which means that the cross section goes from (0,0), (0,10.2), (-6,0), (-6,10.2)
	Matrix4x4 tm = Matrix4x4::identity();

	//Top wall is curved
	std::vector<int> hcind = sceneManager->loadTrianglesFromFile("meshes/anvers/top-i.txt");
	std::vector<float3> hcvert = sceneManager->loadVerticesFromFile("meshes/anvers/top-v.txt");
	std::vector<float4> pd1 = sceneManager->loadPDFromFile("meshes/anvers/top-pd1.txt");	
	std::vector<float4> pd2 = sceneManager->loadPDFromFile("meshes/anvers/top-pd2.txt");	
	std::cout << "Loading Anvers top with indices=" << hcind.size() << ", triangles="<<(hcind.size()/3)<<", vertices=" << hcvert.size() <<" and curvatures pd1="<<pd1.size()<<"pd2="<<pd2.size() << std::endl;
	if (tunnelAsFlat) {
		std::cout<<"Loading top as flat walls."<<std::endl;
		sceneManager->addStaticMesh(hcvert.size(), hcvert.data(), hcind.size(), hcind.data(), tm, emProp1 );
	} else {
		sceneManager->addStaticCurvedMesh(hcvert,  hcind, pd1, pd2, tm, emProp1, true);
	}
	//std::cout<<"Writing top wall" <<std::endl;
	//sceneManager->writeMeshToPLYFile("top.ply", hcvert,hcind, tm);

	//Bottom wall is curved
	 hcind = sceneManager->loadTrianglesFromFile("meshes/anvers/bottom-i.txt");
	 hcvert = sceneManager->loadVerticesFromFile("meshes/anvers/bottom-v.txt");
	 pd1 = sceneManager->loadPDFromFile("meshes/anvers/bottom-pd1.txt");	
	 pd2 = sceneManager->loadPDFromFile("meshes/anvers/bottom-pd2.txt");	
	std::cout << "Loading Anvers bottom with indices=" << hcind.size() << ", triangles="<<(hcind.size()/3)<<", vertices=" << hcvert.size() <<" and curvatures pd1="<<pd1.size()<<"pd2="<<pd2.size() << std::endl;
	if (tunnelAsFlat) {
		std::cout<<"Loading bottom as flat walls."<<std::endl;
		sceneManager->addStaticMesh(hcvert.size(), hcvert.data(), hcind.size(), hcind.data(), tm, emProp1 );
	} else {
		sceneManager->addStaticCurvedMesh(hcvert,  hcind, pd1, pd2, tm, emProp1, true);
	}
	//std::cout<<"Writing bottom wall" <<std::endl;
	//sceneManager->writeMeshToPLYFile("bottom.ply", hcvert,hcind, tm);
	 
	//left and right walls
	hcind = sceneManager->loadTrianglesFromFile("meshes/anvers/left-i.txt");
	 hcvert = sceneManager->loadVerticesFromFile("meshes/anvers/left-v.txt");
	std::cout << "Loading Anvers left with indices=" << hcind.size() << ", triangles="<<(hcind.size()/3)<<", vertices=" << hcvert.size() << std::endl;
	sceneManager->addStaticMesh(hcvert.size(), hcvert.data(), hcind.size(), hcind.data(), tm, emProp1 );
	
        //std::cout<<"Writing left wall" <<std::endl;
	//sceneManager->writeMeshToPLYFile("left.ply", hcvert,hcind, tm);
	
	hcind = sceneManager->loadTrianglesFromFile("meshes/anvers/right-i.txt");
	 hcvert = sceneManager->loadVerticesFromFile("meshes/anvers/right-v.txt");
	std::cout << "Loading Anvers right with indices=" << hcind.size() << ", triangles="<<(hcind.size()/3)<<", vertices=" << hcvert.size() << std::endl;
	sceneManager->addStaticMesh(hcvert.size(), hcvert.data(), hcind.size(), hcind.data(), tm, emProp1 );
        
	//std::cout<<"Writing right wall" <<std::endl;
	//sceneManager->writeMeshToPLYFile("right.ply", hcvert,hcind, tm);
}
void AnversTests::loadStraightAnversTunnel() {
	//Anvers is a slightly curved tunnel on the z direction with a rectangular cross section of 10.2 x 6 m (width,height) a length of around 1081 m and a maximum slope in Z of 4%
	//Here we load a tunnel with a down slope but with straight sections, which would make a peak in the middle
		
	MaterialEMProperties emProp1;
	emProp1.dielectricConstant = make_float2(5.0f, -60.0f*sceneManager->getChannelParameters().waveLength*0.01f);
	emProp1.tattenuation = make_float2(0.1f,-75.f );
	
	//All the meshes are already scaled, so the transform matrices should be identities
	//Origin is at top left corner, which means that the cross section goes from (0,0), (0,10.2), (-6,0), (-6,10.2)
	Matrix4x4 tm = Matrix4x4::identity();

	std::vector<int> hcind = sceneManager->loadTrianglesFromFile("meshes/anvers/top-straight-i.txt");
	std::vector<float3> hcvert = sceneManager->loadVerticesFromFile("meshes/anvers/top-straight-v.txt");
	std::cout << "Loading Anvers straight top with indices=" << hcind.size() << ", triangles="<<(hcind.size()/3)<<", vertices=" << hcvert.size() << std::endl;
	sceneManager->addStaticMesh(hcvert.size(), hcvert.data(), hcind.size(), hcind.data(), tm, emProp1 );

	//Bottom wall is curved
	 hcind = sceneManager->loadTrianglesFromFile("meshes/anvers/bottom-straight-i.txt");
	 hcvert = sceneManager->loadVerticesFromFile("meshes/anvers/bottom-straight-v.txt");
	std::cout << "Loading Anvers straight bottom with indices=" << hcind.size() << ", triangles="<<(hcind.size()/3)<<", vertices=" << hcvert.size() << std::endl;
	sceneManager->addStaticMesh(hcvert.size(), hcvert.data(), hcind.size(), hcind.data(), tm, emProp1 );
	 
	//left and right walls
	hcind = sceneManager->loadTrianglesFromFile("meshes/anvers/left-straight-i.txt");
	 hcvert = sceneManager->loadVerticesFromFile("meshes/anvers/left-straight-v.txt");
	std::cout << "Loading Anvers straight left with indices=" << hcind.size() << ", triangles="<<(hcind.size()/3)<<", vertices=" << hcvert.size() << std::endl;
	sceneManager->addStaticMesh(hcvert.size(), hcvert.data(), hcind.size(), hcind.data(), tm, emProp1 );
	
	hcind = sceneManager->loadTrianglesFromFile("meshes/anvers/right-straight-i.txt");
	 hcvert = sceneManager->loadVerticesFromFile("meshes/anvers/right-straight-v.txt");
	std::cout << "Loading Anvers straight right with indices=" << hcind.size() << ", triangles="<<(hcind.size()/3)<<", vertices=" << hcvert.size() << std::endl;
	sceneManager->addStaticMesh(hcvert.size(), hcvert.data(), hcind.size(), hcind.data(), tm, emProp1 );
}
void AnversTests::loadAnversScenario(std::string path) {
	ScenarioLoader* loader=new ScenarioLoader(sceneManager);
	loader->loadMeshesFromFiles(path);
	loader->loadEdgesFromFiles(path);
}
void AnversTests::loadAnversJsonScenario(std::string path) {
	ScenarioLoader* loader=new ScenarioLoader(sceneManager);
	std::cout<<"Path"<<path<<std::endl;
	loader->loadJSONScenario(path);
}
void AnversTests::runTests6GHz(std::string test, bool useGain) {
	this->useAntennaGain=useGain;
	//this->useRDN=false;
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
		std::cout<<"Running Flat Anvers at 5.9 GHz /runTestsBetweenLeanAntenna()/ with f="<<frequency<<"; tx ="<<tx<<"; rx="<<rx<<"; polTx="<<polarizationTx<<"; polRx="<<polarizationRx<<"field="<<field<<"forward="<<forward<<"sectorized="<<sectorized<<"useReflection="<<useReflection<<"useDiffraction="<<useDiffraction<<"emptyTunnel="<<emptyTunnel<<std::endl;
		runAnversFlatStaticBetweenTrucksRDNIsotropic(true,field, forward, sectorized, useReflection, useDiffraction);
		//runAnversFlatOvertakingTrucksRDNIsotropic(true,field, forward, sectorized, useReflection, useDiffraction);
	} else if (tokens[5]==1){
		//runRDNIsotropic(true, true);
		std::cout<<"Running Flat Anvers at 5.9 GHz overtaking /runTestsBetweenLeanAntenna()/ with f="<<frequency<<"; tx ="<<tx<<"; rx="<<rx<<"; polTx="<<polarizationTx<<"; polRx="<<polarizationRx<<"field="<<field<<"forward="<<forward<<"sectorized="<<sectorized<<"useReflection="<<useReflection<<"useDiffraction="<<useDiffraction<<"emptyTunnel="<<emptyTunnel<<std::endl;
		runAnversFlatOvertakingTrucksRDNIsotropic(true,field, forward, sectorized, useReflection, useDiffraction);
	} else {
		//Curved empty Anvers with LPFlatMeshReflectionSimulation 
		runTunnel();
	}
	
}
void AnversTests::runTestsBetweenLeanAntenna(std::string test,  bool useGain) {
	this->useAntennaGain=useGain;
	std::vector<float3> postx(4);
	std::vector<float3> posrx(1);
	std::vector<float3> pol(4);
	std::vector<float> freq(10);
	std::vector<float> disc(3);

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
      
	//runAnversTunnelNoCurvedSection(sceneManager,sphereRadius,useDepolarization,freq[tokens[0]],postx[tokens[2]],rx[tokens[3]],pol[tokens[1]]);
	//runAnversTunnel(sceneManager,sphereRadius,useDepolarization,freq[tokens[0]],postx[tokens[2]],rx[tokens[3]],pol[tokens[1]]);
	frequency=freq[tokens[0]];
//Polarization
	polarizationTx=pol[tokens[1]];
	polarizationRx=pol[tokens[7]]; //New token...
	tx=postx[tokens[2]];
	rx=posrx[tokens[3]];
	discriminateAngle=disc[tokens[4]];
	//token[4] use also for filtering in RDN... change it
	filtering = tokens[4];
	tunnelAsFlat=false;
	bool field=false;
	if (tokens[6]==0) {
		field=false;
	} else {
		field=true;
	}
	bool forward=false;
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
	std::cout<<"Running Anvers Lean Trucks /runTestsBetweenLeanAntenna()/ Curved Anvers tunnel with f="<<frequency<<"; tx ="<<tx<<"; rx="<<rx<<"; polTx="<<polarizationTx<<"; polRx="<<polarizationRx<<"field="<<field<<"forward="<<forward<<"sectorized="<<sectorized<<"useReflection="<<useReflection<<"useDiffraction="<<useDiffraction<<std::endl;
	if (tokens[5]==0) {
		//runRDN();
		//runAnversStaticTruckRDNIsotropic(true,true);
		runAnversFlatStaticBetweenTrucksRDNIsotropic( true,field, forward, sectorized, useReflection, useDiffraction);
		//runAnversStaticBetweenTrucksRDNIsotropic(true,field, false);
	} else if (tokens[5]==1){
		runRDNIsotropic(true, true);
		//runRDNStraightTunnel(true);
		//runStraightTunnel();
	} else {
		runTunnel();
	}
	
}
void AnversTests::runTestsDepolarizationLeanAntenna(std::string test, bool useGain) {
	this->useAntennaGain=useGain;
	std::vector<float3> postx(4);
	std::vector<float3> posrx(1);
	std::vector<float3> pol(4);
	std::vector<float> freq(10);
	std::vector<float> disc(3);

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
      
	//runAnversTunnelNoCurvedSection(sceneManager,sphereRadius,useDepolarization,freq[tokens[0]],postx[tokens[2]],rx[tokens[3]],pol[tokens[1]]);
	//runAnversTunnel(sceneManager,sphereRadius,useDepolarization,freq[tokens[0]],postx[tokens[2]],rx[tokens[3]],pol[tokens[1]]);
	frequency=freq[tokens[0]];
//Polarization
	polarizationTx=pol[tokens[1]];
	polarizationRx=pol[tokens[7]]; //New token...
	tx=postx[tokens[2]];
	rx=posrx[tokens[3]];
	discriminateAngle=disc[tokens[4]];
	//token[4] use also for filtering in RDN... change it
	filtering = tokens[4];
	if (tokens[6]==0) {
		tunnelAsFlat=false;
	} else {
		tunnelAsFlat=true;
	}
	std::cout<<"Running Anvers Lean Empty /runTestsDepolartizationLeanAntenna()/ Curved Anvers tunnel with f="<<frequency<<"; tx ="<<tx<<"; rx="<<rx<<"; polTx="<<polarizationTx<<"; polRx="<<polarizationRx<<std::endl;
	if (tokens[5]==0) {
		//runRDN();
		runAnversStaticTruckRDNIsotropic(true,true);
	} else if (tokens[5]==1){
		runRDNIsotropic(true, true);
		//runRDNStraightTunnel(true);
		//runStraightTunnel();
	} else {
		runTunnel();
	}
	
}
void AnversTests::runTestsLeanAntenna(std::string test, bool useGain) {
	this->useAntennaGain=useGain;
	std::vector<float3> postx(4);
	std::vector<float3> posrx(1);
	std::vector<float3> pol(4);
	std::vector<float> freq(10);
	std::vector<float> disc(3);

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
      
	//runAnversTunnelNoCurvedSection(sceneManager,sphereRadius,useDepolarization,freq[tokens[0]],postx[tokens[2]],rx[tokens[3]],pol[tokens[1]]);
	//runAnversTunnel(sceneManager,sphereRadius,useDepolarization,freq[tokens[0]],postx[tokens[2]],rx[tokens[3]],pol[tokens[1]]);
	frequency=freq[tokens[0]];
//Polarization
	polarizationTx=pol[tokens[1]];
	polarizationRx=pol[tokens[7]]; //New token...
	tx=postx[tokens[2]];
	rx=posrx[tokens[3]];
	discriminateAngle=disc[tokens[4]];
	//token[4] use also for filtering in RDN... change it
	filtering = tokens[4];
	if (tokens[6]==0) {
		tunnelAsFlat=false;
	} else {
		tunnelAsFlat=true;
	}
	std::cout<<"Running Anvers Lean Empty Curved Anvers tunnel with f="<<frequency<<"; tx ="<<tx<<"; rx="<<rx<<"; polTx="<<polarizationTx<<"; polRx="<<polarizationRx<<std::endl;
	if (tokens[5]==0) {
		runRDN();
	} else if (tokens[5]==1){
		runRDNIsotropic(true, false);
		//runAnversStaticTruckRDNIsotropic(true);
		//runRDNStraightTunnel(true);
		//runStraightTunnel();
	} else {
		runTunnel();
	}
	
}
void AnversTests::runTestsJuly21(std::string test) {
	std::vector<float3> postx(4);
	std::vector<float3> posrx(1);
	std::vector<float3> pol(2);
	std::vector<float> freq(10);
	std::vector<float> disc(3);

	//Origin seem to be at NW corner of the tunnel cross section
	postx[0]=make_float3(0.50f,2.2f-6.0f, 0.0f);
	//postx[0]=make_float3(0.50f,2.52f-6.0f, 0.0f);
	postx[1]=make_float3(0.50f,2.18f-6.0f, 0.0f);
	postx[2]=make_float3(0.50f,1.84f-6.0f, 0.0f);
	postx[3]=make_float3(0.50f,1.50f-6.0f, 0.0f);


	posrx[0]=make_float3(3.2f,2.5f-6.0f, 20.0f);
	
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
      
	//runAnversTunnelNoCurvedSection(sceneManager,sphereRadius,useDepolarization,freq[tokens[0]],postx[tokens[2]],rx[tokens[3]],pol[tokens[1]]);
	//runAnversTunnel(sceneManager,sphereRadius,useDepolarization,freq[tokens[0]],postx[tokens[2]],rx[tokens[3]],pol[tokens[1]]);
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
	rx=posrx[tokens[3]];
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
		//runRDNIsotropic(true, false);
		runAnversStaticTruckRDNIsotropic(true, false);
		//runRDNStraightTunnel(true);
		//runStraightTunnel();
	} else {
		runTunnel();
	}
	
}
void AnversTests::runTests(std::string test) {
	std::vector<float3> postx(4);
	std::vector<float3> posrx(4);
	std::vector<float3> pol(2);
	std::vector<float> freq(10);
	std::vector<float> disc(3);

	//Origin seem to be at NW corner of the tunnel cross section
	postx[0]=make_float3(0.50f,2.52f-6.0f, 0.0f);
	postx[1]=make_float3(0.50f,2.18f-6.0f, 0.0f);
	postx[2]=make_float3(0.50f,1.84f-6.0f, 0.0f);
	postx[3]=make_float3(0.50f,1.50f-6.0f, 0.0f);
	posrx[0]=make_float3(3.72f,2.7f-6.0f, 20.0f);
	posrx[1]=make_float3(3.38f,2.7f-6.0f, 20.0f);
	posrx[2]=make_float3(3.04f,2.7f-6.0f, 20.0f);
	posrx[3]=make_float3(2.7f,2.7f-6.0f, 20.0f);
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
      
	//runAnversTunnelNoCurvedSection(sceneManager,sphereRadius,useDepolarization,freq[tokens[0]],postx[tokens[2]],rx[tokens[3]],pol[tokens[1]]);
	//runAnversTunnel(sceneManager,sphereRadius,useDepolarization,freq[tokens[0]],postx[tokens[2]],rx[tokens[3]],pol[tokens[1]]);
	frequency=freq[tokens[0]];
	polarizationTx=pol[tokens[1]];
	polarizationRx=pol[tokens[1]];
	tx=postx[tokens[2]];
	rx=posrx[tokens[3]];
	discriminateAngle=disc[tokens[4]];
	//token[4] use also for filtering in RDN... change it
	filtering = tokens[4];
	if (tokens[6]==0) {
		tunnelAsFlat=false;
	} else {
		tunnelAsFlat=true;
	}
	if (tokens[5]==0) {
		runRDN();
	} else if (tokens[5]==1){
		runRDNIsotropic(true, false);
	} else {
		runTunnel();
	}
	
}
void AnversTests::runRDNStraightTunnel(bool half) {
	RayDensityNormalizationSimulation* sim = new RayDensityNormalizationSimulation(sceneManager);
	sceneManager->setSimulation(sim);
	ComputeMode mode=ComputeMode::VOLTAGE;
	sim->setComputeMode(mode);
	sceneManager->enableGenerateRaysOnLaunch();
	frequency=5.9e9;
	polarizationRx=V;
	polarizationTx=V;

	sceneManager->initContext(frequency);
	//sceneManager->getSimulation()->setPrintHits(true);
	
	
	//****Cube tunnel dimensions *****
	float width=10.2f;
	float height=6.0f;
	float length=1100.0f;
	//****************
	
	MaterialEMProperties emProp1;
	emProp1.dielectricConstant = make_float2(5.0f, -60.0f*sceneManager->getChannelParameters().waveLength*0.01f);
	emProp1.tattenuation = make_float2(0.1f,-75.f );

	loadSquareTunnel( width, height, length, emProp1);
	
	
	Timer timer;
	std::cout<<"**** Simulating Anvers tunnel with isotropic RDN ***** "<<std::endl;
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
	sceneManager->setMinEpsilon(1e-4f);
	
	timer.start();
	uint nrx=360;	
	float zinit=20.0f;
	float3	posrx=make_float3(-width/2.0f + 3.2f,-height/2.0f + 2.5f,zinit );
	for (int i=1;i<=nrx;++i) {
		//sceneManager->addReceiver(i,make_float3(3.72f,2.7f-6.0f, 20.0f),polarization, sphereRadius, sceneManager->printPower);
		sceneManager->addReceiver(i,posrx,polarizationRx, sphereRadius, sceneManager->printPower);
	}
	
	
	uint launches=0;
	optix::float3 postx = make_float3(-width/2.0f + 0.50f,-height/2.0f +2.52f, 0.0f);
	ResultReport report;
	for (int i=0;i<floor(1080/nrx);++i) {
		for (int j=1;j<=nrx;++j) {


			posrx.z=zinit ;
			if (increaseRadius) {
				sceneManager->updateReceiver(j, posrx,sphereRadius);
				sphereRadius += sphereDelta;
			} else {
				sceneManager->updateReceiver(j, posrx);
			}
			zinit += zStep;
		}
		//First launch
		ResultReport* rp=sceneManager->transmit(0, 1.0f,postx, polarizationTx, false);
		report.merge(*rp);
		++launches;

	} 
	report.toCSV("rdnanv.csv");
	timer.stop();
	std::cout<<"Time="<<timer.getTime()<<". Time/launch="<<(timer.getTime()/launches)<<std::endl;
	
}
void AnversTests::runStraightTunnel() {
	//float freq = 1.31e9f;
	LPFlatMeshReflectionSimulation* sim= new LPFlatMeshReflectionSimulation(sceneManager);
	sceneManager->setSimulation(sim);	
	ComputeMode mode=ComputeMode::VOLTAGE;
	sim->setComputeMode(mode);
	sceneManager->enableGenerateRaysOnLaunch();	
	sceneManager->initContext(frequency);
	//sceneManager->getSimulation()->setPrintHits(false);
	
	
	//****Cube tunnel dimensions *****
	float width=10.2f;
	float height=6.0f;
	float length=1100.0f;
	//****************
	
	MaterialEMProperties emProp1;
	emProp1.dielectricConstant = make_float2(5.0f, -60.0f*sceneManager->getChannelParameters().waveLength*0.01f);
	emProp1.tattenuation = make_float2(0.1f,-75.f );

	loadSquareTunnel( width, height, length, emProp1);
	
	sceneManager->setMinEpsilon(1e-4f);
	sceneManager->finishSceneContext();
	
	Timer timer;
	//Receiver polarization
	optix::float3 polarization = V; 	
	
	//Transmitter
	optix::float3 postx = make_float3(-width/2.0f + 0.50f,-height/2.0f +2.52f, 0.0f);
	//float3 polarizationTx = V; 
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
	
	std::cout<<"**** Simulating Straight tunnel equivalent to Anvers ***** "<<std::endl;
	std::cout<<"Transmitting at "<<postx<<" with polarization="<<polarizationTx<<std::endl;
	std::cout<<"Receiving with radius="<<sphereRadius<<" with polarization="<<polarization<<std::endl;
	std::cout<<"Scanning the sphere with ASElevation="<<asEl<< " and ASAzimuth="<<asAz<<std::endl;
	
	
	
	
	float currentElevation=initElevation;
	float currentAzimuth=initAzimuth;
	//std::cout<<"Tracing angle (el/az)="<<(currentElevation-overlap)<<","<<(currentElevation+deltaEl+overlap)<<"/"<<(currentAzimuth-overlap)<<","<<(currentAzimuth+deltaAz+overlap)<<std::endl;
	//sceneManager->setPrintEnabled(1024*1024*1024, make_uint3(976,552,0));
	timer.start();
	
	sceneManager->createRaySphere2D(currentElevation-overlap,asEl,currentElevation+deltaEl+overlap,currentAzimuth-overlap,asAz,currentAzimuth+deltaAz+overlap);

	
	
	
	/******** With angle discrimination ********/
	/** Sequential test. For a really large number of reflections and receivers  (that is, potential hits) have to do a sequential test, to avoid buffer overflows */
//	sceneManager->addReceiver(1,make_float3(3.72f,2.7f-6.0f, 20.0f),polarization, sphereRadius, sceneManager->printPower);
//	sceneManager->addReceiver(2,make_float3(3.72f,2.7f-6.0f, 20.0f),polarization, sphereRadius, sceneManager->printPower);
//	sceneManager->addReceiver(3,make_float3(3.72f,2.7f-6.0f, 20.0f),polarization, sphereRadius, sceneManager->printPower);
//	sceneManager->addReceiver(4,make_float3(3.72f,2.7f-6.0f, 20.0f),polarization, sphereRadius, sceneManager->printPower);
	
	float zinit=20.0f;
	float3	posrx=make_float3(-width/2.0f + 3.2f,-height/2.0f + 2.5f,zinit );
	for (int i=1;i<=360;++i) {
		//sceneManager->addReceiver(i,make_float3(3.72f,2.7f-6.0f, 20.0f),polarization, sphereRadius, sceneManager->printPower);
		sceneManager->addReceiver(i,posrx,polarization, sphereRadius, sceneManager->printPower);
	}
	
	uint launches=0;
	
	for (int i=0;i<3;++i) {
		float3 posrx;
		for (int j=1;j<=360;++j) {


			posrx.z=zinit;
			if (increaseRadius) {
				sceneManager->updateReceiver(j, posrx,sphereRadius);
				sphereRadius += sphereDelta;
			} else {
				sceneManager->updateReceiver(j, posrx);
			}
			zinit += zStep;
		}
		//First launch
		sceneManager->transmit(0, 1.0f,postx, polarizationTx, true);

		//Now loop to fill the solid angle
		currentAzimuth += deltaAz;
		//Trace all elevations
		while (currentElevation<endElevation) {

			//Trace all azimuth	
			while(currentAzimuth<endAzimuth) {
				//std::cout<<"Tracing angle (el/az)="<<(currentElevation-overlap)<<","<<(currentElevation+deltaEl+overlap)<<"/"<<(currentAzimuth-overlap)<<","<<(currentAzimuth+deltaAz+overlap)<<std::endl;
				sceneManager->createRaySphere2D(currentElevation-overlap,asEl,currentElevation+deltaEl+overlap,currentAzimuth-overlap,asAz,currentAzimuth+deltaAz+overlap);
				sceneManager->transmit(0, 1.0f, postx, polarizationTx, true);
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
void AnversTests::runTunnel() {
	LPCurvedFlatMeshReflectionSimulation* sim= new LPCurvedFlatMeshReflectionSimulation(sceneManager);
	sceneManager->setSimulation(sim);	
	sceneManager->enableGenerateRaysOnLaunch();	
	sceneManager->initContext(frequency);
	sceneManager->getSimulation()->setPrintHits(false);


	loadAnversTunnel();

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
	
	sceneManager->finishSceneContext();
	sceneManager->createRaySphere2D(currentElevation-overlap,asEl,currentElevation+deltaEl+overlap,currentAzimuth-overlap,asAz,currentAzimuth+deltaAz+overlap);
	
	//sceneManager->setPrintEnabled(1024*1024*1024, make_uint3(976,552,0));
	timer.start();
	
	//Angle for separating duplicate rays

	//LPCurvedMeshReflectionSimulation* sim=dynamic_cast<LPCurvedMeshReflectionSimulation*>(sceneManager->getSimulation());
	sim->setMaxAngleForDuplicateRays(discriminateAngle*M_PIf/180.f);
	
	
	
	/******** With angle discrimination ********/
	/** Sequential test. For a really large number of reflections and receivers  (that is, potential hits) have to do a sequential test, to avoid buffer overflows */
//	sceneManager->addReceiver(1,make_float3(3.72f,2.7f-6.0f, 20.0f),polarization, sphereRadius, sceneManager->printPower);
//	sceneManager->addReceiver(2,make_float3(3.72f,2.7f-6.0f, 20.0f),polarization, sphereRadius, sceneManager->printPower);
//	sceneManager->addReceiver(3,make_float3(3.72f,2.7f-6.0f, 20.0f),polarization, sphereRadius, sceneManager->printPower);
//	sceneManager->addReceiver(4,make_float3(3.72f,2.7f-6.0f, 20.0f),polarization, sphereRadius, sceneManager->printPower);
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

	//Separate the 1500 positions in 3 runs (i) and put 500 receivers simultaneously in each run (with CURVEDFLATWALLS, Opal can manage that without memory overflow)	
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
void AnversTests::runNoCurvedSection() {
	LPFlatMeshReflectionSimulation* sim= new LPFlatMeshReflectionSimulation(sceneManager);
	sceneManager->setSimulation(sim);	
	sceneManager->initContext(frequency);
	sceneManager->getSimulation()->setPrintHits(false);


	loadStraightAnversTunnel();

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
	std::cout<<"**** Simulating Anvers tunnel No Curved Sections ***** "<<std::endl;
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
	
	
	
	
	/******** With angle discrimination ********/
	/** Sequential test. For a really large number of reflections and receivers  (that is, potential hits) have to do a sequential test, to avoid buffer overflows */
//	sceneManager->addReceiver(1,make_float3(3.72f,2.7f-6.0f, 20.0f),polarization, sphereRadius, sceneManager->printPower);
//	sceneManager->addReceiver(2,make_float3(3.72f,2.7f-6.0f, 20.0f),polarization, sphereRadius, sceneManager->printPower);
//	sceneManager->addReceiver(3,make_float3(3.72f,2.7f-6.0f, 20.0f),polarization, sphereRadius, sceneManager->printPower);
//	sceneManager->addReceiver(4,make_float3(3.72f,2.7f-6.0f, 20.0f),polarization, sphereRadius, sceneManager->printPower);
	
	for (int i=1;i<=500;++i) {
		//sceneManager->addReceiver(i,make_float3(3.72f,2.7f-6.0f, 20.0f),polarization, sphereRadius, sceneManager->printPower);
		sceneManager->addReceiver(i,rx,polarizationRx, sphereRadius, sceneManager->printPower);
	}
	
	float zinit=20.0f;
	uint launches=0;
	//To get the correct z
	
	for (int i=0;i<3;++i) {
		float3 posrx;
		for (int j=1;j<=500;++j) {

			if (zinit<=568) {
				posrx=make_float3(rx.x,(2.7f-6.0f)-(0.04f*zinit),zinit );
			} else {

				posrx=make_float3(rx.x,(2.7f-28.72f)+(0.04f*(zinit-568)),zinit );
			}
			if (increaseRadius) {
				sceneManager->updateReceiver(j, posrx,sphereRadius);
				sphereRadius += sphereDelta;
			} else {
				sceneManager->updateReceiver(j, posrx);
			}
			zinit += 0.75f;
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
void AnversTests::runRDN() {
	//float freq = 1.31e9f;
	RayDensityNormalizationSimulation* sim = new RayDensityNormalizationSimulation(sceneManager);
	sceneManager->setSimulation(sim);
	ComputeMode mode=ComputeMode::VOLTAGE;
	sim->setComputeMode(mode);
	sceneManager->enableGenerateRaysOnLaunch(); //No memory used to store the rays: more ray density
	sceneManager->initContext(frequency);
	//sceneManager->getSimulation()->setPrintHits(true);


	//loadStraightAnversTunnel();
	loadAnversTunnel();
	//****Cube tunnel dimensions *****
	//float width=10.2f;
	//float height=6.0f;
	//float length=1100.0f;
	////****************
	//MaterialEMProperties emProp1;
	//emProp1.dielectricConstant = make_float2(5.0f, -60.0f*sceneManager->getChannelParameters().waveLength*0.01f);
	//emProp1.tattenuation = make_float2(0.1f,-75.f );
	//loadSquareTunnel( width, height, length, emProp1);
	//tx = make_float3(-width/2.0f + 0.50f,-height/2.0f +2.52f, 0.0f);

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
	std::cout<<"**** Simulating Anvers tunnel with sectorized RDN ***** "<<std::endl;
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
	float deg2rad=M_PIf/180.0f;
	float rayGoal=1e9;
	float solidAngle=deg2rad*deltaAz*(cosf(deg2rad*currentElevation)-cosf(deg2rad*(currentElevation+deltaEl)));
	int rayD=floor(sqrt(rayGoal*solidAngle)); 
	std::cout<<"\t Scanning the sphere with ASElevation="<<asEl<< " and ASAzimuth="<<asAz<<"rays="<<rayD*rayD<<std::endl;
	//int rayD=10000;
	//OpalRaySphereGenerator* gen=sceneManager->getRaySphereGenerator();
	//gen->generateRandomUniformOnDevice(currentElevation,currentElevation+deltaEl,currentAzimuth,currentAzimuth+deltaAz,rayD*rayD);
	//sceneManager->createRaySphereFromExternalBuffer(rayD,rayD,gen->getDevicePointer());
	//RayDensityNormalizationSimulation* sim=dynamic_cast<RayDensityNormalizationSimulation*>(sceneManager->getSimulation());
	//sim->setInitialDensity(sceneManager->getRaySphere().rayCount,currentAzimuth,currentAzimuth+deltaAz,currentElevation, currentElevation+deltaEl);
	//sim->setFiltering(filtering);	
	
	//sceneManager->setPrintEnabled(1024*1024*1024, make_uint3(976,552,0));
	//sceneManager->finishSceneContext();
	//timer.start();
	

	//With ray generation on launch	
	sceneManager->finishSceneContext();
	timer.start();
	
	sceneManager->setRayRange(currentElevation,currentElevation+deltaEl,currentAzimuth,currentAzimuth+deltaAz,rayD,rayD);
	//RayDensityNormalizationSimulation* sim=dynamic_cast<RayDensityNormalizationSimulation*>(sceneManager->getSimulation());
	sim->setInitialDensity(sceneManager->getRaySphere().rayCount,currentAzimuth,currentAzimuth+deltaAz,currentElevation, currentElevation+deltaEl);
	sim->setFiltering(filtering);	
	//sceneManager->addReceiver(1,make_float3(3.72,-22.4871,569.75),polarization, sphereRadius, sceneManager->printPower);
	//sceneManager->addReceiver(1,make_float3(3.72,-22.4788,570.5),polarizationRx, sphereRadius, sceneManager->printPower);
	
	uint nrx=216;	
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
	
	for (int i=0;i<floor(1080/nrx);++i) {
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
			//For square tunnel below, comment above
			//posrx=make_float3(-width/2.0f + 3.72f,-height/2.0f + 2.7f,zinit );
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

				std::cout<<"Tracing angle (el/az)="<<(currentElevation-overlap)<<","<<(currentElevation+deltaEl+overlap)<<"/"<<(currentAzimuth-overlap)<<","<<(currentAzimuth+deltaAz+overlap)<<std::endl;
				//gen->generateRandomUniformOnDevice(currentElevation,currentElevation+deltaEl,currentAzimuth,currentAzimuth+deltaAz,rayD*rayD);
				//sceneManager->createRaySphereFromExternalBuffer(rayD,rayD,gen->getDevicePointer());
				solidAngle=deg2rad*deltaAz*(cosf(deg2rad*currentElevation)-cosf(deg2rad*(currentElevation+deltaEl)));
				rayD=floor(sqrt(rayGoal*solidAngle)); //Missing here the azimuth part of the density...
				sceneManager->setRayRange(currentElevation,currentElevation+deltaEl,currentAzimuth,currentAzimuth+deltaAz,rayD,rayD);
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
void AnversTests::runAnversStaticTruckReceiversRDNIsotropic(bool half, bool computeField) { 
	//float freq = 1.31e9f;
	RayDensityNormalizationSimulation* sim = new RayDensityNormalizationSimulation(sceneManager);
	sceneManager->setSimulation(sim);
	ComputeMode mode=ComputeMode::VOLTAGE;
	if (computeField) {
		std::cout<<"Computing FIELD" <<std::endl;
		mode=ComputeMode::FIELD;
	}
	sim->setComputeMode(mode);
	sceneManager->enableGenerateRaysOnLaunch();	
	//Add diffraction simulation
	SingleDiffraction* simd= new SingleDiffraction(sceneManager);
	sceneManager->setSimulation(simd);
	simd->setComputeMode(mode);
	sceneManager->initContext(frequency);
	simd->setEnableSimulation(false);
	//sceneManager->getSimulation()->setPrintHits(true);


	Timer timer;
	std::cout<<"**** Simulating Anvers tunnel with trucks with isotropic RDN and diffraction ***** "<<std::endl;
	std::cout<<"\tf="<<(frequency/1e6)<< " MHz; tx=" <<tx<<"; polarization="<<polarizationTx<<std::endl;
	std::cout<<"\trx="<<rx<<"; radius="<<sphereRadius;
	if (increaseRadius) {
		std::cout<<"+ "<<sphereDelta;
	}
	std::cout<<"; polarization="<<polarizationRx<<std::endl;
	sceneManager->setMinEpsilon(1e-4f);
        //std::string path("meshes/anvers/trucks");
	//loadAnversScenario(path);
        std::string path("trucks-metal.json");
	loadAnversJsonScenario(path);
	
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
	
	timer.start();
	
	
	
	uint nrx=360;	
	//Receivers, z positions: end of truck18, mid, beginning truck19
	float zpos[] = {998.0f, 1018.0f, 1038.0f};	
	float3 posrx;
	
	uint launches=0;
	//To get the correct z
	float R=2900.0f;
	float l=112.0f; //2*R*sin a/2=2*R*0.04f;
	float h=2.32092875f; //h=R-0.5*sqrt(4*R*R-4*l*l)
	for (int i=0; i<3; i++) { 
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
	
	++launches;
	timer.stop();
	std::cout<<"Time="<<timer.getTime()<<". Time/launch="<<(timer.getTime()/launches)<<std::endl;
	sceneManager->transmit(0, 1.0f, tx, polarizationTx, false);

	
}
void AnversTests::runAnversFlatOvertakingTrucksRDNIsotropic(bool half, bool computeField, bool forward, bool sectorized, bool useReflection, bool useDiffraction) {
	RayDensityNormalizationSimulation* sim = new RayDensityNormalizationSimulation(sceneManager);
	//LPFlatMeshReflectionSimulation* sim= new LPFlatMeshReflectionSimulation(sceneManager);
	//sim->setExecutionMethod(RDNExecutionMode::HITINFO);
	//sim->setEnableTraceLog(true);
	sceneManager->setSimulation(sim);
	ComputeMode mode = ComputeMode::VOLTAGE;
	if (computeField) {
		std::cout << "Computing FIELD" << std::endl;
		mode = ComputeMode::FIELD;
	}
	sim->setComputeMode(mode);
	sim->setEnableSimulation(useReflection);
	sceneManager->enableGenerateRaysOnLaunch();
	//Add diffraction simulation
	SingleDiffraction* simd = new SingleDiffraction(sceneManager);
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
	//std::cout<<"**** Simulating Flat Anvers tunnel with trucks with isotropic RDN and diffraction ***** "<<std::endl;
	//std::cout<<"**** Simulating Flat Anvers tunnel with single middle truck with isotropic RDN and diffraction ***** "<<std::endl;
	//std::cout<<"**** Simulating Flat Anvers empty tunnel with same rx positions and  single iddle truck with isotropic RDN  ***** "<<std::endl;
	//std::cout<<"**** Simulating Flat Anvers tunnel with trucks and cars on free-lane with isotropic RDN and diffraction ***** "<<std::endl;
	if (emptyTunnel) {
		std::cout << "**** Simulating Flat Anvers tunnel empty and cars on free-lane with isotropic RDN and diffraction ***** " << std::endl;
	}
	else {
		std::cout << "**** Simulating Flat Anvers tunnel with trucks and overtaking on free-lane with isotropic RDN and diffraction ***** " << std::endl;
	}
	std::cout << "\tf=" << (frequency / 1e6) << " MHz; tx=" << tx << "; polarization=" << polarizationTx << std::endl;
	std::cout << "\trx=" << rx << "; radius=" << sphereRadius;

	//Set sphere radius for this scenario
	this->sphereDelta = 2e-4;
	if (increaseRadius) {
		std::cout << "+ " << sphereDelta;
	}
	std::cout << "; polarization=" << polarizationRx << std::endl;
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

	float width = 10.2f;
	float height = 6.0f;
	float length = 1070.0f;
	Matrix4x4 tm;
	tm.setRow(0, make_float4(width, 0, 0, width / 2));
	tm.setRow(1, make_float4(0, height, 0, -height / 2));
	tm.setRow(2, make_float4(0, 0, length, length / 2.0f));
	tm.setRow(3, make_float4(0, 0, 0, 1));
	MaterialEMProperties emProp1;
	emProp1.dielectricConstant = make_float2(5.0f, -60.0f * sceneManager->getChannelParameters().waveLength * 0.01f);
	emProp1.tattenuation = make_float2(0.1f, -75.f);
	loadTransformedSquareTunnel(emProp1, tm);
	int rayD = 10000;
	//With ray generation on launch
	sceneManager->finishSceneContext();

	//With LPFlatMeshReflectionSimulation
	//sceneManager->createRaySphere2D(0,0.1,180,-90,0.1,90);

	if (!sectorized) {
		if (half) {
			std::cout << "**** Anvers Tunnel with  Half Sphere RDN ***" << std::endl;
			sceneManager->setRayRange(0.0, 180.0, -90.0, 90.0, rayD, rayD);
		}
		else {
			std::cout << "**** Anvers Tunnel with Isotropic RDN ***" << std::endl;
			sceneManager->setRayRange(0.0, 180.0, 0.0, 360.0, rayD, rayD);
		}

		//sceneManager->createRaySphereFromExternalBuffer(rayD,rayD,gen->getDevicePointer());
		//RayDensityNormalizationSimulation* sim=dynamic_cast<RayDensityNormalizationSimulation*>(sceneManager->getSimulation());
		if (half) {
			sim->setInitialDensity(((float)sceneManager->getRaySphere().rayCount) / (2 * M_PIf));
		}
		else {
			sim->setInitialDensity(((float)sceneManager->getRaySphere().rayCount) / (4 * M_PIf));
		}
		sim->setFiltering(filtering);
	}

	timer.start();
	int gainId;
	if (this->useAntennaGain && !emptyTunnel) {
		if (forward) {
			AntennaGain gains = sceneManager->loadGainsFromFileIndBPower("forward.txt", true);
			gainId = sceneManager->registerAntennaGain(gains);
		}
		else {
			AntennaGain gains = sceneManager->loadGainsFromFileIndBPower("backward.txt", true);
			gainId = sceneManager->registerAntennaGain(gains);
		}
		//sceneManager->registerTransmitterGain(0,gainId);
	}


	//sceneManager->addReceiver(1,make_float3(3.2,(2.7-6.0), 808),polarizationRx,0.5, sceneManager->printPower);
	//if (this->useAntennaGain) {
	//	sceneManager->registerReceiverGain(1,gainId);
	//}
	//sceneManager->transmit(0, 1.0f, tx, polarizationTx, false);

//uint nrx=380;	

//Free lane
	uint nrx = 1040;
	for (int i = 1; i <= nrx; ++i) {
		//sceneManager->addReceiver(i,make_float3(3.72f,2.7f-6.0f, 20.0f),polarization, sphereRadius, sceneManager->printPower);
		sceneManager->addReceiver(i, rx, polarizationRx, sphereRadius, sceneManager->printPower);
		if (this->useAntennaGain && !emptyTunnel) {
			sceneManager->registerReceiverGain(i, gainId);
		}
	}

	float zinit = 498.0f;
	//this->zStep = 0.1f;
	this->zStep = 0.07f;
	uint launches = 0;
	//To get the correct z
	float R = 2900.0f;
	float l = 112.0f; //2*R*sin a/2=2*R*0.04f;
	float h = 2.32092875f; //h=R-0.5*sqrt(4*R*R-4*l*l)
	int zgroup = 13;
	int k = 0;
	float tl = 0.0f;
	for (int i = 0; i < floor(1040 / nrx); ++i) {
		float3 posrx;
		for (int j = 1; j <= nrx; ++j) {
			//if (k==39) {
			//	zinit += zgroup;
			//	k=0;
			//}
			//Cars on Free lane
			posrx = make_float3(3.2f, (2.7f - 6.0f), zinit);
			//Cars between trucks
			//posrx=make_float3(7.0f,(2.7f-6.0f),zinit );

			//if (zinit<=426) {
			//	posrx=make_float3(7.0f,(2.7f-6.0f)-(0.04f*zinit),zinit );
			//} else if ((zinit>426) && (zinit<=658)) {
			//	float zp=zinit-538; //z-426+112
			//	float yp=sqrt((R*R)-(zp*zp))-R+h; 
			//	float fl=-(23.04+yp); //Floor Y coordinate on curved section
			//	posrx=make_float3(7.0f, fl+2.7f,zinit );
			//} else {

			//	posrx=make_float3(7.0f,(2.7f-23.04f)+(0.04f*(zinit-658)),zinit );
			//}
			if (increaseRadius) {
				sceneManager->updateReceiver(j, posrx, sphereRadius);
				sphereRadius += sphereDelta;
			}
			else {
				sceneManager->updateReceiver(j, posrx);
			}
			zinit += zStep;
			k++;
		}
		if (sectorized) {
			SphereScanConfiguration c;
			tl += runSectorizedLaunch(c, 0, 1.0, tx, polarizationTx, sim);
		}
		else {
			//First launch
			//sceneManager->transmit(0, 1.0f, postx, polarizationTx, false);
			sceneManager->transmit(0, 1.0f, tx, polarizationTx, false);

		}

		++launches;
	}
	if (sectorized) {
		std::cout << "Total time=" << tl << ". Time/launch=" << (tl / launches) << std::endl;
	}
	else {
		timer.stop();
		std::cout << "Time=" << timer.getTime() << ". Time/launch=" << (timer.getTime() / launches) << std::endl;
	}

}

void AnversTests::runAnversFlatStaticBetweenTrucksRDNIsotropic(bool half, bool computeField, bool forward,  bool sectorized, bool useReflection, bool useDiffraction) { 
	
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
	//std::cout<<"**** Simulating Flat Anvers tunnel with trucks with isotropic RDN and diffraction ***** "<<std::endl;
	//std::cout<<"**** Simulating Flat Anvers tunnel with single middle truck with isotropic RDN and diffraction ***** "<<std::endl;
	//std::cout<<"**** Simulating Flat Anvers empty tunnel with same rx positions and  single iddle truck with isotropic RDN  ***** "<<std::endl;
	//std::cout<<"**** Simulating Flat Anvers tunnel with trucks and cars on free-lane with isotropic RDN and diffraction ***** "<<std::endl;
	if (emptyTunnel ) {
		std::cout<<"**** Simulating Flat Anvers tunnel empty and cars on free-lane with isotropic RDN and diffraction ***** "<<std::endl;
	} else {
		std::cout<<"**** Simulating Flat Anvers tunnel with trucks and cars on free-lane with isotropic RDN and diffraction ***** "<<std::endl;
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
	uint nrx=252;	
	if (useRDN) {
	 	nrx=1008;
	} else {
		sphereDelta=2e-4;
	}	
	for (int i=1;i<=nrx;++i) {
		//sceneManager->addReceiver(i,make_float3(3.72f,2.7f-6.0f, 20.0f),polarization, sphereRadius, sceneManager->printPower);
		sceneManager->addReceiver(i,rx,polarizationRx, sphereRadius, sceneManager->printPower);
		if (this->useAntennaGain && !emptyTunnel) {
			sceneManager->registerReceiverGain(i,gainId);
		}
	}
	
	float zinit=63.0f;
	uint launches=0;
	//To get the correct z
	float R=2900.0f;
	float l=112.0f; //2*R*sin a/2=2*R*0.04f;
	float h=2.32092875f; //h=R-0.5*sqrt(4*R*R-4*l*l)
	int zgroup=13;
	int k=0;
	float tl=0.0f;
	for (int i=0;i<floor(1008/nrx);++i) {
		float3 posrx;
		for (int j=1;j<=nrx;++j) {
			//if (k==39) {
			//	zinit += zgroup;
			//	k=0;
			//}
			//Cars on Free lane
			posrx=make_float3(3.2f,(2.7f-6.0f),zinit );
			//Cars between trucks
			//posrx=make_float3(7.0f,(2.7f-6.0f),zinit );

			//if (zinit<=426) {
			//	posrx=make_float3(7.0f,(2.7f-6.0f)-(0.04f*zinit),zinit );
			//} else if ((zinit>426) && (zinit<=658)) {
			//	float zp=zinit-538; //z-426+112
			//	float yp=sqrt((R*R)-(zp*zp))-R+h; 
			//	float fl=-(23.04+yp); //Floor Y coordinate on curved section
			//	posrx=make_float3(7.0f, fl+2.7f,zinit );
			//} else {

			//	posrx=make_float3(7.0f,(2.7f-23.04f)+(0.04f*(zinit-658)),zinit );
			//}
			if (increaseRadius) {
				sceneManager->updateReceiver(j, posrx,sphereRadius);
				sphereRadius += sphereDelta;
			} else {
				sceneManager->updateReceiver(j, posrx);
			}
			zinit += zStep;
			k++;
		}
		if (sectorized) {
				SphereScanConfiguration c;
				//throw  opal::Exception("AnversTests::runAnversFlatStaticBetweenTrucksRDNIsotropic(): not implemented sectorized for LPFlatMeshReflectionSimulation yet");
				tl +=runSectorizedLaunch(c, 0, 1.0, tx, polarizationTx, sim);
		} else {
			//First launch
			//sceneManager->transmit(0, 1.0f, postx, polarizationTx, false);
			sceneManager->transmit(0, 1.0f, tx, polarizationTx, false);

		}

		++launches;
	} 
	if (sectorized) {
		std::cout<<"Total time="<<tl<<". Time/launch="<<(tl/launches)<<std::endl;
	} else {	
		timer.stop();
		std::cout<<"Time="<<timer.getTime()<<". Time/launch="<<(timer.getTime()/launches)<<std::endl;
	}
	
}
void AnversTests::runAnversStaticBetweenTrucksRDNIsotropic(bool half, bool computeField, bool forward) { 
	//float freq = 1.31e9f;
	RayDensityNormalizationSimulation* sim = new RayDensityNormalizationSimulation(sceneManager);
	//sim->setExecutionMethod(RDNExecutionMode::HITINFO);
	//sim->setEnableTraceLog(true);
	sceneManager->setSimulation(sim);
	ComputeMode mode=ComputeMode::VOLTAGE;
	if (computeField) {
		std::cout<<"Computing FIELD" <<std::endl;
		mode=ComputeMode::FIELD;
	}
	sim->setComputeMode(mode);
	sceneManager->enableGenerateRaysOnLaunch();	
	//Add diffraction simulation
	SingleDiffraction* simd= new SingleDiffraction(sceneManager);
	sceneManager->setSimulation(simd);
	simd->setComputeMode(mode);
	if (this->useAntennaGain) {
		sceneManager->setUseAntennaGain(true);
	}
	sceneManager->initContext(frequency);
	//simd->setEnableSimulation(false);
	//sceneManager->getSimulation(0)->setPrintHits(true);


	Timer timer;
	std::cout<<"**** Simulating Anvers tunnel with trucks with isotropic RDN and diffraction ***** "<<std::endl;
	std::cout<<"\tf="<<(frequency/1e6)<< " MHz; tx=" <<tx<<"; polarization="<<polarizationTx<<std::endl;
	std::cout<<"\trx="<<rx<<"; radius="<<sphereRadius;
	if (increaseRadius) {
		std::cout<<"+ "<<sphereDelta;
	}
	std::cout<<"; polarization="<<polarizationRx<<std::endl;
	sceneManager->setMinEpsilon(1e-4f);
        //std::string path("meshes/anvers/trucks");
	//loadAnversScenario(path);
        std::string path("trucks-metal.json");
	loadAnversJsonScenario(path);
	
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
	
	timer.start();
	int gainId;	
	if (this->useAntennaGain) {
		if (forward) {
			AntennaGain gains=sceneManager->loadGainsFromFileIndBPower("forward.txt", true);
			gainId=sceneManager->registerAntennaGain(gains);
		} else {
			AntennaGain gains=sceneManager->loadGainsFromFileIndBPower("backward.txt", true);
			gainId=sceneManager->registerAntennaGain(gains);
		}
		//sceneManager->registerTransmitterGain(0,gainId);
	}
	
	
	//	sceneManager->addReceiver(1,make_float3(7,-3.9, 1069),polarizationRx, sphereRadius, sceneManager->printPower);
	//	if (this->useAntennaGain) {
	//		sceneManager->registerReceiverGain(1,gainId);
	//	}
	//	sceneManager->transmit(0, 1.0f, tx, polarizationTx, false);
	
	uint nrx=380;	
	for (int i=1;i<=nrx;++i) {
		//sceneManager->addReceiver(i,make_float3(3.72f,2.7f-6.0f, 20.0f),polarization, sphereRadius, sceneManager->printPower);
		sceneManager->addReceiver(i,rx,polarizationRx, sphereRadius, sceneManager->printPower);
		if (this->useAntennaGain) {
			sceneManager->registerReceiverGain(i,gainId);
		}
	}
	
	float zinit=63.0f;
	uint launches=0;
	//To get the correct z
	float R=2900.0f;
	float l=112.0f; //2*R*sin a/2=2*R*0.04f;
	float h=2.32092875f; //h=R-0.5*sqrt(4*R*R-4*l*l)
	int zgroup=13;
	int k=0;
	for (int i=0;i<floor(760/nrx);++i) {
		float3 posrx;
		for (int j=1;j<=nrx;++j) {
			if (k==39) {
				zinit += zgroup;
				k=0;
			}
			

			if (zinit<=426) {
				posrx=make_float3(7.0f,(2.7f-6.0f)-(0.04f*zinit),zinit );
			} else if ((zinit>426) && (zinit<=658)) {
				float zp=zinit-538; //z-426+112
				float yp=sqrt((R*R)-(zp*zp))-R+h; 
				float fl=-(23.04+yp); //Floor Y coordinate on curved section
				posrx=make_float3(7.0f, fl+2.7f,zinit );
			} else {

				posrx=make_float3(7.0f,(2.7f-23.04f)+(0.04f*(zinit-658)),zinit );
			}
			if (increaseRadius) {
				sceneManager->updateReceiver(j, posrx,sphereRadius);
				sphereRadius += sphereDelta;
			} else {
				sceneManager->updateReceiver(j, posrx);
			}
			zinit += zStep;
			k++;
		}
		//First launch
		//sceneManager->transmit(0, 1.0f, postx, polarizationTx, false);
		sceneManager->transmit(0, 1.0f, tx, polarizationTx, false);

		++launches;

	} 
	
	timer.stop();
	std::cout<<"Time="<<timer.getTime()<<". Time/launch="<<(timer.getTime()/launches)<<std::endl;
	
}
void AnversTests::runAnversStaticTruckRDNIsotropic(bool half, bool computeField) { 
	//float freq = 1.31e9f;
	RayDensityNormalizationSimulation* sim = new RayDensityNormalizationSimulation(sceneManager);
	sceneManager->setSimulation(sim);
	ComputeMode mode=ComputeMode::VOLTAGE;
	if (computeField) {
		std::cout<<"Computing FIELD" <<std::endl;
		mode=ComputeMode::FIELD;
	}
	sim->setComputeMode(mode);
	sceneManager->enableGenerateRaysOnLaunch();	
	//Add diffraction simulation
	SingleDiffraction* simd= new SingleDiffraction(sceneManager);
	sceneManager->setSimulation(simd);
	simd->setComputeMode(mode);
	sceneManager->initContext(frequency);
	simd->setEnableSimulation(false);
	//sceneManager->getSimulation()->setPrintHits(true);


	Timer timer;
	std::cout<<"**** Simulating Anvers tunnel with trucks with isotropic RDN and diffraction ***** "<<std::endl;
	std::cout<<"\tf="<<(frequency/1e6)<< " MHz; tx=" <<tx<<"; polarization="<<polarizationTx<<std::endl;
	std::cout<<"\trx="<<rx<<"; radius="<<sphereRadius;
	if (increaseRadius) {
		std::cout<<"+ "<<sphereDelta;
	}
	std::cout<<"; polarization="<<polarizationRx<<std::endl;
	sceneManager->setMinEpsilon(1e-4f);
        //std::string path("meshes/anvers/trucks");
	//loadAnversScenario(path);
        std::string path("trucks-metal.json");
	loadAnversJsonScenario(path);
	
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
	
	timer.start();
	
	
	
	uint nrx=360;	
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
	
	for (int i=0;i<floor(1080/nrx);++i) {
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
		//sceneManager->transmit(0, 1.0f, postx, polarizationTx, false);
		sceneManager->transmit(0, 1.0f, tx, polarizationTx, false);

		++launches;

	} 
	
	timer.stop();
	std::cout<<"Time="<<timer.getTime()<<". Time/launch="<<(timer.getTime()/launches)<<std::endl;
	
}

//Curved empty Anvers with RDN
void AnversTests::runRDNIsotropic(bool half, bool computeField) {
	//float freq = 1.31e9f;
	RayDensityNormalizationSimulation* sim = new RayDensityNormalizationSimulation(sceneManager);
	sceneManager->setSimulation(sim);
	ComputeMode mode=ComputeMode::VOLTAGE;
	if (computeField) {
		std::cout<<"Computing FIELD" <<std::endl;
		mode=ComputeMode::FIELD;
	}
	sim->setComputeMode(mode);
	sceneManager->enableGenerateRaysOnLaunch();	
	if (this->useAntennaGain) {
		sceneManager->setUseAntennaGain(true);
	}
	sceneManager->initContext(frequency);
	//sceneManager->getSimulation()->setPrintHits(true);


	//loadStraightAnversTunnel(sceneManager);
	loadAnversTunnel();
	Timer timer;
	std::cout<<"**** Simulating Anvers tunnel with isotropic RDN ***** "<<std::endl;
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
	

	int gainId;	
	if (this->useAntennaGain) {
		AntennaGain gains=sceneManager->loadGainsFromFileIndBPower("dipole.txt");
		gainId=sceneManager->registerAntennaGain(gains);
		sceneManager->registerTransmitterGain(0,gainId);
	}
	uint nrx=360;	
	for (int i=1;i<=nrx;++i) {
		//sceneManager->addReceiver(i,make_float3(3.72f,2.7f-6.0f, 20.0f),polarization, sphereRadius, sceneManager->printPower);
		sceneManager->addReceiver(i,rx,polarizationRx, sphereRadius, sceneManager->printPower);
		if (this->useAntennaGain) {
			sceneManager->registerReceiverGain(i,gainId);
		}
	}
	
	
	float zinit=20.0f;
	uint launches=0;
	//To get the correct z
	float R=2900.0f;
	float l=112.0f; //2*R*sin a/2=2*R*0.04f;
	float h=2.32092875f; //h=R-0.5*sqrt(4*R*R-4*l*l)
	
	for (int i=0;i<floor(1080/nrx);++i) {
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
		//sceneManager->transmit(0, 1.0f, postx, polarizationTx, false);
		sceneManager->transmit(0, 1.0f, tx, polarizationTx, false);

		++launches;

	} 
	
	timer.stop();
	std::cout<<"Time="<<timer.getTime()<<". Time/launch="<<(timer.getTime()/launches)<<std::endl;
	
}
void AnversTests::runRDNSingleRx(bool half) {
	float freq = 1.31e9f;
	RayDensityNormalizationSimulation* sim = new RayDensityNormalizationSimulation(sceneManager);
	sceneManager->setSimulation(sim);
	sceneManager->initContext(freq);
	//sceneManager->getSimulation()->setPrintHits(true);


	//uuuuuuuaightAnversTunnel(sceneManager);
	loadAnversTunnel();

	Timer timer;
	//Receiver polarization
	//optix::float3 polarization = H; 	
	//optix::float3 polarization = pol; 	
	
	//Transmitter
	//optix::float3 postx = make_float3(0.50f,2.52f-6.0f, 0.0f);
	//float3 polarizationTx = H; 
	//float3 polarizationTx = pol; 
	//float overlap=0.0f;
	std::cout<<"**** Simulating Anvers tunnel with isotropic RDN ***** "<<std::endl;
	std::cout<<"Transmitting at "<<tx<<" with polarization="<<polarizationTx<<std::endl;
	std::cout<<"Receiving with radius="<<sphereRadius<<" with polarization="<<polarizationRx<<std::endl;
	sceneManager->setMinEpsilon(1e-4f);
	
	int rayD=10000;
	OpalRaySphereGenerator* gen=sceneManager->getRaySphereGenerator();
	if (half) { 
		gen->generateRandomUniformOnDevice(0.0,180.0,-90.0,90.0,rayD*rayD); 
		std::cout <<"**** Anvers Tunnel with  Half Sphere RDN ***"<<std::endl;	
	} else {
		gen->generateRandomUniformSphereOnDevice(rayD*rayD);
		std::cout <<"**** Anvers Tunnel with Isotropic RDN ***"<<std::endl;	
	}
	sceneManager->createRaySphereFromExternalBuffer(rayD,rayD,gen->getDevicePointer());
	//RayDensityNormalizationSimulation* sim=dynamic_cast<RayDensityNormalizationSimulation*>(sceneManager->getSimulation());
	if (half) { 
		sim->setInitialDensity(((float)sceneManager->getRaySphere().rayCount)/(2*M_PIf));
	} else {
		sim->setInitialDensity(((float)sceneManager->getRaySphere().rayCount)/(4*M_PIf));
	}

	
	
	//sceneManager->setPrintEnabled(1024*1024*1024, make_uint3(976,552,0));
	sceneManager->finishSceneContext();
	timer.start();
	
	uint launches=1;

	//float sr=0.1460013241f;
	//float sr=0.1188008636f;
	float sr=0.3050040007;
	float3 postx=make_float3(0.50f,2.52f-6.0f, 0.0f);
	//float3 posrx=make_float3(3.72f,-22.6499424f,530 );
	//float3 posrx=make_float3(3.72f,-20.574014f,428 );
	float3 posrx=make_float3(3.72f,-1.61f,1126.25f );
	sceneManager->addReceiver(1,posrx,polarizationRx, sr, sceneManager->printPower);
	sceneManager->transmit(0, 1.0f, postx, polarizationTx, false);
	
	timer.stop();
	std::cout<<"Time="<<timer.getTime()<<". Time/launch="<<(timer.getTime()/launches)<<std::endl;
	
}


