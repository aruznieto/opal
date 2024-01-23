/***************************************************************/
//
//Copyright (c) 2021 Esteban Egea-Lopez http://girtel.upct.es/~simulations
//
/**************************************************************/


#ifndef LORA_H
#define LORA_H
#include "../Opal.h"
#include "tests.h"
#include <iostream>
#include "../basicSimulation.h"
#include "../flatSimulation.h"
#include "../singleDiffraction.h"
#include "../rayDensityNormalizationSimulation.h"
using namespace opal;
namespace opal {
	class Lora : public BasicTests {
		protected:
			float freq;
			float txPower;
			float3 baseStation;
			std::string gsufix;
			uint replications;
			bool useGain;
			RayDensityNormalizationSimulation* rdnSim ;
			LPFlatMeshReflectionSimulation* lpflatSim;
			BasicFlatMeshReflectionSimulation* basicSim;
			SingleDiffraction* diffSim;
			int gIdTx;
			int gIdRx;
		public:
			Lora(OpalSceneManager* sceneManager, float sphereRadius, bool useDepolarization); 
			void runPoint(bool addRandom, bool useRDN, std::string path); 
			void loadScenarioLora(bool addRandom, bool useRDN, std::string path); 
			std::vector<float3> loadReceiversFromFile(std::string file); 
			void runTests(std::string test);
			void runAll(bool empty,std::string test, std::string folder) ;
			void runAllMoving(bool empty,std::string test, std::string folder) ;
			void test();
			void  buildAllScenario(bool useRDN,std::string path ) ;
			void  runScenarioVariant(bool addRandom, bool useRDN,  std::string folder) ;
			void  runScenarioVariantMovingTransmitter(bool addRandom, bool useRDN,  std::string folder) ;

	};

}
#endif


