/***************************************************************/
//
//Copyright (c) 2021 Esteban Egea-Lopez http://girtel.upct.es/~simulations
//
/**************************************************************/


#ifndef LORADOPPLER_H
#define LORADOPPLER_H
#include "../Opal.h"
#include "tests.h"
#include "../configuration.h"
#include <iostream>
#include "../basicSimulation.h"
#include "../flatSimulation.h"
#include "../singleDiffraction.h"
#include "../rayDensityNormalizationSimulation.h"
using namespace opal;
namespace opal {
	class LoraDoppler : public BasicTests {
		protected:
			float freq;
			float txPower;
			float3 baseStation;
			std::string gsufix;
			uint replications;
			bool useGain;
			float minEpsilon;
            bool useDepolarization;
            bool useRandom;
            float sphereRadius;
			RayDensityNormalizationSimulation* rdnSim ;
			LPFlatMeshReflectionSimulation* lpflatSim;
			BasicFlatMeshReflectionSimulation* basicSim;
			SingleDiffraction* diffSim;
			int gIdTx;
			int gIdRx;
            ConfigurationParser* config;
		public:
			LoraDoppler(OpalSceneManager* sceneManager, ConfigurationParser* config);
			void runPoint(bool addRandom, bool useRDN, std::string path);
			void loadScenarioLora(bool addRandom, bool useRDN, std::string path);
			std::vector<float3> loadReceiversFromFile(std::string file);
            std::vector<float3> loadReceiversFromJSON();
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


