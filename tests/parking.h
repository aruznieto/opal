/***************************************************************/
//
//Copyright (c) 2021 Esteban Egea-Lopez http://girtel.upct.es/~simulations
//
/**************************************************************/


#ifndef PARKING_H
#define PARKING_H
#include "../Opal.h"
#include "tests.h"
#include <iostream>
using namespace opal;
namespace opal {
	class Parking : public BasicTests {
		protected:
			float freq;
			float3 baseStation;
			std::string gsufix;
			uint replications;
			bool useGain;
		public:
			//Horizontal plane test. To validate against a two-ray model
			Parking(OpalSceneManager* sceneManager, float sphereRadius, bool useDepolarization); 
			void  loadScenarioParkingMultifrequency(bool addRandom, bool useRDN, std::string path, std::vector<float>& frequencies, std::vector<float>& radius, std::string outputFile) ;
			void loadScenarioParking(bool addRandom, bool useRDN, std::string path); 
			std::vector<float3> loadReceiversFromFile(std::string file); 
			void runTests(std::string test, std::string outputFile);
			void test();

	};

}
#endif


