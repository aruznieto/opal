/***************************************************************/
//
//Copyright (c) 2021 Esteban Egea-Lopez http://girtel.upct.es/~simulations
//
/**************************************************************/


#ifndef CARTAGENAGSM_H
#define CARTAGENAGSM_H
#include "../Opal.h"
#include "tests.h"
#include <iostream>
using namespace opal;
namespace opal {
	class CartagenaGSM : public BasicTests {
		public:
			//Horizontal plane test. To validate against a two-ray model
			CartagenaGSM(OpalSceneManager* sceneManager, float sphereRadius, bool useDepolarization); 
			void loadScenarioCartagena(bool addRandom, bool useRDN); 
			std::vector<float2> loadReceiversFromFile(std::string file); 
			void runTests(std::string test);

	};

}
#endif


