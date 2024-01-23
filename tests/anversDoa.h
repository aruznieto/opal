/***************************************************************/
//
//Copyright (c) 2021 Esteban Egea-Lopez http://girtel.upct.es/~simulations
//
/**************************************************************/


#ifndef ANVERSDOA_H
#define ANVERSDOA_H
#include "anvers.h"
#include <string>
namespace opal {
	class AnversDoaTests: public AnversTests {
		public:
			AnversDoaTests(OpalSceneManager* sceneManager, float sphereRadius, bool useDepolarization); 
			virtual void runTests(std::string test) override;
			virtual void runTunnel() override;
			virtual void runIsotropicTunnel(bool curvedSim );
			virtual void runRDNIsotropic(bool half, bool computeField) override;
			virtual void runTests6GHz(std::string test, bool useGain) override;
			void runDoaFlatTrucks(bool half, bool computeField, bool forward,  bool sectorized, bool useReflection, bool useDiffraction); 
	};
}
#endif



