/***************************************************************/
//
//Copyright (c) 2019 Esteban Egea-Lopez http://girtel.upct.es/~simulations
//
/**************************************************************/


#ifndef ANVERSMIMO_H
#define ANVERSMIMO_H
#include "anvers.h"
#include <string>
namespace opal {
	class AnversMimo: public AnversTests {
		protected:
			float delta_tx;
			float delta_rx;
			int mimotx;
			int mimotx_x;
			int mimotx_y;
			int mimorx;
		public:
			AnversMimo(OpalSceneManager* sceneManager, float sphereRadius, bool useDepolarization); 
			virtual void runLille(std::string test, std::string outputPath);
			virtual void runTests(std::string test, std::string outputPath);
			virtual void runFlatEmpty(bool half, bool computeField, bool useReflection, bool useDiffraction, bool sectorized, bool forward, bool multitransmitter); 
			virtual void runFlatTrucks(bool half, bool computeField, bool useReflection, bool useDiffraction, bool sectorized, bool forward, bool multitransmitter, std::string filePath); 
			virtual void runFlatCarTrucksLille(bool half, bool computeField, bool useReflection, bool useDiffraction, bool sectorized, bool forward, bool multitransmitter, std::string filePath); 
			virtual void runFlatTrucksLille(bool half, bool computeField, bool useReflection, bool useDiffraction, bool sectorized, bool forward, bool multitransmitter, std::string filePath); 
	};
}
#endif



