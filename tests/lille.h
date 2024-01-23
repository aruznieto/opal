/***************************************************************/
//
//Copyright (c) 2021 Esteban Egea-Lopez http://girtel.upct.es/~simulations
//
/**************************************************************/


#ifndef LILLE_H
#define LILLE_H
#include "../Opal.h"
#include "tests.h"
#include <iostream>
using namespace opal;
namespace opal {
	struct Basestation {
		int id;
		int antenna_id;
		float3 postx;
		float azimuth;
		float freq;
		float freq_end;
		int site;
		std::string tech;
		float lat;
		float lon;
		int antenna_index;
		std::string toString() const;
	};
	class Lille : public BasicTests {
		protected:
			float freq;
			float txPower;
			float3 baseStation;
			std::string gsufix;
			std::string rx_file;
			uint replications;
			bool useGain;
			float delta_rx;
			float delta_tx_x;
			float delta_tx_y;
			uint mimotx;
			uint mimorx;
		public:
			//Horizontal plane test. To validate against a two-ray model
			Lille(OpalSceneManager* sceneManager, float sphereRadius, bool useDepolarization); 
			void loadScenario(bool addRandom, bool useRDN, std::string path, std::string outputFile); 
			void loadWazem(bool addRandom, bool useRDN, std::string path, std::string outputFile);
			void runBasestations(std::string test, std::string outputFile);
			void runWaz(std::string test, std::string outputFile);
			void  executeSingleBasestation(bool addRandom, bool useRDN, std::string scenarioPath, std::string outputFile, std::string txPath) ;
			void  executeBasestation(bool addRandom, bool useRDN, std::string scenarioPath, std::string outputFile, std::string txPath) ;
			void  executeBasestationBeamforming(bool addRandom, bool useRDN, std::string scenarioPath, std::string outputFile, std::string txPath) ;
			std::vector<float3> loadReceiversFromFile(std::string file); 
			virtual std::vector<Basestation> loadTransmittersFromFile(std::string file);
			void testAntennaOrientation(); 
			virtual void runTests(std::string test, std::string outputFile);
			void test();
			void  angleOfDeparture(bool addRandom, bool useRDN, std::string path, std::string outputFile) ;

	};

}
#endif


