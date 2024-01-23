/***************************************************************/
//
//Copyright (c) 2023 Esteban Egea-Lopez http://girtel.upct.es/~simulations
//
/**************************************************************/


#ifndef EXPOSURE_H
#define EXPOSURE_H
#include "../Opal.h"
#include "lille.h"
#include <iostream>
using namespace opal;
namespace opal {
	class Exposure : public Lille {
		protected:
			float peakFraction;
			float lowFraction;
			float fraction;
			std::string sensors_file;
			int maxUsersPerCell;
			int maxFreq;
		public:
			Exposure(OpalSceneManager* sceneManager, float sphereRadius, bool useDepolarization); 
			void configureSimulations(bool useRDN ) ;
			void setLaunchSizeAndFinish(bool useRDN ) ;
			void  executeExposure(bool addRandom, bool useRDN, std::string scenarioPath, std::string outputFile, std::string txPath) ;
			std::vector<float3> sampleActiveReceivers(std::vector<float3>& grid, Basestation base, float distance, float fraction); 
			void run5GBeamSelection(std::string test, std::string outputFile);
			std::vector<float3> getPointsAtDistance(std::vector<float3>& grid,float distance, float3 tx); 
			std::vector<int> assignActiveReceiversToBeam(std::vector<float3>& au, Basestation base) ; 
			virtual std::vector<Basestation> loadBasestationsFromFile(std::string file, bool replicatePerDiagram);
			float getMeanBasestationDistace(std::vector<Basestation>& bases); 
			void saveData(std::string outputFile, std::string header, std::string data); 
			void 	run5GOneOperator(std::string test, std::string outputFile); 
			void  executeExposureOneOperator(bool addRandom, bool useRDN, std::string scenarioPath, std::string outputFile, std::string txPath); 

	};

}
#endif


