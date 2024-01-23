/***************************************************************/
//
//Copyright (c) 2019 Esteban Egea-Lopez http://girtel.upct.es/~simulations
//
/**************************************************************/


#ifndef ANVERS_H
#define ANVERS_H
#include "tunnelsBase.h"
#include <string>
namespace opal {
	class AnversTests: public TunnelsBase {
		protected:
			bool useCurvedSection;
			bool useAntennaGain;
			bool increaseRadius;
			bool tunnelAsFlat;
			bool emptyTunnel;
			bool useRDN;
			float sphereDelta;
			float frequency;
			float zStep;
			float3 rx;
			float3 tx;
			float3 polarizationTx;
			float3 polarizationRx;
			float discriminateAngle; //For CURVEDFLATWALLS
			unsigned int filtering; //For RDN
			void loadStraightAnversTunnel();
			void loadAnversTunnel();	
			void loadAnversScenario(std::string path);
			void loadAnversJsonScenario(std::string path);
		public:
			AnversTests(OpalSceneManager* sceneManager, float sphereRadius, bool useDepolarization); 
			virtual void runTests(std::string test);
			void runTestsJuly21(std::string test);
			void runTestsDepolarizationLeanAntenna(std::string test, bool useGain);
			void runTestsLeanAntenna(std::string test, bool useGain);
			void runTestsBetweenLeanAntenna(std::string test, bool useGain) ;
			//The Anvers tunnel with the curved section of the slope
			virtual void runTunnel();
			//A cubic tunnel with the dimensions of Anvers (no down slope)
			void runStraightTunnel();
			void runRDNStraightTunnel(bool half);
			void runAnversStaticTruckReceiversRDNIsotropic(bool half, bool computeField);
			void runAnversStaticTruckRDNIsotropic(bool half, bool computeField);
			void runAnversFlatStaticBetweenTrucksRDNIsotropic(bool half, bool computeField, bool forward, bool sectorized, bool useReflection, bool useDiffraction);
			void runAnversStaticBetweenTrucksRDNIsotropic(bool half, bool computeField, bool forward);
			//The Anvers tunnels with the slope but all sections straight
			void runNoCurvedSection();
			void runRDN();
			virtual void runRDNIsotropic(bool half, bool computeField);
			void runRDNSingleRx(bool half);
			virtual void runTests6GHz(std::string test, bool useGain); 
			virtual void runAnversFlatOvertakingTrucksRDNIsotropic(bool half, bool computeField, bool forward, bool sectorized, bool useReflection, bool useDiffraction);
	};
}
#endif



