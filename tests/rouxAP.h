/***************************************************************/
//
//Copyright (c) 2019 Esteban Egea-Lopez http://girtel.upct.es/~simulations
//
/**************************************************************/


#ifndef ROUXAP_H
#define ROUXAP_H
#include "tunnelsBase.h"
#include <string>
namespace opal {
	class RouxAPTests: TunnelsBase {
		protected:
			bool increaseRadius;
			float sphereDelta;
			float frequency;
			float3 rx;
			float3 tx;
			float3 polarizationTx;
			float3 polarizationRx;
			float discriminateAngle; //For CURVEDFLATWALLS
			unsigned int filtering; //For RDN
			float sphereDistanceError;
			void loadRouxTunnel();
			void loadCubeEquivalentTunnel();	
		public:
			RouxAPTests(OpalSceneManager* sceneManager, float sphereRadius, bool useDepolarization); 
			void runADSingleRay();
			void runADCube();
			void runAD();
			void runRDNCube();
			void runRDNIsotropic();
			void runRDN();
			void runRDNFixedRay();
			void runTests(std::string test);
	};
}
#endif



