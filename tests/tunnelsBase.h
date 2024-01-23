
/***************************************************************/
//
//Copyright (c) 2019 Esteban Egea-Lopez http://girtel.upct.es/~simulations
//
/**************************************************************/
#ifndef TUNNELSBASE_H
#define TUNNELSBASE_H
#include "../rayDensityNormalizationSimulation.h"
#include "../Opal.h"
#include <string>
namespace opal {
	struct SphereScanConfiguration {
		float initElevation;
		float initAzimuth;
		float endElevation;
		float endAzimuth;
		float deltaEl;
		float deltaAz;
		float asEl;
		float asAz;
		float overlap;
		float rayGoal;//For RDN
		int filtering; //For RDN
		SphereScanConfiguration() {
		//Default sensible values for tunnel	
			this->initElevation=0.0f;
			this->initAzimuth=-90.0f;
			this->endElevation=180.0f;
			this->endAzimuth=90.0f;
			this->deltaEl=10.0f;
			this->deltaAz=10.0f;
			this->asEl=0.01f;
			this->asAz=0.01f;
			this->overlap=0.0f;
			this->rayGoal=1e9;
			this->filtering=2;

		} 
	};
	class TunnelsBase {
		protected:
	//Sphere scanning
			OpalSceneManager* sceneManager;
			float sphereRadius;
			bool useDepolarization;
			void loadCircularTunnel(float radius,  float length, MaterialEMProperties emProp1);
			void loadSquareTunnel(float width, float height, float length, MaterialEMProperties emProp1);
			void loadTransformedSquareTunnel(MaterialEMProperties emProp1, optix::Matrix4x4 tm);
			void loadHalfCylinder(float radius,  float length, float height, MaterialEMProperties emProp1);
			float runSectorizedLaunch(SphereScanConfiguration conf, int txId, float txPower, float3 txpos, float3 polarizationTx, OpalSimulation* sim ); 
		public:
			TunnelsBase(OpalSceneManager* sceneManager, float sphereRadius, bool useDepolarization); 
			std::vector<int> parseTestString(std::string test);
	};
}
#endif
