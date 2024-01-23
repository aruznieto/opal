/***************************************************************/
//
//Copyright (c) 2019 Esteban Egea-Lopez http://girtel.upct.es/~simulations
//
/**************************************************************/


#ifndef TESTS_H
#define TESTS_H
#include "../Opal.h"
#include "../configuration.h"
#include <iostream>

//Basic tests to show and validate the features of Opal and how to use them. 
//In the implementation file you can find additional comments and explanations
//This class is also used as base class for additional tests
using namespace opal;
namespace opal {
	class BasicTests {
		protected:
			OpalSceneManager* sceneManager;
			ConfigurationParser* config;
			float sphereRadius;
			bool useDepolarization;
			

		public:
			BasicTests(OpalSceneManager* sceneManager, float sphereRadius, bool useDepolarization);
		        virtual ~BasicTests();	
			//Set the parameters from a configuration file
			void setConfiguration(ConfigurationParser* config);
			//Horizontal plane test. To validate against a two-ray model
			void planeTest(int mode);
			//Horizontal plane test with different frequencies. Validate multichannel configuration
			void planeTestMultichannel(int mode) ;
			//Free space to test that the filtering procedure works
			void freeSpace();
			//This function shows how to load and use general static scenarios from files generated by Veneris 
			//It shows  the basic template for use of Opal
			void loadScenario();
			//This function also shows how to load and use general static scenarios from files generated by Veneris 
			//but now reading the configuration from a file
			void  loadScenarioLoraFromConfigFile();
			//Adding compund dynamic meshes. A number of meshes share a common transformation matrix
			//that can be used to move and rotate them around the scenario
			void addCompoundDynamicMeshes(); 
			//Use quad that is moved
			void quadTest( bool print, bool subSteps) ;
			//Simple update of receivers
			void moveReceivers();
			//How to remove dynamic meshes from the scenario
			void addRemoveDynamicMeshes( bool print, bool subSteps);
			//Test penetration in walls (transmission) in addition to reflection
			void penetrationTest(bool print, bool subSteps);
			void penetrationPlane( bool print, bool subSteps);
			//Street crossing test. Cubes are intended to be buildings and a plane is the floor
			void crossingTest( bool print, bool subSteps, bool useDepolarization, float radius);
			//Street crossing test again but now in a more efficient way.
			void crossingTestEfficient( bool print, bool subSteps, bool useDepolarization, float radius);
			//Crossing test now with a complex dynamic mesh  as vehicle
			void crossingTestAndVehicle(bool useDepolarization, float radius) ;
			//The crossing test again but now using multiple transmitters at the same time
			void crossingTestMulti( bool print, bool subSteps, bool useDepolarization, float radius);
			//How to use callbacks to get more information from Opal
			void useCallbacks(int mode) ;
			void exampleCallback(ResultRecord r, int txId);
			void anotherCallback(ResultRecord r, int txId);
			void hitCallback(int txId, int rxId, HitRecord r);
			//Utility functions
			std::vector<int> parseTestString(std::string test);
			std::vector<float3> loadReceiversFromFile(std::string file) ;
			std::vector<float3> loadReceiversFromJSON();

	};

}
//Basic tests

//Old versions. Uncomment if necessary to try
//std::unique_ptr<OpalSceneManager> addCompoundDynamicMeshes(std::unique_ptr<OpalSceneManager> sceneManager); 
////Adding, moving and removing dynamic meshes
//std::unique_ptr<OpalSceneManager> addRemoveDynamicMeshes(std::unique_ptr<OpalSceneManager> sceneManager, bool print, bool subSteps) ;
////Adding and removing dynamic meshes
//std::unique_ptr<OpalSceneManager> addRemoveReceivers(std::unique_ptr<OpalSceneManager> sceneManager); 
////moving receivers
//std::unique_ptr<OpalSceneManager> moveReceivers(std::unique_ptr<OpalSceneManager> sceneManager, bool useDepolarization) ;
////Street crossing test. Cubes are intended to be buildings and a plane is the floor
//std::unique_ptr<OpalSceneManager> crossingTest(std::unique_ptr<OpalSceneManager> sceneManager, bool print, bool subSteps, bool useDepolarization, float radius) ;
////Street crossing test. Cubes are intended to be buildings and a plane is the floor
//std::unique_ptr<OpalSceneManager> crossingTestMulti(std::unique_ptr<OpalSceneManager> sceneManager, bool print, bool subStepsi, bool useDepolarization, float radius) ;
////Two quads as walls and two overlapping receivers
//std::unique_ptr<OpalSceneManager> quadTest(std::unique_ptr<OpalSceneManager> sceneManager, bool print, bool subSteps) ;
////Street crossing with vehicle mesh 
//std::unique_ptr<OpalSceneManager> crossingTestAndVehicle(std::unique_ptr<OpalSceneManager> sceneManager) ;
////Penetration tests
////Penetration test. One cube, transmitter and receiver
//std::unique_ptr<OpalSceneManager> penetrationTest(std::unique_ptr<OpalSceneManager> sceneManager, bool print, bool subSteps) ;
////Penetration test. Plane 
//std::unique_ptr<OpalSceneManager> penetrationPlane(std::unique_ptr<OpalSceneManager> sceneManager, bool print, bool subSteps) ;
////Free space: to validate filtering
//std::unique_ptr<OpalSceneManager> freeSpaceRDN(std::unique_ptr<OpalSceneManager> sceneManager,  bool useDepolarization, float radius);
//std::unique_ptr<OpalSceneManager> freeSpace(std::unique_ptr<OpalSceneManager> sceneManager,  bool useDepolarization, float radius);
////Horizontal plane test. To validate against a two-ray model
//std::unique_ptr<OpalSceneManager> planeTest(std::unique_ptr<OpalSceneManager> sceneManager, float radius, bool useDepolarization );
//std::unique_ptr<OpalSceneManager> planeTestProgressive(std::unique_ptr<OpalSceneManager> sceneManager, float radius, bool useDepolarization) ;

#endif

