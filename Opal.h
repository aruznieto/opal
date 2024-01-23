/***************************************************************/
//
//Copyright (c) 2019 Esteban Egea-Lopez http://girtel.upct.es/~simulations
//
/**************************************************************/

#ifndef OPAL_H
#define OPAL_H



//Do not change the order of includes, unless you want to fix many dependencies
#include "tutils.h"
#include "Common.h"
#include <map>
#include <fstream>
#include <sstream>
#include <functional>
#include <vector>
#include <memory>
#include <tuple>
#include <set>
#include "opalSimulation.h" 
#include "transmitterManager.h"
#include "raySphere.h"
#include <optix_world.h>
#include "results.h"
namespace opal {

 
 
	//Forward declarations
	class TransmitterManager;
	class OpalSimulation;
	class OpalRaySphereGenerator;
	class ResultRecord;
	class ResultReport;
	class HitCollector;
	class Exception : public std::exception {
		public:
			/// Create exception
			Exception(const std::string& message)
				: m_message(message) {}

			/// Virtual destructor (needed for virtual function calls inherited from
			/// std::exception).
			virtual ~Exception() throw() {}

			/// Retrieve the error message
			const std::string& getErrorString() const { return m_message; }

			/// From std::exception
			virtual const char* what() const throw() { return getErrorString().c_str(); }
		private:
			std::string m_message;

	};


	struct OpalMesh
	{
		optix::GeometryInstance      geom_instance;
		optix::float3                bbox_min;
		optix::float3                bbox_max;
		int                          num_triangles;
	};

	class OpalDynamicMeshGroup
	{
		public:
			optix::GeometryGroup		geom_group;
			optix::Transform	transform;
			unsigned int childIndex;
			std::vector<std::pair<Edge,Edge*>> edges;	
	};


	struct ChannelParameters {
		float frequency;
		float waveLength;
		float k;
		float eA;
		float c;
	};

	struct RaySphere {
		optix::Buffer raySphereBuffer;
		optix::uint elevationSteps;
		optix::uint azimuthSteps;
		optix::uint  rayCount;

	};

	struct EComponents {
			optix::float2 E;
			optix::float2 Ex;
			optix::float2 Ey;
			optix::float2 Ez;
			unsigned int index;
			int refHits;
			int difHits;
	};

	class SphereReceiver {
		public:
			optix::float3 position;
			optix::float3 polarization;
			float radius;
			optix::GeometryInstance geomInstance;
			//optix::Program closestHitProgram;
			std::function<void(float, int)> callback;
			std::vector<std::function<void(ResultRecord, int)>> extendedCallbacks;
			int externalId;
			bool dirty;
			std::vector<std::tuple<int,float3>> lastReceivedPower;
			std::vector<std::tuple<int,EComponents>> lastReceivedE;
			EComponents getLastReceivedE(int txId);  
			float getLastReceivedPower(int txId);  
			void setLastReceivedResult(int txId, float3 r);
			void setLastReceivedE(int txId, EComponents e);
			void clearLastResult();
			int antennaGainId;	
			optix::Matrix<4,4> transformToPolarization;
			
	};


	//Compare mesh faces. Note that for us a face is not just a mesh triangle. Each face is used to define an equal interaction element used for ray filtering.
	//For example a wall (made of many triangles). All rays hitting that same wall (equal face id) will be considered the same EM wave for filtering purposes
	//Differences below epsilon are consider equal to avoid precision problems
	struct face_compare
	{
		const float epsilon = 1e-6f;
                bool esentiallyEqual(float a,  float b) const {
			return fabs(a - b) <= ( (fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * epsilon);
		}
		bool operator() (const optix::float3 first, const optix::float3 second) const {
			//const float xa = fabs(first.x - second.x);
			//const float ya = fabs(first.y - second.y);
			//const float za = fabs(first.z - second.z);
			if (esentiallyEqual(first.x, second.x)  && esentiallyEqual(first.y, second.y) && esentiallyEqual(first.z, second.z)  ) {
			//if (xa <= epsilon && ya <= epsilon && za <= epsilon) {
				return false;
			}
			else {
				return (std::tie(first.x, first.y, first.z) < std::tie(second.x, second.y, second.z));
			}

		}
	};


	typedef std::map<optix::float3, unsigned int, face_compare> FaceMap; //Used to extract the faces of the meshes


	typedef std::vector<std::vector<float> > AntennaGain; //Azimuth/Elevation gain in linear units


	//System information
	
	struct DeviceInformation {
		char name[256];
		int computeCapability[2] = {0, 0};
		RTsize totalMemory = 0;
		int clockRate = 0;
		int maxThreadsPerBlock = 0;
		int smCount = 0;
		int executionTimeoutEnabled = 0;
		int maxHardwareTextureCount = 0 ;
		int tccDriver = 0;
		int cudaDeviceOrdinal = 0;
	};
	struct SystemInformation {
		unsigned int optixVersion;
		unsigned int major = optixVersion / 1000; // Check major with old formula.
		unsigned int minor;
		unsigned int micro;
		unsigned int numberOfDevices = 0;
		std::vector<DeviceInformation> devices;		
		void getSystemInformation();
		std::string printSystemInformation();
	};
	
	struct ConfigurationOptions {
		bool useExactSpeedOfLight; //speed of light (can be exact or approximated)
		bool usePenetration; //Allow penetration of rays through walls TODO: not thoroughly tested
		bool useDepolarization; //Allow general linear polarization
		bool useMultiGPU; //Use multiple GPUs if available
		bool useFastMath; //Enable CUDA fast math
		bool useMultiChannel; //Allow the use of different transmission frequencies
		bool useAntennaGain; //Apply radiation patterns to rays
		bool executeCallbacks; //Call callbacks on receivers
                bool printRecords; //Print individual hits
		bool recordHits; // Record individual hits, can only be used with certain types of simulations	
		//Improves performance if ray sphere changes from launch to launch. Otherwise create once a ray sphere
		bool generateRaysOnLaunch;
		std::tuple<bool,int,bool,optix::uint3> printEnabled; //Print debug CUDA info
	};
		
	class FieldInfo {
		protected:
			//External tx id to [external rx id, Ecomponents]
		std::map<int, std::map<int,EComponents>*>  Emap;
		//std::map<int, std::map<int,EComponents>> Ext;   // Complex
		public:
		void updateTransmitters(std::vector<Transmitter*>& activeTransmitters);
		void reset();
		void updateField(float2 E, optix::float2 Ex, optix::float2 Ey, optix::float2 Ez,  int rxId, int txId, unsigned int index,  int refHits, int difHits);
		
		std::map<int, std::map<int,EComponents>*> getE();	
	};	

	class OpalSceneManager {
		protected:
			//Source directory in SDK, used by sutil to find .cu programs
			std::string baseDir;
			std::string optixProgramsDir;
			std::string cudaProgramsDir;
			std::string currentDir;
			
			optix::Context context;


			optix::Buffer txOriginBuffer;
			optix::Buffer rayRangeBuffer;
			optix::Buffer raySphereParametersBuffer;
			
			//Gains
			std::map<int, optix::Buffer > antennaGainBuffers;
			//Scene graph variables
			optix::GeometryGroup receiversGroup;
			optix::Group rootGroup;
			optix::GeometryGroup staticMeshesGroup;
			std::vector<OpalMesh> staticMeshes;

			//Diffraction support
			std::vector<Edge*> edges;	


			//External to internal info
			std::map<int, unsigned int> receiverExtToBufferId; //External id to receptionInfo buffer Id
			std::vector<SphereReceiver*> receivers; //receivers info (mapped to receptionInfo buffer)
			
			//Transmitters
			TransmitterManager* transmitterManager;
		
			bool sceneGraphCreated;
			bool sceneFinished;
			bool contextInitialized;
		
			bool usingCurvedMeshes;	
			RaySphere raySphere;
			ChannelParameters defaultChannel;
			ChannelParameters currentChannelParameters;
			std::map<int, OpalDynamicMeshGroup*> dynamicMeshes; //externalId to dynamic mesh

			OpalRaySphereGenerator* rayGenerator;
			
			
			//Configuration parameters
			ConfigurationOptions options;
			
			unsigned int numberOfFaces;
			std::set<unsigned int> globalFaces;
			float attenuationLimit;	
			//PtxUtil*  ptxHandler;

			unsigned int maxReflections; //Default 10
			float minEpsilon; //Default 1.e-3f
			std::ostringstream configInfo;

#ifdef OPALDEBUG
			std::ofstream outputFile;
#endif // OPALDEBUG


			//Launch variables
			Transmitter currentTransmitter;
			
			
			//Default programs

			std::map<std::string,optix::Program> defaultPrograms;
			optix::Material defaultMeshMaterial;
			optix::Material defaultReceiverMaterial;
			optix::Material defaultCurvedMeshMaterial;
			//Default programs can be overriden in simulation classes
			void createDefaultPrograms();	
			std::vector<OpalSimulation*> simulations;
			float deg2rad;   //Degrees to radians
			PtxUtil*  ptxHandler;
			
			#ifdef OPAL_USE_TRI
			virtual optix::Program createTriangleAttributesProgram();
			virtual optix::Program createCurvedTriangleAttributesProgram();
			#endif

			virtual optix::Program createIntersectionTriangle();
			virtual optix::Program createBoundingBoxSphere();
			virtual optix::Program createIntersectionSphere();
			virtual optix::Program createBoundingBoxTriangle();
			virtual optix::Program createExceptionReflectionProgram();
					

	
			//System information
			SystemInformation sysInfo;
			std::vector<int> enabledDevices;
			
			//Accumulate results	
			FieldInfo* info;
			HitCollector* hitCollector;

		public:
		//Init
			OpalSceneManager();
			OpalSceneManager(float f,  bool useExactSpeedOfLight=true);
			virtual ~OpalSceneManager();
			virtual void initContext(float f);
			virtual void initMembers();
		//State
			ChannelParameters getChannelParameters() const;
			std::string getBaseDirectory() const;
			RaySphere getRaySphere() const; 
			optix::Context getContext() const;
			std::vector<int> getEnabledDevices() const;
			SystemInformation getSystemInformation() const;
			void setEnabledDevices();
			
			std::vector<Transmitter*> getActiveTransmitters();
			unsigned int getNumberOfActiveTransmitters() const;
			std::vector<SphereReceiver*> getReceivers();
			SphereReceiver* getReceiver(int externalId);
			Transmitter* getTransmitter(int externalId);
			unsigned int getTransmitterIndex(int externalId);
			unsigned int getReceiverIndex(int externalId);
			unsigned int getNumberOfReceivers() const;
			OpalSimulation* getSimulation() const;
			OpalSimulation* getSimulation(uint index) const;
			optix::Program getDefaultProgram(std::string p);	
			unsigned int getNumberOfEdges () const;
			void removeEdge(Edge* e);
			std::vector<Edge*> getEdges() ; 	
			Edge* addEdge(optix::float3 p, optix::float3 v, optix::uint2 faces, optix::float3 face_a, optix::float3 face_b, optix::float3 normal_a, optix::float3 normal_b, MaterialEMProperties emProp, int id=-1);
		//Static meshes	
			void setMeshEMProperties(optix::GeometryInstance geom_instance, MaterialEMProperties emProp);
			OpalMesh createMesh(int meshVertexCount, optix::float3* meshVertices, int meshTriangleCount, int* meshTriangles, optix::Program intersectionProgram, optix::Program boundingBoxProgram,  bool makeSingleFace = false);
			
			//Static mesh with default closest hit programs
			OpalMesh addStaticMesh(std::vector<optix::float3>& meshVertices, std::vector<int>& meshTriangles, optix::Matrix4x4 transformationMatrix, MaterialEMProperties emProp, bool makeSingleFace=false); 
			OpalMesh addStaticMesh(int meshVertexCount, optix::float3* meshVertices, int meshTriangleCount, int* meshTriangles, optix::Matrix4x4 transformationMatrix, MaterialEMProperties emProp, bool makeSingleFace = false);
			OpalMesh addStaticMeshWithFaces(std::vector<optix::float3> &meshVertices,std::vector<std::pair<optix::int3, unsigned int>> &triangleIndexFaceBuffer, optix::Matrix4x4 transformationMatrix, MaterialEMProperties emProp); 
			OpalMesh addStaticCurvedMesh(std::vector<optix::float3>& meshVertices, std::vector<int>& meshTriangles, std::vector<optix::float4>& pd1, std::vector<optix::float4>& pd2, optix::Matrix4x4 transformationMatrix, MaterialEMProperties emProp, bool makeSingleFace, int faceId = -1);
			void addStaticMesh(OpalMesh mesh);
			void setMeshFaceId(OpalMesh mesh, uint id);
			OpalMesh setMeshOnDevice(int meshVertexCount, optix::float3* meshVertices,std::vector<std::pair<optix::int3, unsigned int>> &triangleIndexBuffer, optix::Program intersectionProgram, optix::Program boundingBoxProgram); 
			MaterialEMProperties ITUparametersToMaterial(float a, float b, float c, float d, bool perfectConductor=false);

			bool checkFaceIds(std::vector<std::pair<optix::int3, unsigned int>> &triangleIndexFaceBuffer, int &uniqueFaces);
		//Moving meshes functions
			OpalDynamicMeshGroup* addDynamicMeshGroup(int groupId);
			void removeDynamicMeshGroup(int groupId);
			void  addMeshToGroup(int id, int meshVertexCount, optix::float3* meshVertices, int meshTriangleCount, int* meshTriangles,  MaterialEMProperties emProp);
			void addMeshWithFacesToGroup(int groupId, std::vector<optix::float3> &meshVertices,std::vector<std::pair<optix::int3, unsigned int>> &triangleIndexFaceBuffer,  MaterialEMProperties emProp); 
			void updateTransformInGroup(int groupId, optix::Matrix4x4 transformationMatrix);
			void finishDynamicMeshGroup(int groupId);
			const std::map<int, OpalDynamicMeshGroup*> &  getDynamicMeshes() const;
			void transformEdge(const Edge& original, Edge* e, optix::Matrix4x4 t); 
			Edge* addEdgeToGroup(optix::float3 p, optix::float3 v, optix::uint2 faces, optix::float3 face_a, optix::float3 face_b, optix::float3 normal_a, optix::float3 normal_b, MaterialEMProperties emProp, int id, int groupId); 
			
		//Ray spheres
		
			OpalRaySphereGenerator* getRaySphereGenerator() const;	
			void createRaySphereFromExternalBuffer(int elevationSteps, int azimuthSteps, optix::float3*  bpointer); 
			//Create and arbitrary 2D ray sphere, with ray directions provided by the user
			void createRaySphere2D(int elevationSteps, int azimutSteps, optix::float3 * rayDirections);
			//Create a 2D ray sphere in discrete steps of elevation and azimuth
			void createRaySphere2D(int elevationDelta, int azimuthDelta);
			//Create a ray sphere with fractions of degree
			//Now elevation delta and azimuthDelta are a decimal fraction of degree, that is, every unit is 0.1 degree, so elevationDelta=1 means 0.1 degree, elevationDelta=2 means 0.2 degree
			void createRaySphere2DSubstep(int elevationDelta, int azimuthDelta);
			//Create a 2D ray sphere with arbitray angular separation and solid angle coverage
			void createRaySphere2D(float initElevation, float elevationDelta, float endElevation, float initAzimuth, float azimuthDelta, float endAzimuth);
		
		//Antenna Gains
	// Consider our assumptions
	//- the orientation of the electric field vector on the ray when it is launched depends on the polarization of the antenna, since
	//  we are considering only linear polarizations and assume that the polarization is given by the antenna spatial orientation, which we actually call polarization,
	//  that is, we assume that a linear antenna (a dipole, for instance) oriented in Y axis (0,1,0) has a linear polarization=(0,1,0)
	//  If this antenna is rotated to (1,0,0), we assume the polarization is now (1,0,0) too.
	//- The antenna pattern file is assumed to be given with respect to a standard azimuth (Z is 0 degrees) and elevation (relative to Y). It may be independent of the polarization.
	// Imagine a directive antenna, let us say, with a beam pointing in the Z direction (0 azimuth). It may have linear polarization (0,1,0). If we orientate the antenna 
	//  90 degrees azimuth the beam is now pointing to X, but the linear polarization is still (0,1,0) it has not changed. We can tilt the beam -45 degrees to point to "ground" 
	//  and the polarization may be now (0,1,1) if we have physically rotated the antenna, or  still (0,1,0) if we have electrically rotated the beam 
	//  A dipole has a zero in the radiation pattern along the antenna orientation. If we rotate the antenna we rotate both the  polarization and the radiation pattern
	//  In Opal, we assume that if we change the polarization, we have physically rotated the antenna, and so the radiation pattern is also rotated accordingily,
	//  so when we use updateReceiver or transmit with a different polarization, the diagram is automatically rotated by Opal
	//  If one just want to rotate the radiation pattern without changing the polarization, one has to use 
	//    transmit(int txId, float txPower,  float3 origin, float3 polarization, float frequency, bool partial, optix::Matrix<4,4> antennaPatternTransform) ;
	//    and pass the transformation matrix for the diagram or externally orient the antenna (see antennaGain examples) and pass orientedGainAxes to true
	//
			
		//Register once a radiation pattern with Opal, it can be reused by multiple transmitters or receivers by associating it (see below) to them
		//Returns the gainId that is used by registerReceiverGain or registerTransmitterGain
			int registerAntennaGain(AntennaGain& gain);
		//Associate a given radiation pattern to transmitters or receivers 
			void registerReceiverGain(int rxId, int gainId); 
			void registerTransmitterGain(int txId, int gainId); 
			optix::Buffer getAntennaGainBuffer(int gainId);
		//Orientate the radiation pattern according to the antenna polarization (physical rotation)
			optix::Matrix<4,4> computeMatrixFromWorldToPolarization(float3 pol);
		//Orientate arbitrarily the radiation pattern	
			optix::Matrix<4,4> orientateAntennaPattern(float azimuthDegrees, float tilt);
			//optix::Matrix<4,4> orientateAntennaPattern(float azimuthDegrees, float tilt, float3 azimuthRotation, float3 zeroAzimuth);
			optix::Matrix<4,4> pointAntennaPatternTo(optix::float3 antennaPosition, optix::float3 point) ;


		//Receivers (actually, transceivers) a receiver is any point in space where we want to measure the
		//electromagnetic field. It can be a real antenna or just a location in space. Transmissions can be done from the location 
		//of the receiver, and in that case, the transmitter will not receive anything, if the id of trasnmitter is equal to the id of receiver
			void addReceiver(int id, optix::float3  position, optix::float3 polarization, float radius, std::function<void(float, int)>  callback, optix::Matrix<4,4>* antennaGainTransform=nullptr);
			void removeReceiver(int id);
			void updateReceiverAntennaTransform(int id, optix::Matrix<4,4> antennaGainTransform);
			void updateReceiver(int id, optix::float3 position);
			void updateReceiver(int id, optix::float3 position, optix::float3 polarization, float radius);
			void updateReceiver(int id, optix::float3 position, float radius);

			//Set function to be called  for detailed info of hits on a given receiver. Multiple callbacks can be set per receiver
			void setReceiverHitCallback(int id, std::function<void(ResultRecord,int)> extendedCallback);
			void clearReceiverHitCallbacks(int id);
			//Set function to be called   with the final power  on a given receiver. Only one function can be set per receiver
			void setReceiverCallback(int id, std::function<void(float,int)> callback);
			void clearReceivers();
		//Transmitters
			TransmitterManager* getTransmitterManager() ;
		//Transmit functions
			//Most general function so far, it establishes all the parameters, including the transform for the orientation of the antenna pattern (gain)
			ResultReport* transmit(int txId, float txPower,  float3 origin, float3 polarization, float frequency, bool partial, optix::Matrix<4,4> antennaPatternTransform) ;
			//As above, but antena pattern orientation may have been established externally and it may be set to false
			ResultReport* transmit(int txId, float txPower, optix::float3 origin, optix::float3 polarization, float frequency, bool partial=false, bool orientedGainAxes=false);
			//As above, but default frequency is used 
			ResultReport* transmit(int txId, float txPower, optix::float3 origin, optix::float3 polarization, bool partial=false, bool orientedGainAxes=false);

			//A utility function to transmit from a receiver (transceiver) position. It will take the missing
			//parameters (position, polarization, orientedGainAxes) from the receiver, which mya have been updated previously with updateReceiver
			ResultReport* transmitFromReceiver(int rxId, float txPower, float frequency);

			//Multitransmitter: transmit with a number of transmitters in the same launch, it is necessary to add the 
			//transmitters to the transmit group before
			ResultReport* groupTransmit(bool partial=false, bool orientedGainAxes=false);
			
	
		//Finish building scene
			void finishSceneContext();
			void printConfigInfo() const;
		//Configure launch
			ConfigurationOptions getConfigurationOptions() const;
			void setSimulation(OpalSimulation* sim);

			//Use antenna radiation patterns
			void setUseAntennaGain(bool use);
			//Call receiver callbacks when results are available
			void setExecuteCallback(bool execute);
			//Print all ray hits 
 			void setPrintRecords(bool printRecords);
			//Record and call calbacks for each ray hit on the sphere
			void setRecordHits(bool record);
			//Launch with multiple transmitters simultaneously
			void enableMultitransmitter();
			//Allow the use of different frequencies for different transmitters
			void enableMultiChannel();
			void disableMultiChannel();
			//Use multiple GPUs when available
			void enableMultiGPU();
			void disableMultiGPU();
			//Signals can propagate through walls (with attenuation). Does not work with RDN or curved surfaces
			//TODO: require further tests are carried out
			void enablePenetration();
			void disablePenetration();
			//Use linear polarization
			void enableDepolarization(); 
			void disableDepolarization();
			//Set default transmission frequency (when multichannel is not used)
			void setFrequency(float f);
			//Set speed of light to 3e8 m/s
			void useApproximateSpeedLight();
			//Set max attenuation for penetration in dB. Rays with a higher attenuation are not traced for penetration 
			void setAttenuationLimit(float f);
			//Set maximum number of reflections that will be traced
			void setMaxReflections(unsigned int m);
			unsigned int  getMaxReflections() const {return this->maxReflections;};
			//Set minimum separation distance from last geometry element hit by a ray
			//When rays hit an element and are reflected, their origin has to be separated from the element, 
			//Otherwise, due to numeric precision, they may hit again the same element
			void setMinEpsilon(float f);
			//Disabling fast math results in higher accuracy but worse performance
			void enableFastMath();
			void disableFastMath();
			//Set the ranges of the ray sphere
			void setRayRange(float initElevation, float elevationDelta,  float initAzimuth, float azimuthDelta, unsigned int elevationSteps,unsigned int azimuthSteps);
			void setRaySphereParameters(unsigned int elevationSteps, unsigned int azimuthSteps, unsigned int standardSphere);
			
			// To generate rays directly on the launch.  Improves performance when rays change from launch to launch. It is also necessary when the number of rays is high,
			// since no memory is allocated to them
			// Otherwise, create once a ray sphere. TODO: test whether using a ray sphere is better than using always generation on launch 
			void enableGenerateRaysOnLaunch();
			void disableGenerateRaysOnLaunch();
			
			void setBaseDir(std::string b);
			
			ResultReport* endPartialLaunch(uint numTransmitters);
			
			void enableExceptions();

			uint getNumberOfFaces();
			void setInitialHash(uint h);
			//void setMaxAngleForDuplicateRays(float angle); //In radians, all rays whose direction at departure is separated less than angle are considered duplicates when filtering
			std::map<std::string,optix::Program>& getDefaultPrograms();
			optix::Material getDefaultMeshMaterial() ;
			optix::Material getDefaultCurvedMeshMaterial() ;
			optix::Material getDefaultReceiverMaterial() ;


			


		//Log
			void setPrintEnabled(int bufferSize);
			void setPrintEnabled(int bufferSize, optix::uint3 index);
			void setUsageReport();
			virtual std::string printSceneReport();
			virtual std::string printContextInformation();
			virtual std::string printInternalBuffersState();
		//Util
			optix::float2 getAngles(optix::float3 const ray);
			//This should probably be somewhere else, but VS refuses to link it in other places and it is at so many places now that it is annoying to refactor, so we put it here to 
			std::vector<float3>  loadVerticesFromFile(const char* file); 
			std::vector<int>  loadTrianglesFromFile(const char* file) ;
			std::vector<float4>  loadPDFromFile(const char* file); 
			std::vector<float3>  loadRaysFromFile(const char* file);
			void writeMeshToPLYFile(std::string fileName, std::vector<optix::float3>& vertices, std::vector<int>& indices, optix::Matrix4x4& transformationMatrix); 
			void writeMeshToCustomFile(std::string fileName, std::vector<optix::float3>& vertices, std::vector<int>& indices, optix::Matrix4x4& transformationMatrix); 
			static void printPower(float power, int txId ); 
			AntennaGain loadGainsFromFileIndBPower(const char* file, bool linear=false);
			//Sign function
			template <typename T> int sgn(T val) {
			    return (T(0) < val) - (val < T(0));
			}; 
			float signedAngle(optix::float3 from, optix::float3 to, optix::float3 axis); 
			bool isUsingCurvedMeshes() const;
		//Results
			ResultReport* generateResultReport(bool callCallback);
			ResultRecord getFieldRecord(optix::float2 E, optix::float2 Ex, optix::float2 Ey, optix::float2 Ez, unsigned int index,unsigned int txIndex, int refHits, int difHits);
			FieldInfo* getFieldInfo() {return info;};
			float getReceivedPower(int rxId, int txId);
			EComponents getReceivedE(int rxId, int txId);


				
		protected:
		//Scene
			void extractFaces(optix::float3* meshVertices, std::vector<std::pair<optix::int3, unsigned int>> &triangleIndexBuffer);


			void createSceneContext();
			void buildSceneGraph();
			ResultReport* launchSingleTransmitter( bool partial, bool orientedGainAxes) ;



		//Internal buffers
			void setInternalBuffers();
			void setRayBuffers();
			optix::Buffer setTransmitterBuffer(optix::uint tx);
			optix::Buffer setRayRangeBuffer();
			optix::Buffer setRaySphereParametersBuffer();

		//Utils
			static void callbackUsageReport(int level, const char* tag, const char* msg, void* cbdata);
	};



} //namespace opal
#endif

