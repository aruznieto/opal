
/***************************************************************/
//
//Copyright (c) 2019 Esteban Egea-Lopez http://girtel.upct.es/~simulations
//
/**************************************************************/

#ifndef SINGLEDIFFRACTION_H
#define SINGLEDIFFRACTION_H
#include "opalSimulation.h"

namespace opal {
	class  SingleDiffraction: public  OpalSimulation {
		protected:

			virtual void setOtherDirectory() override;
			virtual optix::Buffer setEdgeBuffer();
			virtual optix::Buffer setReceiverPositionsBuffer(uint nrx);
			virtual optix::Buffer setReceiverPolarizationBuffer(uint nrx);
			virtual optix::Buffer setHitBuffer(uint hits);
			virtual optix::Buffer setVisibilityBuffer(uint ne, uint nrx,uint ntx);
			virtual optix::Buffer setLoSBuffer(uint hits);

			virtual void fillEdgeBuffer();
			virtual void fillReceiverPositionsBuffer();
			virtual void fillReceiverPositionsBuffer(std::vector<SphereReceiver*>& rx);
			virtual void fillReceiverPolarizationBuffer();
			virtual void fillReceiverPolarizationBuffer(std::vector<SphereReceiver*>& rx);
			virtual void fillAntennaGainIdBuffer(); 
			virtual void fillAntennaGainIdBuffer(std::vector<SphereReceiver*>& rx); 
			optix::Buffer edgeBuffer;
			optix::Buffer receiverPositionsBuffer;
			optix::Buffer receiverPolarizationBuffer;
			optix::Buffer hitBuffer;
			optix::Buffer losBuffer;
			optix::Buffer visibilityBuffer;
			optix::Buffer antennaGainIdBuffer;
			optix::Buffer transformToPolarizationBuffer;
			bool updateReceiverBuffer;
			bool updateEdgeBuffer;
			unsigned int diffractionEntryIndex;
			unsigned int diffractionRayIndex;
			unsigned int visibilityEntryIndex;

			//Callback to work directly with hits recorded
			//std::vector<std::function<void(int,int,HitRecord)>> callbacks;
			void clearCallbacks();
			optix::Program  createVisibilityProgram();
			optix::Program  createComputeSimpleDiffractionProgram();
			optix::Program createMissDiffractionProgram() ;
			optix::Program createExceptionDiffractionProgram() ;
			optix::Program createClosestHitMeshDiffractionProgram(); 
			optix::Program createClosestHitCurvedMeshDiffractionProgram(); 
			virtual void createClosestHitPrograms() override;
			virtual void addReceiver(int id, float3  position, float3 polarization, float radius, std::function<void(float, int)>  callback, std::vector<optix::Material>& materials) override;
			//virtual void addReceiver(int id, float3  position, float3 polarization, float radius, std::function<void(float, int)>  callback) override;
			virtual void removeReceiver(int id) override;
			virtual void updateReceiver(int id, float3 position, float3 polarization, float radius) override;
			//virtual void processDiffractionLaunch(uint totalHits);
			virtual void processDiffractionLaunch(std::vector<uint3>& totalHits);
			virtual void saveTraceToFile(thrust::host_vector<LogTraceHitInfo> trace, std::string fileName) override;
			void processTraceLog(unsigned int maxTraceSize);
			void checkBuffersSize(uint nrx, uint ntx, uint ne);
			std::vector<uint3> processVisibilityBuffer();
			float maxMemoryFraction;
			bool adjustToMemorySize;
		       	
		public:	
			SingleDiffraction(OpalSceneManager*  m);
			virtual std::string printConfigInfo() const override;
			virtual void init() override; 
			//virtual void setDefaultPrograms(std::map<std::string,optix::Program>& defaultPrograms, optix::Material& defaultMeshMaterial) override;
			virtual void setInternalBuffers() override;
			virtual void executeTransmitLaunch(uint numTransmitters, bool partial) override;
			virtual void endPartialLaunch(uint numTransmitters) override;
			virtual void finishSceneContext() override;
			virtual std::string printInternalBuffersState() override;
			virtual void addStaticMesh(opal::OpalMesh& mesh, std::vector<optix::Material>& materials) override;
			virtual void addStaticCurvedMesh(opal::OpalMesh& mesh, std::vector<optix::Material>& materials) override ;
			virtual void registerReceiverGain(int rxId, int gainId) override; 
			virtual void transformEdge(Edge* e, optix::Matrix4x4 t) override;	
			virtual void clearReceivers();
		        //virtual void setCallback(std::function<void(int, int, HitRecord)> callback);	
			void setMaxMemoryFraction(float f) ;
			void setAdjustToMemorySize(bool adjust);
		        HitRecord toHitRecord(RDNHit* h );	

	};
}
#endif //SIMPLEDIFFRACTION_H
