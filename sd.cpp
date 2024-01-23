/***************************************************************/
//
//Copyright (c) 2019 Esteban Egea-Lopez http://girtel.upct.es/~simulations
//
/**************************************************************/

#include "singleDiffraction.h"
#include <iomanip>
#include <optixu/optixu_matrix_namespace.h>
#include "timer.h"
namespace opal {
	SingleDiffraction::SingleDiffraction(OpalSceneManager* m) : OpalSimulation(m) {
		this->withPolarization=true;
		this->withRaySphere = false;
		this->acceptEdge = true;
		this->acceptCurved =true; //We accept curved just to allow curved meshes to be used in the scenario. Curved diffraction is not computed
		this->simType =OpalSimulationTypes::SINGLEDIFFRACTION;
		this->updateReceiverBuffer = true;
		this->updateEdgeBuffer = true;
		this->edgeBuffer=nullptr;
		this->hitBuffer = nullptr;
		this->transformToPolarizationBuffer = nullptr;
		this->receiverPositionsBuffer=nullptr;
		this->receiverPolarizationBuffer=nullptr;
		this->mode=ComputeMode::VOLTAGE;
		this->diffractionRayIndex=0u;
		this->diffractionEntryIndex=0u;
		this->resultsMessage="SingleDiffraction";
		//this->callback=nullptr;
		this->adjustToMemorySize=false;
		this->maxMemoryFraction=0.5f;
	}
	void SingleDiffraction::setAdjustToMemorySize(bool adjust) {
		this->adjustToMemorySize=adjust;
	}
	void SingleDiffraction::setMaxMemoryFraction(float f) {
		this->maxMemoryFraction=f;
	}
	std::string SingleDiffraction::printConfigInfo() const  {
		std::ostringstream stream;
		stream<<"--- Simulation Type---"<<std::endl;
		stream<<"\tUsing SingleDiffraction"<<std::endl;
	        if (!isEnabled()) {
			stream<<"\t SingleDiffraction DISABLED"<<std::endl;
		}	
		return stream.str();
	}
	std::string SingleDiffraction::printInternalBuffersState() {
		std::ostringstream stream;
	 	RTsize w,h,d;
		unsigned long long totalBytes = 0;
		unsigned long long sb;
		if (hitBuffer) {
			hitBuffer->getSize(w,h,d);
			sb = sizeof(RDNHit)*w*h*d;
			totalBytes += sb;
			stream << "\t hitBuffer=(" << w <<"," <<h<<","<<d<<"). size=" << (sb / (1024.f*1024.f)) << " MiB" << std::endl;
		}
		if (receiverPositionsBuffer) {
			receiverPositionsBuffer->getSize(w);
			sb = sizeof(float4)*w;
			totalBytes += sb;
			stream << "\t receiverPositionsBuffer=(" << w <<"). size=" << (sb / (1024.f*1024.f)) << " MiB" << std::endl;
		}
		if (edgeBuffer) {
			edgeBuffer->getSize(w);
			sb = sizeof(Edge)*w;
			totalBytes += sb;
			stream << "\t edgeBuffer=(" << w <<"). size=" << (sb / (1024.f*1024.f)) << " MiB" << std::endl;
		}
		//Check memory usage
		stream << "Total memory in SingleDiffraction internal buffers:  " << (totalBytes / (1024.f*1024.f)) << " MiB" << std::endl;
		return stream.str();
	}
	void SingleDiffraction::setOtherDirectory()  {
		this->currentDir=(cudaProgramsDir+"/polarization/diffraction");
	}
	void SingleDiffraction::fillEdgeBuffer() {
		//Fill the GPU edge buffer
		std::vector<Edge*> edges=myManager->getEdges();
		Edge* edges_gpu=reinterpret_cast<Edge*>(edgeBuffer->map());
		for (size_t i=0; i<edges.size(); ++i) {
			edges_gpu[i]=(*edges[i]);
			//std::cout<<"ede["<<i<<"]<<"<<edges_gpu[i].id<<std::endl;
		}
		edgeBuffer->unmap();	
	}

	optix::Buffer SingleDiffraction::setHitBuffer(uint ne, uint nrx, uint ntx) {
		if (hitBuffer) {
			hitBuffer->destroy();
		}
		//std::cout<<"Setting difBuffer to ("<<ne<<","<<nrx<<","<<ntx<<")"<<std::endl;
		optix::Buffer b = myManager->getContext()->createBuffer(RT_BUFFER_OUTPUT, RT_FORMAT_USER, ne,nrx,ntx);
		//Reuse RDN hit form
		b->setElementSize(sizeof(RDNHit));
		myManager->getContext()["difBuffer"]->set(b);
		return b;
	}
	optix::Buffer SingleDiffraction::setEdgeBuffer() {
		if (edgeBuffer) {
			edgeBuffer->destroy();
		}
		uint nEdges=myManager->getNumberOfEdges();
		if (nEdges==0) {
			nEdges=1u;
		}
		optix::Buffer b = myManager->getContext()->createBuffer(RT_BUFFER_INPUT, RT_FORMAT_USER, nEdges);
		b->setElementSize(sizeof(Edge));
		myManager->getContext()["edgeBuffer"]->set(b);
		return b;
	}
	optix::Buffer SingleDiffraction::setReceiverPositionsBuffer(uint nrx) {
		if (receiverPositionsBuffer) {
			receiverPositionsBuffer->destroy();
		}
		optix::Buffer b = myManager->getContext()->createBuffer(RT_BUFFER_INPUT, RT_FORMAT_FLOAT4, nrx);
		myManager->getContext()["receiverPositionsBuffer"]->set(b);
		return b;
	}
	optix::Buffer SingleDiffraction::setReceiverPolarizationBuffer(uint nrx) {
		if (receiverPolarizationBuffer) {
			receiverPolarizationBuffer->destroy();
		}
		optix::Buffer b = myManager->getContext()->createBuffer(RT_BUFFER_INPUT, RT_FORMAT_FLOAT3, nrx);
		myManager->getContext()["receiverPolarizationBuffer"]->set(b);
		return b;
	}
	void SingleDiffraction::fillAntennaGainIdBuffer(std::vector<SphereReceiver*>&  rx) {
		antennaGainIdBuffer->setSize(static_cast<unsigned int>(rx.size()));
		transformToPolarizationBuffer->setSize(static_cast<unsigned int>(rx.size()));
		int* b = static_cast<int*>(antennaGainIdBuffer->map());
		optix::Matrix<4,4>* m=static_cast<optix::Matrix<4,4>*>(transformToPolarizationBuffer->map());
		for (int i=0; i<rx.size(); i++) {
			int gainId = rx[i]->antennaGainId;
			if (gainId>=0) {
				optix::Buffer buffer=myManager->getAntennaGainBuffer(gainId);
				b[i]=buffer->getId();	
				optix::Matrix<4,4> pol_t=myManager->computeMatrixFromWorldToPolarization(rx[i]->polarization);
				m[i]=pol_t;
				
			} else {
				b[i]=RT_BUFFER_ID_NULL;	
				m[i]=optix::Matrix<4,4>::identity();
			}
		}
		antennaGainIdBuffer->unmap();
		transformToPolarizationBuffer->unmap();
	}
	void SingleDiffraction::fillAntennaGainIdBuffer() {
		std::vector<SphereReceiver*>  rx=myManager->getReceivers();
		fillAntennaGainIdBuffer(rx);
	}
	void SingleDiffraction::transformEdge(Edge* e, optix::Matrix4x4 t) {
		updateEdgeBuffer=true;
	}	
	void SingleDiffraction::fillReceiverPositionsBuffer() {
		std::vector<SphereReceiver*>  rx=myManager->getReceivers();
		fillReceiverPositionsBuffer(rx);
	}
	void SingleDiffraction::fillReceiverPolarizationBuffer() {
		std::vector<SphereReceiver*>  rx=myManager->getReceivers();
		fillReceiverPolarizationBuffer(rx);
	}
	void SingleDiffraction::fillReceiverPositionsBuffer(std::vector<SphereReceiver*>& rx) {
		float4* rx_gpu=nullptr;
		rx_gpu=reinterpret_cast<float4*>(receiverPositionsBuffer->map());
		for (int i=0; i<rx.size(); i++) {
			//float4 sphere=make_float4(rx[i]->position,rx[i]->radius);
			//We do not really need the sphere radius, but we need the id of the receiver, so we use the last position for id
			float4 sphere=make_float4(rx[i]->position,rx[i]->externalId);
			rx_gpu[i]=sphere;
		}
		receiverPositionsBuffer->unmap();	
	}
	void SingleDiffraction::fillReceiverPolarizationBuffer(std::vector<SphereReceiver*>& rx) {
		float3* rx_gpu=nullptr;
		rx_gpu=reinterpret_cast<float3*>(receiverPolarizationBuffer->map());
		for (int i=0; i<rx.size(); i++) {
			//float4 sphere=make_float4(rx[i]->position,rx[i]->radius);
			//We do not really need the sphere radius, but we need the id of the receiver, so we use the last position for id
			//float4 sphere=make_float4(rx[i]->position,rx[i]->externalId);
			rx_gpu[i]=rx[i]->position;
		}
		receiverPolarizationBuffer->unmap();	
	}
	void SingleDiffraction::setInternalBuffers() {
		uint nrx=myManager->getNumberOfReceivers();
		if (nrx==0) {
			nrx=1u;
		}
		uint ne=myManager->getNumberOfEdges();
		if (ne==0) {
			ne=1u;
		}
		uint ntx=myManager->getNumberOfActiveTransmitters();
		if (ntx==0) {
			ntx=1u;
		}
		edgeBuffer=setEdgeBuffer();
		if (adjustToMemorySize) {
			std::cout<<"Adjust to memory set: Delaying buffer setup"<<std::endl;
			nrx=1u;
			receiverPositionsBuffer=setReceiverPositionsBuffer(nrx);
			receiverPolarizationBuffer=setReceiverPolarizationBuffer(nrx);
		} else {
			receiverPositionsBuffer=setReceiverPositionsBuffer(nrx);
			receiverPolarizationBuffer=setReceiverPolarizationBuffer(nrx);
			fillEdgeBuffer();
		}
		hitBuffer=setHitBuffer(ne, nrx, ntx);
		//Buffers to be resized later
		traceBuffer = myManager->getContext()->createBuffer(RT_BUFFER_INPUT_OUTPUT | RT_BUFFER_GPU_LOCAL, RT_FORMAT_USER, 1u);
		traceBuffer->setElementSize(sizeof(LogTraceHitInfo));
		myManager->getContext()["traceBufferDiffraction"]->set(traceBuffer);
		traceAtomicIndexBuffer = myManager->getContext()->createBuffer(RT_BUFFER_INPUT_OUTPUT | RT_BUFFER_GPU_LOCAL, RT_FORMAT_UNSIGNED_INT, 1u);
		myManager->getContext()["traceAtomicIndexDiffraction"]->set(traceAtomicIndexBuffer);
		antennaGainIdBuffer = myManager->getContext()->createBuffer(RT_BUFFER_INPUT, RT_FORMAT_BUFFER_ID,1u);
		myManager->getContext()["antennaGainIdBuffer"]->set(antennaGainIdBuffer);	
		transformToPolarizationBuffer = myManager->getContext()->createBuffer(RT_BUFFER_INPUT, RT_FORMAT_USER,1u);
		transformToPolarizationBuffer->setElementSize(sizeof(optix::Matrix<4,4>));
		myManager->getContext()["transformToPolarizationBuffer"]->set(transformToPolarizationBuffer);
		if (generateTraceLog) {
			myManager->getContext()["traceDiffraction"]->setUint(1u);
		} else {
			myManager->getContext()["traceDiffraction"]->setUint(0u);
		}	
	}
	void SingleDiffraction::setReducedBuffersSize(std::vector<SphereReceiver*> rx, uint ntx, uint ne) {
		RTsize w,h,d;
		uint nrx=rx.size();
		edgeBuffer->getSize(w);
		if (w!=ne) {
			edgeBuffer=setEdgeBuffer();
			fillEdgeBuffer();
			updateEdgeBuffer=false;
		}
		if (updateEdgeBuffer) {
			fillEdgeBuffer();
			updateEdgeBuffer=false;
		}
		hitBuffer=setHitBuffer(ne, nrx, ntx);
		receiverPositionsBuffer=setReceiverPositionsBuffer(nrx);
		receiverPolarizationBuffer=setReceiverPolarizationBuffer(nrx);
		fillReceiverPositionsBuffer(rx);
		fillReceiverPolarizationBuffer(rx);
		if (myManager->getConfigurationOptions().useAntennaGain) {
				fillAntennaGainIdBuffer(rx);
		}
		updateReceiverBuffer=false;
	
		
	}
	void SingleDiffraction::checkBuffersSize(uint nrx, uint ntx, uint ne) {
		RTsize w,h,d;
		edgeBuffer->getSize(w);
		if (w!=ne) {
			edgeBuffer=setEdgeBuffer();
			fillEdgeBuffer();
			updateEdgeBuffer=false;
		}
		if (updateEdgeBuffer) {
			fillEdgeBuffer();
			updateEdgeBuffer=false;
		}
		hitBuffer->getSize(w,h,d);
		if ((w!=ne) || (h != nrx) || (d != ntx)) {
			hitBuffer=setHitBuffer(ne, nrx, ntx);
		}
		
		receiverPositionsBuffer->getSize(w);
		if (updateReceiverBuffer || (w!=nrx)) {
			receiverPositionsBuffer=setReceiverPositionsBuffer(nrx);
			receiverPolarizationBuffer=setReceiverPolarizationBuffer(nrx);
			fillReceiverPositionsBuffer();
			fillReceiverPolarizationBuffer();
			if (myManager->getConfigurationOptions().useAntennaGain) {
				fillAntennaGainIdBuffer();
			}
			updateReceiverBuffer=false;
		}
		if (generateTraceLog) {
			RTsize bf=2*ne*nrx*ntx;
			traceBuffer->setSize(bf);
		}
		
	}
	void SingleDiffraction::executeTransmitLaunch(uint numTransmitters, bool partial) {
		if (partial) {
			//TODO:Partial launches are not well defined so far for diffraction. We just compute all the edges in one launch, do not discriminate by zones
			//this is something that can be done in the future. At the moment we just do not do nothing and delay until the end of partial launch
			std::cout<<"Partial launch for diffraction. Delaying computation until end of partial launch"<<std::endl;
			return;
		}
		uint ne=myManager->getNumberOfEdges();
		uint nrx=myManager->getNumberOfReceivers();
		if (ne==0) {
			std::cout<<"No edges found on scene. Not running single diffraction simulation" <<std::endl;
			return;
			
		
		}
		Timer launchTime;
		launchTime.start();

		if (adjustToMemorySize) {
			if (generateTraceLog) {
				std::cout<<"Trace log not created for separated launches"<<std::endl;
			}
			separateLaunches(numTransmitters,ne);
		} else {
			checkBuffersSize(nrx, numTransmitters, ne);
			//Diffraction launch
			std::cout<<"Diffraction. Launch ["<<ne<<","<<nrx<<","<<numTransmitters<<"]"<<std::endl;
			std::cout<<printInternalBuffersState()<<std::endl;
			checkLaunchSize(ne,nrx,numTransmitters);
			myManager->getContext()->launch(diffractionEntryIndex, ne,nrx,numTransmitters); //Launch 3D (edges, receivers, transmitters);
			processDiffractionLaunch();
			if (generateTraceLog) {
				uint maxTraceSize=2*ne*nrx*numTransmitters;
				processTraceLog(maxTraceSize);
			}
		}
		launchTime.stop();
		std::cout<<"Diffraction launch time=" <<launchTime.getTime()<<std::endl;
		
	}
	void SingleDiffraction::separateLaunches(uint ntx, uint ne) {
		SystemInformation si = myManager->getSystemInformation();
		unsigned long long totalMemory=0;
		for (auto d: si.devices) {
			 totalMemory+=(unsigned long long) d.totalMemory;
		}
		std::vector<SphereReceiver*>  rx=myManager->getReceivers();
		uint totalRx=rx.size();


		unsigned long long requiredMemory = sizeof(RDNHit)*ne*ntx*totalRx;
		uint maxNrx=floor((maxMemoryFraction*totalMemory)/(ne*ntx*sizeof(RDNHit)));

		uint batches = ceil(totalRx/maxNrx);
		std::cout<<"Required memory for diffraction launch="<< requiredMemory<<"; totalMemory="<<totalMemory<<"; max fraction="<<maxMemoryFraction*totalMemory<<";max number of receivers="<<maxNrx<<"batches="<<batches<<std::endl;
		if (batches==0) {
			checkBuffersSize(totalRx, ntx, ne);
			//Diffraction launch
			std::cout<<"Diffraction. Launch ["<<ne<<","<<totalRx<<","<<ntx<<"]"<<std::endl;
			checkLaunchSize(ne,totalRx,ntx);
			myManager->getContext()->launch(diffractionEntryIndex, ne,totalRx,ntx); //Launch 3D (edges, receivers, transmitters);
			processDiffractionLaunch();
			return;
		}
		for (int b=0; b<=batches; ++b) {
			//std::vector<SphereReceiver*>::const_iterator first=rx.begin()+b*maxNrx;
			//std::vector<SphereReceiver*>::const_iterator last=rx.begin()+(b+1)*maxNrx;
			std::vector<SphereReceiver*> receivers;
			std::cout<<"Batch="<<b<<std::endl;
			for (int i=0; i<maxNrx; ++i) {
				uint index=(b*maxNrx)+i;
				if (index>(totalRx-1)){
					break;
				}
				receivers.push_back(rx[index]);
			}

			setReducedBuffersSize(receivers, ntx, ne);
			//Diffraction launch
			std::cout<<"Reduced diffraction. Launch ["<<ne<<","<<receivers.size()<<","<<ntx<<"]"<<std::endl;
			checkLaunchSize(ne,receivers.size(),ntx);
			myManager->getContext()->launch(diffractionEntryIndex, ne,receivers.size(),ntx); //Launch 3D (edges, receivers, transmitters);
			processDiffractionLaunch(receivers);

		}
	

	}
	void SingleDiffraction::processDiffractionLaunch() {
		std::vector<SphereReceiver*>	receivers=myManager->getReceivers();
	        processDiffractionLaunch(receivers);	
	}
		
	void SingleDiffraction::processDiffractionLaunch(std::vector<SphereReceiver*>& receivers) {
		RTsize w,he,d;
		hitBuffer->getSize(w,he,d);	
		RDNHit* h=reinterpret_cast<RDNHit*>(hitBuffer->map());
		//uint nrx=myManager->getNumberOfReceivers();
		uint nrx=receivers.size();
		std::vector<Transmitter*> activeTransmitters = myManager->getActiveTransmitters();
		FieldInfo* info=myManager->getFieldInfo();
		std::cout<<"processDiffractionLaunch() nrx="<<nrx<<std::endl;
		int realHits=0;
		for (unsigned int x = 0; x < w; ++x) // Edges 
		{ 
			for (unsigned int y = 0; y < he; ++y) //Receivers
			{
				for (unsigned int z = 0; z < d; ++z) //transmitters
				{
					unsigned int i=(z*w*he)+(y*w)+x;
					if (mode==ComputeMode::VOLTAGE) {
						float2 E=make_float2(h[i].EEx.x,h[i].EEx.y);
						//std::cout<<"processDiffractionLaunch() E("<<x<<","<<receivers[y]->externalId<<","<<activeTransmitters[z]->externalId<<")="<<E<<std::endl;
					        	
						if (h[i].EEx.z==1.0f) {
							//std::cout<<"["<<x<<","<<y<<","<<z<<"]="<<E<<std::endl;
							//computeReceivedPower(E,y,z, 1u);
							for (auto cb: callbacks) {	
								HitRecord hr=toHitRecord(&h[i]);
								cb( activeTransmitters[z]->externalId,receivers[y]->externalId, hr);
								//cb( activeTransmitters[z]->externalId,receivers[y]->externalId, &h[i]);
							}	
							++realHits;
							if (printHits) {
								float4 doad=h[i].doaD;
								float4 dod=h[i].doDu;
							//DOA
								std::cout<<std::setprecision(15)<<"DIFD\t"<<E.x<<"\t"<<E.y<<"\t"<<doad.x<<"\t"<<doad.y<<"\t"<<doad.z<<"\t"<<doad.w<<"\t"<<receivers[y]->externalId<<"\t"<< activeTransmitters[z]->externalId<<"\t"<<dod.x<<"\t"<<dod.y<<"\t"<<dod.z<<std::endl;
							
								//std::cout<<std::setprecision(15)<<"DIFD\t"<<E.x<<"\t"<<E.y<<"\t"<<receivers[y]->externalId<<"\t"<< activeTransmitters[z]->externalId<<"\t"<<doad.w<<"\t"<<dod.w<<"\t"<<doad.x<<"\t"<<doad.y<<"\t"<<doad.z<<std::endl;
								//std::cout<<std::setprecision(15)<<"DIFD\t"<<E.x<<"\t"<<E.y<<"\t"<<receivers[y]->externalId<<"\t"<< x<<"\t"<<doad.w<<"\t"<<dod.w<<"\t"<<doad.x<<"\t"<<doad.y<<"\t"<<doad.z<<std::endl;
								//std::cout<<std::setprecision(15)<<"DIFD\t"<<E.x<<"\t"<<E.y<<"\t"<< x<<"\t"<<y<<"\t"<<z<<std::endl;
							}
						}
						if (info) {
							info->updateField(E,receivers[y]->externalId,activeTransmitters[z]->externalId,y,0,h[i].EEx.z); 			
						}
					} else {
						float2 Ex=make_float2(h[i].EEx.z,h[i].EEx.w);
						float2 Ey=make_float2(h[i].EyEz.x,h[i].EyEz.y);
						float2 Ez=make_float2(h[i].EyEz.z,h[i].EyEz.w);
						if (h[i].EEx.x==1.0f) {
							if (printHits) {
								float4 doad=h[i].doaD;
								float4 dod=h[i].doDu;
								std::cout<<std::setprecision(15)<<"DIFD\t"<<Ex.x<<"\t"<<Ex.y<<"\t"<<Ey.x<<"\t"<<Ey.y<<"\t"<<Ez.x<<"\t"<<Ez.y<<"\t"<<doad.x<<"\t"<<doad.y<<"\t"<<doad.z<<"\t"<<doad.w<<"\t"<<receivers[y]->externalId<<"\t"<< activeTransmitters[z]->externalId<<"\t"<<dod.x<<"\t"<<dod.y<<"\t"<<dod.z<<std::endl;
							}
							for (auto cb: callbacks) {
								HitRecord hr=toHitRecord(&h[i]);
								cb( activeTransmitters[z]->externalId,receivers[y]->externalId, hr);
								//cb( activeTransmitters[z]->externalId,receivers[y]->externalId, &h[i]);
							}	
							//std::cout<<"["<<x<<","<<y<<","<<z<<"]"<<std::endl;
							//computeReceivedPower(Ex,Ey,Ez,y, z, 1u); 
							++realHits;
						}
						if (info) {	
							info->updateField(Ex,Ey,Ez,receivers[y]->externalId,activeTransmitters[z]->externalId,y,0, h[i].EEx.x); 			
						}
					}
				}
			}
		}
		hitBuffer->unmap(); 
		std::cout<<"\t Diffraction hits="<<realHits<<std::endl;
	}
	HitRecord SingleDiffraction::toHitRecord(RDNHit* h) {
		HitRecord hr;
		hr.diffraction=true;
		hr.E=make_float2(h->EEx.x,h->EEx.y);
		hr.Ex=make_float2(h->EEx.z,h->EEx.w);
		hr.Ey=make_float2(h->EyEz.x,h->EyEz.y);
		hr.Ez=make_float2(h->EyEz.z,h->EyEz.w);
		hr.directionOfArrival=make_float3(h->doaD.x,h->doaD.y,h->doaD.z);
		hr.directionOfDeparture=make_float3(h->doaD.x,h->doaD.y,h->doaD.z);
		hr.unfoldedPath=h->doaD.w;
		
		
	}
	void SingleDiffraction::processTraceLog(unsigned int maxTraceSize) {
		std::cout<<"Processing Diffraction Trace Log with maxTraceSize="<<maxTraceSize<<std::endl;
		std::vector<int> enabledDevices= myManager->getEnabledDevices();
		thrust::host_vector<LogTraceHitInfo> trace=opalthrustutils::getLogTraceOrderer(traceBuffer,traceAtomicIndexBuffer,enabledDevices, maxTraceSize);
		if (trace.size()>0) {
			saveTraceToFile(trace,"dif-trace.txt");
		} else {
			std::cout<<"No trace for diffraction generated="<<maxTraceSize<<std::endl;
		}	
	}
	void SingleDiffraction::saveTraceToFile(thrust::host_vector<LogTraceHitInfo> trace, std::string fileName) {
		std::cout<<"Saving diffraction log trace to "<<fileName<<std::endl;
		std::ofstream file(fileName.c_str(),std::ofstream::out);
		uint currentIndex=0;
		for (int i=0; i<trace.size(); i +=2) {
			LogTraceHitInfo l=trace[i];
			LogTraceHitInfo o=trace[i+1];
			file<<currentIndex<<":"<<l.hitp.x<<"\t"<<l.hitp.y<<"\t"<<l.hitp.z;
			file<<"|"<<o.hitp.x<<"\t"<<o.hitp.y<<"\t"<<o.hitp.z;	
			file<<std::endl;
			++currentIndex;
		}
		file.close();

	}
	void SingleDiffraction::updateReceiver(int id, float3 position, float3 polarization, float radius) {
		updateReceiverBuffer=true;
	}
	void SingleDiffraction::removeReceiver(int id) {
		updateReceiverBuffer=true;
	}
	void SingleDiffraction::addReceiver(int id, float3  position, float3 polarization, float radius, std::function<void(float, int)> callback, std::vector<optix::Material>& materials) {
		updateReceiverBuffer=true;
		//std::cout<<"SingleDiffraction::addReceiver() "<<std::endl;
	}
	void SingleDiffraction::registerReceiverGain(int rxId, int gainId) {
		updateReceiverBuffer=true;
	}
        void SingleDiffraction::clearReceivers() {
		updateReceiverBuffer=true;
	}
	void SingleDiffraction::endPartialLaunch(uint numTransmitters) {
		//Compute now delayed diffraction launc
		std::cout<<"Called endPartialLaunch for diffraction. Executing diffraction launch now"<<std::endl;
		executeTransmitLaunch(numTransmitters, false);
	}
	//void SingleDiffraction::setDefaultPrograms() {
	////void SingleDiffraction::setDefaultPrograms(std::map<std::string,optix::Program>& defaultPrograms, optix::Material& defaultMeshMaterial) {
	//	//Create compiler options here
	//	std::vector<const char *> nvccOptions;

	//	nvccOptions.push_back("-arch");
	//	nvccOptions.push_back("compute_30");
	//	ConfigurationOptions options=myManager->getConfigurationOptions(); 
	//	if (options.useFastMath) {
	//		nvccOptions.push_back("-use_fast_math");
	//	}
	//	nvccOptions.push_back("-lineinfo");
	//	nvccOptions.push_back("-default-device");
	//	nvccOptions.push_back("-rdc");
	//	nvccOptions.push_back("true");
	//	nvccOptions.push_back("-D__x86_64");

	//	ptxHandler = new PtxUtil(nvccOptions);
	//	
	//	setOtherDirectory();
	//	
	//	defaultPrograms.insert(std::pair<std::string, optix::Program>("computeSimpleDiffraction",createComputeSimpleDiffractionProgram()));
	//	defaultPrograms.insert(std::pair<std::string, optix::Program>("meshClosestHitDiffraction", createClosestHitMeshDiffractionProgram()));
	//	defaultPrograms.insert(std::pair<std::string, optix::Program>("missDiffraction", createMissDiffractionProgram()));
	//	defaultMeshMaterial->setClosestHitProgram(diffractionEntryIndex, defaultPrograms.at("meshClosestHitDiffraction")); //Add a program for visibility rays for diffraction
	//}
	void SingleDiffraction::createClosestHitPrograms() {
		//std::cout<<"SingleDiffraction::createClosestHitPrograms()"<<std::endl;
		std::map<std::string, optix::Program>& defaultPrograms=myManager->getDefaultPrograms();
		defaultPrograms.insert(std::pair<std::string, optix::Program>("computeSimpleDiffraction",createComputeSimpleDiffractionProgram()));
		//std::cout<<"SingleDiffraction::createClosestHitPrograms() meshClosestHitDiffraction"<<std::endl;
		defaultPrograms.insert(std::pair<std::string, optix::Program>("meshClosestHitDiffraction", createClosestHitMeshDiffractionProgram()));
		//std::cout<<"SingleDiffraction::createClosestHitPrograms() meshClosestHitCurvedDiffraction"<<std::endl;
		defaultPrograms.insert(std::pair<std::string, optix::Program>("meshClosestHitCurvedDiffraction", createClosestHitCurvedMeshDiffractionProgram())); //Only to be able to hit on curved meshes
		defaultPrograms.insert(std::pair<std::string, optix::Program>("missDiffraction", createMissDiffractionProgram()));
		defaultPrograms.insert(std::pair<std::string, optix::Program>("exceptionDiffraction", createExceptionDiffractionProgram()));
	}
	void SingleDiffraction::addStaticMesh(OpalMesh& mesh, std::vector<optix::Material>& materials) {
		std::map<std::string, optix::Program>& defaultPrograms=myManager->getDefaultPrograms();
		//std::cout<<"SingleDiffraction::addStaticMesh(): Setting closest hit for mesh with diffractionRayIndex="<<diffractionRayIndex<<std::endl;
		materials[0]->setClosestHitProgram(diffractionRayIndex,defaultPrograms.at("meshClosestHitDiffraction")); //Visibility rays for diffraction
	}
	void SingleDiffraction::addStaticCurvedMesh(OpalMesh& mesh, std::vector<optix::Material>& materials) {
		std::map<std::string, optix::Program>& defaultPrograms=myManager->getDefaultPrograms();
		std::cout<<"SingleDiffraction::addStaticCurvedMesh(): Setting closest hit for mesh with diffractionRayIndex="<<diffractionRayIndex<<std::endl;
		materials[0]->setClosestHitProgram(diffractionRayIndex,defaultPrograms.at("meshClosestHitCurvedDiffraction")); //Visibility rays for diffraction
	}
	void SingleDiffraction::init() {
		unsigned int entryPointCount= myManager->getContext()->getEntryPointCount();
		//Add this ray and entry type
		diffractionEntryIndex = entryPointCount; 
		myManager->getContext()->setEntryPointCount(diffractionEntryIndex+1);
		diffractionRayIndex= myManager->getContext()->getRayTypeCount();
		myManager->getContext()->setRayTypeCount(diffractionRayIndex+1u);
		std::cout<<"SingleDiffraction::init() diffractionEntryIndex="<<diffractionEntryIndex<<"diffractionRayIndex="<<diffractionRayIndex<<std::endl;
		if (mode==ComputeMode::VOLTAGE) {
			//myManager->getContext()["computeMode"]->setUint(0u);
			checkAndSetComputeMode(0u);
		} else {
			//myManager->getContext()["computeMode"]->setUint(1u);
			checkAndSetComputeMode(1u);
		}
	}
	void SingleDiffraction::finishSceneContext() {
		myManager->getContext()->setRayGenerationProgram(diffractionEntryIndex, myManager->getDefaultProgram("computeSimpleDiffraction")); //Diffraction
		//Miss programs require ray_type_index instead of entry_point_index
		myManager->getContext()->setMissProgram(diffractionRayIndex, myManager->getDefaultProgram("missDiffraction"));
		myManager->getContext()->setExceptionProgram(diffractionEntryIndex, myManager->getDefaultProgram("exceptionDiffraction"));
	}
	optix::Program  SingleDiffraction::createComputeSimpleDiffractionProgram()
	{

		optix::Program prog= myManager->getContext()->createProgramFromPTXString(ptxHandler->getPtxString(currentDir.c_str(), "singleDiffraction.cu"), "computeSingleDiffraction");
		//Set the corresponding ray type
		prog["rayTypeIndex"]->setUint(diffractionRayIndex);
		return prog;

	}	
	optix::Program SingleDiffraction::createMissDiffractionProgram() 
	{
		return myManager->getContext()->createProgramFromPTXString(ptxHandler->getPtxString(currentDir.c_str(), "singleDiffraction.cu"), "missDiffraction");

	}
	optix::Program SingleDiffraction::createClosestHitMeshDiffractionProgram() {
		return myManager->getContext()->createProgramFromPTXString(ptxHandler->getPtxString(currentDir.c_str(), "singleDiffraction.cu"), "closestHitTriangleDiffraction");
	}
	optix::Program SingleDiffraction::createClosestHitCurvedMeshDiffractionProgram() {
		return myManager->getContext()->createProgramFromPTXString(ptxHandler->getPtxString(currentDir.c_str(), "singleDiffraction.cu"), "closestHitCurvedTriangleDiffraction");
	}

	optix::Program SingleDiffraction::createExceptionDiffractionProgram() 
	{
		return myManager->getContext()->createProgramFromPTXString(ptxHandler->getPtxString(currentDir.c_str(), "singleDiffraction.cu"), "exception");

	}

	void SingleDiffraction::clearCallbacks() {
		this->callbacks.clear();
	}
}

