/***************************************************************/
//
//Copyright (c) 2021 Esteban Egea-Lopez http://girtel.upct.es/~simulations
//
/**************************************************************/

#include "flatSimulationDoa.h"
#include <iomanip>
namespace opal {
	DoALPFlatMeshReflectionSimulation::DoALPFlatMeshReflectionSimulation(OpalSceneManager* m) : LPFlatMeshReflectionSimulation(m) {
	}
	std::string DoALPFlatMeshReflectionSimulation::printConfigInfo() const  {
		std::ostringstream stream;
		stream<<"--- Simulation Type---"<<std::endl;
		stream<<"\tUsing DoALPFlatMeshReflectionSimulation "<<std::endl;
		stream<< "\t - Using depolarization. Transmitter and receiver can have any linear polarizations (linear antenna may be oriented in any direction) and can be different. DoA and unfoldedPathLength will be logged "<<std::endl;
		return stream.str();
	}
	void  DoALPFlatMeshReflectionSimulation::printHitInfo(HitInfo* host_hits, uint hits) {
		//This implementation assumes that hits come ordered by transmitter and receiver previously
		if (hits==0) {
			return;
		}
		//Get first transmitter 			
		uint currentTx=host_hits->thrd.x;
		//Get first receiver
		uint index=host_hits->thrd.z;
		std::vector<Transmitter*> activeTransmitters = myManager->getActiveTransmitters();
		std::vector<SphereReceiver*>	receivers=myManager->getReceivers(); 
		for (uint i=0; i<hits; i++) {
			currentTx=host_hits->thrd.x; 				
			//New receiver,  start new accumulation 				
			index=host_hits->thrd.z; 				
			//std::cout<<i<<"\t refhash="<<(host_hits)->thrd.y<<std::endl;
			float4 doad=host_hits->doaD;
			float4 dod=host_hits->rdud;
			Transmitter* tx=activeTransmitters[currentTx];
			if (mode==ComputeMode::FIELD) {
				const float2 Ex=make_float2(host_hits->EEx.z,host_hits->EEx.w) ;
				const float2 Ey=make_float2(host_hits->EyEz.x,host_hits->EyEz.y) ;
				const float2 Ez=make_float2(host_hits->EyEz.z,host_hits->EyEz.w) ;
				//std::cout<<i<<"\tEx="<<make_float2(host_hits->EEx.z,host_hits->EEx.w) <<std::endl;
				//std::cout<<i<<"\tEy="<<make_float2(host_hits->EyEz.x,host_hits->EyEz.y) <<std::endl;
				//std::cout<<i<<"\tEz="<<make_float2(host_hits->EyEz.z,host_hits->EyEz.w) <<std::endl;
				std::cout<<std::setprecision(15)<<"DOA\t"<<i<<"\t"<<Ex.x<<"\t"<<Ex.y<<"\t"<<Ey.x<<"\t"<<Ey.y<<"\t"<<Ez.x<<"\t"<<Ez.y<<"\t"<<doad.x<<"\t"<<doad.y<<"\t"<<doad.z<<"\t"<<doad.w<<"\t"<<receivers[index]->externalId<<"\t"<< tx->externalId<<"\t"<<dod.x<<"\t"<<dod.y<<"\t"<<dod.z<<std::endl;
				//DEBUG REFLECTIONS
				//std::cout<<std::setprecision(15)<<"DOA\t"<<i<<"\t"<<host_hits->info.x<<"\t"<<Ex.x<<"\t"<<Ex.y<<"\t"<<Ey.x<<"\t"<<Ey.y<<"\t"<<Ez.x<<"\t"<<Ez.y<<"\t"<<doad.x<<"\t"<<doad.y<<"\t"<<doad.z<<"\t"<<doad.w<<"\t"<<receivers[index]->externalId<<"\t"<< tx->externalId<<"\t"<<dod.x<<"\t"<<dod.y<<"\t"<<dod.z<<std::endl;
			} else {
				//std::cout<<i<<"\tE="<<(host_hits)->E<<std::endl;
				float2 E=make_float2(host_hits->EEx.x,host_hits->EEx.y);	
			//Print ray number, E, DOA, unfolded distance, rxId, txId, DOD
				std::cout<<std::setprecision(15)<<"DOA\t"<<i<<"\t"<<E.x<<"\t"<<E.y<<"\t"<<doad.x<<"\t"<<doad.y<<"\t"<<doad.z<<"\t"<<doad.w<<"\t"<<receivers[index]->externalId<<"\t"<< tx->externalId<<"\t"<<dod.x<<"\t"<<dod.y<<"\t"<<dod.z<<std::endl;
				//std::cout<<std::setprecision(15)<<"DOA\t"<<i<<"\t"<<E.x<<"\t"<<E.y<<"\t"<<doad.x<<"\t"<<doad.y<<"\t"<<doad.z<<"\t"<<doad.w<<"\t"<<receivers[index]->externalId<<"\t"<<receivers[index]->position.x<<"\t"<<receivers[index]->position.y<<"\t"<<receivers[index]->position.z<<"\t"<< tx->externalId<<"\t"<<tx->origin_p.x<<"\t"<<tx->origin_p.y<<"\t"<<tx->origin_p.z<<std::endl;
	//			std::cout<<std::setprecision(15)<<"DOD\t"<<i<<"\t"<<dod.x<<"\t"<<dod.y<<"\t"<<dod.z<<"\t"<<dod.w<<"\t"<<receivers[index]->externalId<<"\t"<<receivers[index]->position.x<<"\t"<<receivers[index]->position.y<<"\t"<<receivers[index]->position.z<<"\t"<< tx->externalId<<"\t"<<tx->origin_p.x<<"\t"<<tx->origin_p.y<<"\t"<<tx->origin_p.z<<std::endl;
			}
			//std::cout<<i<<"\t dist="<<(host_hits)->rdud.w<<std::endl;
			//std::cout<<std::setprecision(15)<<i<<"\t dir="<<(host_hits)->rdud<<std::endl;
			++host_hits;
		}
	}
}
