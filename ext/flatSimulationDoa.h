/***************************************************************/
//
//Copyright (c) 2021 Esteban Egea-Lopez http://girtel.upct.es/~simulations
//
/**************************************************************/

#ifndef FLATSIMULATIONDOA_H
#define FLATSIMULATIONDOA_H
#include "../flatSimulation.h"

namespace opal {



	//Simulation for scenarios with only flat elements (walls) and arbitrary linear polarization
	class  DoALPFlatMeshReflectionSimulation : public  LPFlatMeshReflectionSimulation {
		protected:
			virtual void  printHitInfo(HitInfo* host_hits, uint hits) override;
		public:	
			DoALPFlatMeshReflectionSimulation(OpalSceneManager*  m);
			virtual std::string printConfigInfo() const override; 
	};
}
#endif


