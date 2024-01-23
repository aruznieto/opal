/***************************************************************/
//
//Copyright (c) 2021 Esteban Egea-Lopez http://girtel.upct.es/~simulations
//
/**************************************************************/

#ifndef CURVEDFLATMESHSIMULATIONDOA_H
#define CURVEDFLATMESHSIMULATIONDOA_H
#include "../curvedFlatMeshSimulation.h"

namespace opal {
	//Simulation for scenarios with a mix of curved surfaces and  flat elements (walls) and arbitrary linear polarization
	class  DoALPCurvedFlatMeshReflectionSimulation : public LPCurvedFlatMeshReflectionSimulation {
		protected:
			virtual void  printHitInfo(HitInfo* host_hits, uint hits) override;
		public:
			DoALPCurvedFlatMeshReflectionSimulation(OpalSceneManager*  m);
			virtual std::string printConfigInfo() const override; 

	};


}


#endif
