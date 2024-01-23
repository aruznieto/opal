/***************************************************************/
//
//Copyright (c) 2021 Esteban Egea-Lopez http://girtel.upct.es/~simulations
//
/**************************************************************/
#ifndef CONFIGURATION_H
#define CONFIGURATION_H
#include <optix.h>
using namespace optix;
//Configuration
rtDeclareVariable(uint, usePenetration, , );
rtDeclareVariable(float, attenuationLimit, , );

rtDeclareVariable(uint, useMultichannel, , );
rtDeclareVariable(float, speedOfLight, , );
rtDeclareVariable(uint, useAntennaGain, , );

#endif
