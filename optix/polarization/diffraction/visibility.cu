/***************************************************************/
//
//Copyright (c) 2019 Esteban Egea-Lopez http://girtel.upct.es/~simulations
//
/**************************************************************/
#include "../../../Common.h"
#include <optix_world.h>
#include "../linearPolarizationFunctions.h"
#include "../../configuration.h"
#include "../../reflectionFunctions.h"
#include "../../receiverFunctions.h"
#include "diffractionFunctions.h"
using namespace optix;

//Constants
#define D_2_Pi 6.283185307179586f
#define D_SQRT_2_Pi 2.506628274631f
#define D_SQRT_2ByPi 0.797884560802865

//Scene root
//rtDeclareVariable(rtObject, root, , ); //Already defined in trace functions
//Static meshes root
//rtDeclareVariable(rtObject, staticMeshesRoot, , );

//Edge buffer
rtBuffer<Edge, 1> edgeBuffer;

//TODO:These  should be redundant if we used Receiver buffers
//Receiver position buffer 
rtBuffer<float4, 1> receiverPositionsBuffer;
//rtBuffer<float3, 1> receiverPolarizationBuffer;
//rtBuffer<rtBufferId<float,2>, 1> antennaGainIdBuffer;
//typedef optix::Matrix<4,4> TransMat; 
//rtBuffer<TransMat, 1> transformToPolarizationBuffer;


//Visibility buffer
rtBuffer<uint, 3> visibilityBuffer; //Buffer to store all the hits
//Hit buffer
//rtBuffer<RDNHit, 3> difBuffer; //Buffer to store all the hits

//Local variables
//rtDeclareVariable(float, k, , ); //wavenumber 2pi/lambda
//rtDeclareVariable(uint, computeMode, ,);
//rtDeclareVariable(uint, traceDiffraction, ,);

//Transmitter buffer
rtBuffer<Transmitter, 1> txBuffer;

//Launch variables
rtDeclareVariable(uint3, launchIndex, rtLaunchIndex, );

//Visibility ray payload
rtDeclareVariable(VisibilityPayload, rayPayload, rtPayload, );
rtDeclareVariable(uint, rayTypeIndex, , );

//For debug only
rtDeclareVariable(TriangleHit, ch_triangle_data, attribute triangle_hit_data, );
rtDeclareVariable(CurvedTriangleHit, curved_triangle_data, attribute curved_triangle_hit_data, );
rtDeclareVariable(optix::Ray, ray_hit, rtCurrentRay, );


//Diffraction Launch program
RT_PROGRAM void computeDiffractionVisibility() {

	//3D launch [edges,receivers,transmitters]
		//rtPrintf("%u\t%u\t%u Launch \n",launchIndex.x,launchIndex.y,launchIndex.z);
	uint3 difBufferIndex=launchIndex;
//	//Initialize buffer to make sure it does not carry values from previous launches
//	RDNHit aHit;
//	aHit.EEx=make_float4(0.0f,0.0f, 0.0f,0.0f);
//	aHit.EyEz=make_float4(0.0f,0.0f,0.0f,0.0f);
//	difBuffer[difBufferIndex]=aHit;
//
	visibilityBuffer[difBufferIndex]=0u;
	Transmitter tx = txBuffer[launchIndex.z];

	const float3 origin = make_float3(tx.origin_p);
	const float4 sphere = receiverPositionsBuffer[launchIndex.y];
	//Check if ray is hitting his own tx (transmitter are also receivers usually) A transmitter cannot receive while it is transmitting, unless other channel is used.
	if (static_cast<int>(sphere.w)==tx.externalId) {
		return;	
	}	
	const float3 destination =make_float3(sphere.x,sphere.y,sphere.z);
	Edge e = edgeBuffer[launchIndex.x];
	if (!isDiffractingEdge(origin,e)) {
		//rtPrintf("%u\t%u\t%u Not diffracting edge %d\n",launchIndex.x,launchIndex.y,launchIndex.z,e.id);
		return;
	}

	//Compute diffraction point (DP) between transmitter, receiver and edge
	float3 dp; 
	if (computeDiffractionPoint(origin,destination,e,dp)){
			//rtPrintf("%u\t%u\t%u e=%u dp=(%f,%f,%f) \n",launchIndex.x,launchIndex.y,launchIndex.z,e.id,dp.x,dp.y,dp.z);
		VisibilityPayload visibilityRayPayload;
		visibilityRayPayload.polarization_k = tx.polarization_k; 
		visibilityRayPayload.result.x=OPAL_DIFFRACTION_LOS;
		visibilityRayPayload.faces=e.faces; //Faces of the edge, where the ray can hit to compute the diffraction. If it hits any other face, there is no LoS
		visibilityRayPayload.result.y=0;
		////trace visibility from transmitter to DP
		float3 originToDP=dp-origin;
		float dist_originToDp=length(originToDP);
		float3 txRayDirection = originToDP/dist_originToDp;
		//optix::Ray visibilityRay(origin,txRayDirection , rayTypeIndex, 0.0f,dist_originToDp); //Visibility ray type = 1
		optix::Ray visibilityRay(origin,txRayDirection , rayTypeIndex, min_t_epsilon,dist_originToDp-min_t_epsilon); //Visibility ray type = 1
			//rtPrintf("%u\t%u\t%u e=%u dp=(%f,%f,%f) tx ray=(%f,%f,%f) \n",launchIndex.x,launchIndex.y,launchIndex.z,e.id, dp.x,dp.y,dp.z,visibilityRay.direction.x,visibilityRay.direction.y,visibilityRay.direction.z);
		//TODO: Only check visibility with static meshes so far. Change if we want to consider  moving meshes (such as vehicles)
		//WARNING: ONLY THIS METHODS WORKS. ANY OF THE ONE BELOW GIVES WRONG RESULTS, IT MAY BE A OPTIX BU
		rtTrace(root, visibilityRay, visibilityRayPayload,OPAL_STATIC_MESH_MASK,RT_RAY_FLAG_DISABLE_ANYHIT);

		//rtTrace(root, visibilityRay, visibilityRayPayload,RT_VISIBILITY_ALL);
		//rtTrace(root, visibilityRay, visibilityRayPayload,RT_VISIBILITY_ALL,RT_RAY_FLAG_DISABLE_ANYHIT);
		//rtTrace(root, visibilityRay, visibilityRayPayload);
		if (visibilityRayPayload.result.x!=OPAL_DIFFRACTION_BLOCKED) {
				//rtPrintf("%u\t%u\t%u e=%u dp=(%f,%f,%f) tx not blocked \n",launchIndex.x,launchIndex.y,launchIndex.z,e.id, dp.x,dp.y,dp.z);
			//trace visibility from receiver to DP
			float3 destinationToDP=dp-destination;
				//rtPrintf("%u\t%u\t%u e=%u destinationToDP=(%f,%f,%f) tx not blocked \n",launchIndex.x,launchIndex.y,launchIndex.z,e.id, destinationToDP.x,destinationToDP.y,destinationToDP.z);
			float dist_destinationToDp=length(destinationToDP);
			visibilityRay.origin=destination;
			visibilityRay.direction=destinationToDP/dist_destinationToDp;
			visibilityRay.tmin=0.0f;
			visibilityRay.tmax=dist_destinationToDp; 
			float3 rxRayDir=destinationToDP/dist_destinationToDp;
			//optix::Ray visibilityRayRx(destination, rxRayDir , rayTypeIndex, 0.0f,dist_destinationToDp); //Visibility ray type = 1
			optix::Ray visibilityRayRx(destination, rxRayDir , rayTypeIndex, min_t_epsilon,dist_destinationToDp-min_t_epsilon); //Visibility ray type = 1

			VisibilityPayload visibilityRayPayloadRx;
			visibilityRayPayloadRx.polarization_k = tx.polarization_k; 
			visibilityRayPayloadRx.result.x=OPAL_DIFFRACTION_LOS;
			visibilityRayPayloadRx.faces=e.faces; //Faces of the edge, where the ray can hit to compute the diffraction. If it hits any other face, there is no LoS
			visibilityRayPayloadRx.result.y=1;
			//visibilityRayPayload.result.x=OPAL_DIFFRACTION_LOS;
			//visibilityRayPayload.result.y=1;
				//rtPrintf("%u\t%u\t%u e=%u dp=(%f,%f,%f) rx ray=(%f,%f,%f) d=%f \n",launchIndex.x,launchIndex.y,launchIndex.z,e.id, dp.x,dp.y,dp.z,visibilityRay.direction.x,visibilityRay.direction.y,visibilityRay.direction.z, dist_destinationToDp	);
			float3 rxRay=visibilityRay.direction;
				//rtPrintf("%u\t%u\t%u e=%u rxRay=(%f,%f,%f) rx not blocked \n",launchIndex.x,launchIndex.y,launchIndex.z,e.id, rxRay.x,rxRay.y,rxRay.z);
				//rtPrintf("%u\t%u\t%u e=%u diff sangles(beta, beta',phi, phi')=(%f,%f,%f,%f) spar(s,s')=(%f,%f) face_0=(%f,%f,%f)  \n",launchIndex.x,launchIndex.y,launchIndex.z,e.id, sangles.x,sangles.y,sangles.z,sangles.w,spar.x,spar.y, face_0.x,face_0.y,face_0.z);

			rtTrace(root, visibilityRayRx, visibilityRayPayloadRx,OPAL_STATIC_MESH_MASK,RT_RAY_FLAG_DISABLE_ANYHIT);
			//rtTrace(root, visibilityRayRx, visibilityRayPayload,RT_VISIBILITY_ALL,RT_RAY_FLAG_DISABLE_ANYHIT);
			//rtTrace(root, visibilityRay, visibilityRayPayload,RT_VISIBILITY_ALL);
			//rtTrace(root, visibilityRayRx, visibilityRayPayloadRx);
			if (visibilityRayPayloadRx.result.x!=OPAL_DIFFRACTION_BLOCKED) {
				//LoS, so compute diffraction
				float2 spar; //Distance parameters [s,s']
				float3 n_iplane; //Normal vector of the incidence plane
				float3 n_dplane; //Normal vector of the diffraction plane
				float4 R_0; 
				float4 R_n;
				float4 angles=getDiffractionParameters<VisibilityPayload>(visibilityRayPayload,origin, destination,e,dp,spar, n_iplane, n_dplane, R_0, R_n);
				float4 sangles = angles*180.0f/M_PIf;
				
					//rtPrintf("%u\t%u\t%u e=%u  sangles(beta, beta',phi, phi')=(%f,%f,%f,%f) spar(s,s')=(%f,%f)   \n",launchIndex.x,(launchIndex.y),launchIndex.z,e.id, sangles.x,sangles.y,sangles.z,sangles.w,spar.x, spar.y);
				if (angles.z>=(M_PIf*e.pn.w)) {
					//Receiver is between face_0 and face_n (inside the edge). It cannot receive even if there is no blocking (actually blocking may not be detected by visibility)
						//rtPrintf("%u\t%u\t%u e=%u receiver is inside the wedge  sangles(beta, beta',phi, phi')=(%f,%f,%f,%f) spar(s,s')=(%f,%f)   \n",launchIndex.x,launchIndex.y,launchIndex.z,e.id, sangles.x,sangles.y,sangles.z,sangles.w,spar.x);
					return;
				}
				visibilityBuffer[difBufferIndex]=1;
			} 
		}
	} else {
		//rtPrintf("%u\t%u\t%u Not diffracting point on edge %d\n",launchIndex.x,launchIndex.y,launchIndex.z,e.id);
	}
}
//Closest hit program for triangles and visibility rays
//It is recommended to use CH instead of AH since AH forces to use the SM, see https://developer.download.nvidia.com/video/gputechconf/gtc/2019/presentation/s9768-new-features-in-optix-6.pdf
RT_PROGRAM void closestHitTriangleDiffraction() {
	//float3 hp=ray_hit.origin + ray_hit.direction*ch_triangle_data.geom_normal_t.w;
	//float3 hp=ray_hit.origin + ray_hit.direction*rtIntersectionDistance();
	//	rtPrintf("%u\t%u\t%u hit on face=%u hp=(%f,%f,%f) \n",launchIndex.x,launchIndex.y,launchIndex.z,ch_triangle_data.faceId,hp.x,hp.y,hp.z);
	//rtPrintf("%u\t%u\t%u\t%d   hit on face=%u hp=(%f,%f,%f) \n",launchIndex.x,launchIndex.y,launchIndex.z,rayPayload.result.y,ch_triangle_data.faceId,hp.x,hp.y,hp.z);
	if (ch_triangle_data.faceId!=rayPayload.faces.x && ch_triangle_data.faceId!=rayPayload.faces.y) { 
		//rtPrintf("%u\t%u\t%u\t%d blocked on face=%u ray.faces=%u,%u hp=(%f,%f,%f) \n",launchIndex.x,launchIndex.y,launchIndex.z,rayPayload.result.y,ch_triangle_data.faceId,rayPayload.faces.x, rayPayload.faces.y,hp.x,hp.y,hp.z);
		//LoS is blocked, set flag 
		rayPayload.result.x=OPAL_DIFFRACTION_BLOCKED;
	} else {
		//rtPrintf("%u\t%u\t%u\t%d LoS on face=%u ray.faces=%u,%u hp=(%f,%f,%f) \n",launchIndex.x,launchIndex.y,launchIndex.z,rayPayload.result.y,ch_triangle_data.faceId,rayPayload.faces.x, rayPayload.faces.y,hp.x,hp.y,hp.z);
	}

}
RT_PROGRAM void closestHitCurvedTriangleDiffraction() {
	//float3 hp=ray_hit.origin + ray_hit.direction*curved_triangle_data.geom_normal_t.w;
	if (curved_triangle_data.faceId!=rayPayload.faces.x && curved_triangle_data.faceId!=rayPayload.faces.y) { 
		//rtPrintf("%u\t%u\t%u\t%d blocked on curved face=%u ray.faces=%u,%u hp=(%f,%f,%f) \n",launchIndex.x,launchIndex.y,launchIndex.z,rayPayload.result.y,curved_triangle_data.faceId,rayPayload.faces.x, rayPayload.faces.y,hp.x,hp.y,hp.z);
		//LoS is blocked, set flag 
		rayPayload.result.x=OPAL_DIFFRACTION_BLOCKED;
	} else {
		//rtPrintf("%u\t%u\t%u\t%d LoS on curved face=%u ray.faces=%u,%u hp=(%f,%f,%f) \n",launchIndex.x,launchIndex.y,launchIndex.z,rayPayload.result.y,curved_triangle_data.faceId,rayPayload.faces.x, rayPayload.faces.y,hp.x,hp.y,hp.z);
	}

}
//Miss program for visibility rays
RT_PROGRAM void missDiffraction() {
	//rtPrintf("%u\t%u\t%u miss \n", launchIndex.x,launchIndex.y,launchIndex.z); 
	//rtPrintf("%u\t%u\t%u\t%d miss \n", launchIndex.x,launchIndex.y,launchIndex.z,rayPayload.result.y); 
	rayPayload.result.x=OPAL_DIFFRACTION_MISS;
}

RT_PROGRAM void exception()
{
	const unsigned int code = rtGetExceptionCode();
	if (RT_EXCEPTION_USER <= code)
	{
		printf("Diffraction visibility user exception %d at (%d, %d)\n", code - RT_EXCEPTION_USER, launchIndex.x, launchIndex.y);
	}
	else
	{
		printf("Diffraction visibility Exception code 0x%X at (%d, %d)\n", code, launchIndex.x, launchIndex.y);
	}

}

