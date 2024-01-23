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


//Scene root
//rtDeclareVariable(rtObject, root, , ); //Already defined in trace functions
//Static meshes root
//rtDeclareVariable(rtObject, staticMeshesRoot, , );

//Edge buffer
rtBuffer<Edge, 1> edgeBuffer;

//TODO:These  should be redundant if we used Receiver buffers
//Receiver position buffer 
rtBuffer<float4, 1> receiverPositionsBuffer;
rtBuffer<float3, 1> receiverPolarizationBuffer;
rtBuffer<rtBufferId<float,2>, 1> antennaGainIdBuffer;
typedef optix::Matrix<4,4> TransMat; 
rtBuffer<TransMat, 1> transformToPolarizationBuffer;

//Hit buffer
rtBuffer<RDNHit, 1> difBuffer; //Buffer to store all the hits
//Real hits buffer
rtBuffer<uint3, 1> losBuffer; 

//Local variables
//rtDeclareVariable(float, k, , ); //wavenumber 2pi/lambda
rtDeclareVariable(uint, computeMode, ,);
rtDeclareVariable(uint, traceDiffraction, ,);

//Transmitter buffer
rtBuffer<Transmitter, 1> txBuffer;

//Launch variables
rtDeclareVariable(uint, launchIndex, rtLaunchIndex, );

//Visibility ray payload
rtDeclareVariable(VisibilityPayload, rayPayload, rtPayload, );
rtDeclareVariable(uint, rayTypeIndex, , );

//For debug only
rtDeclareVariable(TriangleHit, ch_triangle_data, attribute triangle_hit_data, );
rtDeclareVariable(CurvedTriangleHit, curved_triangle_data, attribute curved_triangle_hit_data, );
rtDeclareVariable(optix::Ray, ray_hit, rtCurrentRay, );


rtBuffer<LogTraceHitInfo, 1> traceBufferDiffraction;
rtBuffer<uint, 1> traceAtomicIndexDiffraction; //Buffer to store the current log trace buffer index 



__forceinline__ __device__ void computeElectricFieldAtReceiver(const float4& angles, const float2& spar, const float3& txRay, const float4& polarization_k, const float3& n_iplane, const float3& n_dplane, float n, float3& rxRay, float4& R_0, float4& R_n , const float gain ) {
	uint3 index=losBuffer[launchIndex];
	const uint edgeIndex=	index.x;
	const uint rxIndex = index.y;
	const uint txIndex = index.z;
	const float k = polarization_k.w;
	const float3 pol=make_float3(polarization_k);
	const float s=spar.x;
	const float s_prime=spar.y;
	//Parameter A for spherical waves, see Balanis TODO: add support for other waves in the future
	const float A=sqrtf(s_prime/(s*(s_prime+s))); 
	const float sinbetap=sinf(angles.y);
	//L parameter for spherical and conical incidence, see Balanis. TODO: add support for other waves in the future 
	const float L=(s*s_prime*sinbetap*sinbetap)/(s_prime+s);
	//Get the electrical field vector for this ray
	const float3 Ev=getLinearPolarizationForRaySimple(pol,txRay); //It is normalized

		//rtPrintf("%u\tray=(%f,%f,%f) pol=(%f,%f,%f) Ev=(%f,%f,%f) A=%f L=%f)\n",launchIndex.y,txRay.x,txRay.y,txRay.z,pol.x,pol.y,pol.z,Ev.x,Ev.y,Ev.z,A,L);
		//rtPrintf("%u\tray=(%f,%f,%f)  Ev=(%f,%f,%f) A=%f L=%f s=%f,s'=%f)\n",launchIndex.x,txRay.x,txRay.y,txRay.z,Ev.x,Ev.y,Ev.z,A,L,s,s_prime);
		//rtPrintf("rt=(%f,%f,%f) pol=(%f,%f,%f) Ev2=(%f,%f,%f) Evsimple2=(%f,%f,%f)\n",rt.x,rt.y,rt.z,pol.x,pol.y,pol.z,Ev3.x,Ev3.y,Ev3.z,Ev4.x,Ev4.y,Ev4.z);

	//Get the unit vector parallel to the plane of incidence
	const float3 phi_p = n_iplane;
	const float3 beta_p=normalize(cross(txRay,phi_p)); 
	const float3 ss=normalize(cross(phi_p,beta_p)); 
		//rtPrintf("%u\t beta_p=(%f,%f,%f) phi_p=(%f,%f,%f) s_p=(%f,%f,%f) \n",launchIndex.x,beta_p.x,beta_p.y,beta_p.z,phi_p.x,phi_p.y,phi_p.z,ss.x,ss.y,ss.z);

	//Compute incident electric field at the point of diffraction (complex)	
	float2 z = make_float2(0.0f, -k*s_prime);
	float2 zexp = complex_exp_only_imaginary(z);
	const float2 Ei = sca_complex_prod((gain/s_prime), zexp); //TODO: Assuming the initial amplitude is 1.0. To Change if antenna radiation patter or similar is used 

	//Geometric parts
	const float Ei_beta=dot(beta_p,Ev); //component of the incident E field parallel to the plane of incidence at the point of diffraction
	const float Ei_phi=dot(phi_p,Ev); //component of the incident E field perpendicular to the plane of incidence at the point of diffraction



	//Incident complex amplitude components
	float2 Ei_beta_q=sca_complex_prod(Ei_beta,Ei);
	float2 Ei_phi_q=sca_complex_prod(Ei_phi,Ei);
		//rtPrintf("%u\t Ev=(%f,%f,%f) Ei_beta=%f Ei_phi=%f Ei_beta_q=(%f,%f) Ei_phi_q=(%f,%f) an=%f\n",launchIndex.y,Ev.x,Ev.y,Ev.z,Ei_beta,Ei_phi, Ei_beta_q.x,Ei_beta_q.y, Ei_phi_q.x,Ei_phi_q.y,(atan2(Ei_phi_q.y,Ei_phi_q.x)*180/M_PIf));


	//Diffraction coefficients
	float4 D=computeLuebbersDiffractionCoefficient(k,n, angles.z,angles.w,angles.y,L, R_0, R_n);


	//Test: store diffraction coefficient
	//		RDNHit aHit;
	//		aHit.EEx=make_float4(0.0f,0.0f,0.0f,0.0f);
	//		aHit.EyEz=D;
	//		difBuffer[index]=aHit;
	//		return;
	//Split for better readability
	const float2 Ds=make_float2(D.x,D.y);
	const float2 Dh=make_float2(D.z,D.w);


	//Attenuation term at receiver due to propagation 
	z = make_float2(0.0f, -k*s);
	zexp = complex_exp_only_imaginary(z);
	float2 E_at_r = sca_complex_prod(-1.0f*A, zexp);

	//Some tests	
	//float2 Vbi=computeViB(zexp,Ds,Dh,s); 
	//float2 Vbr=computeVrB(zexp,Ds,Dh,s); 
	//float2 Vbir=Vbi+Vbr;
	//float2 Dhb=sca_complex_prod(-1.0*sqrtf(s),complex_prod(zexp,Vbir));
	//
	//	RDNHit aHit;
	//	aHit.EEx=make_float4(Vbi.x,Vbi.y, 0.0f,0.0f);
	//	aHit.EyEz=make_float4(0.0f,0.0f,0.0f,0.0f);
	//	difBuffer[index]=aHit;
	//return;
	//Complex amplitude of the diffracted E field at the receiver. Eq. [13-88] Balanis
	float2 Er_beta=complex_prod(E_at_r,complex_prod(Ds,Ei_beta_q)); //component of the diffracted E field parallel to the plane of diffraction at the receiver
	float2 Er_phi=complex_prod(E_at_r,complex_prod(Dh,Ei_phi_q)); //component of the diffracted E field perpendicular to the plane of diffraction at the  receiver
		//rtPrintf("%u \t Er_beta=(%6e,%6e) |Er_beta|=%6e  Er_phi=(%6e,%6e) |Er_phi|=%6e Ds=(%6e,%6e) Dh=(%6e,%6e))\n",launchIndex.y,Er_beta.x,Er_beta.y,length(Er_beta),Er_phi.x,Er_phi.y,length(Er_phi),D.x,D.y,D.z,D.w);


	//float2 Ee=Er_beta+Er_phi;

		//rtPrintf("%u \t Dhb=(%6e,%6e) Dh=(%6e,%6e)) angle(Vi+Vr)=%f angle(Er_phi)=%f a(E)=%f\n",launchIndex.y,Dhb.x,Dhb.y,Dh.x,Dh.y, (atan2(Vbir.y,Vbir.x)*180/M_PIf),(atan2(-Er_phi.y,-Er_phi.x)*180/M_PIf),(atan2(-Ee.y,-Ee.x)*180/M_PIf));

	//Get the unit vectors for the plane of diffraction. The above complex amplitude multiply the corresponding (beta and phi) unit vectors in the diffraction plane
	const float3 phi_u = n_dplane; 
	//s_unit vector is defined from DP to receiver, so we have to reverse ray here
	const float3 beta_o_u=normalize(cross(-rxRay,phi_u)); 
	//const float3 sss = normalize(cross(phi_u,beta_o_u));
		//rtPrintf("%u\t beta_o_u=(%f,%f,%f) phi_u=(%f,%f,%f) s=(%f,%f,%f) \n",launchIndex.x,beta_o_u.x,beta_o_u.y,beta_o_u.z,phi_u.x,phi_u.y,phi_u.z,sss.x,sss.y,sss.z);

	//Compute FIELD
		float3 xu=make_float3(1.0f,0.0f,0.0f);
		float3 yu=make_float3(0.0f,1.0f,0.0f);
		float3 zu=make_float3(0.0f,0.0f,1.0f);
		float2 Ex=sca_complex_prod(dot(beta_o_u,xu),Er_beta) + sca_complex_prod(dot(phi_u,xu),Er_phi);
		float2 Ey=sca_complex_prod(dot(beta_o_u,yu),Er_beta) + sca_complex_prod(dot(phi_u,yu),Er_phi);
		float2 Ez=sca_complex_prod(dot(beta_o_u,zu),Er_beta) + sca_complex_prod(dot(phi_u,zu),Er_phi);
			//printf("%u\t Ex=(%f,%f) |Ex|=%f Ey=(%f,%f) |Ey|=%f Ez=(%f,%f) |Ez|=%f \n",launchIndex.x,Ex.x, Ex.y,length(Ex),Ey.x,Ey.y,length(Ey),Ez.x,Ez.y,length(Ez));
		//float4 sangles=angles*57.2968f;
			//printf("%u\t%u\t%u  sangles(beta, beta',phi, phi')=(%f,%f,%f,%f) dif=%6e \n",launchIndex.x,launchIndex.y,launchIndex.z, sangles.x,sangles.y,sangles.z,sangles.w, (sangles.x-sangles.y));
		RDNHit aHit;
		float g=1;
		if (useAntennaGain) {
		
			g=getAntennaGain(rxRay, antennaGainIdBuffer[rxIndex],transformToPolarizationBuffer[rxIndex]);	
			Ex=sca_complex_prod(g,Ex);
			Ey=sca_complex_prod(g,Ey);
			Ez=sca_complex_prod(g,Ez);
			//printf("%u\t HDIF Ex=(%f,%f) |Ex|=%f Ey=(%f,%f) |Ey|=%f Ez=(%f,%f) |Ez|=%f g=%f \n",launchIndex.x,Ex.x, Ex.y,length(Ex),Ey.x,Ey.y,length(Ey),Ez.x,Ez.y,length(Ez),g);
		}

		//aHit.EEx=make_float4(0.0f,0.0f,Ex.x,Ex.y);
		//Additional output
		//float unfoldedPathLength = s+s_prime;
		//aHit.doaD = make_float4(rxRay.x, rxRay.y,rxRay.z, unfoldedPathLength);
		//aHit.doDu = make_float4(txRay.x, txRay.y,txRay.z, s);
	
		//difBuffer[index]=aHit;
	
	//Compute VOLTAGE

		//****************************
		//To get the induced voltage, we need to 
		//apply the dot product with the effective lenght of the received antenna. 

		//float3 ver_o; //Receiver vertical field vector
		//float3 hor_o; //Receiver horizontal field vector

		//Get polarization for receiver for this ray rxRay is already in the direction receiver to DP 
		//getLinearPolarizationInRayBasis(pol, rxRay,  hor_o,ver_o);

		//Get the  components of received field for the normal and parallel field vectors (geometric projection on receiver polarization vectors times reflection coefficients)
		//This would be equivalent to a dot product with the effective length (not the conjugated beacuse we already reversed and it is a linear polarization anyway)
		//const float2 Einorm=sca_complex_prod(dot(beta_o_u,hor_o),Er_beta) + sca_complex_prod(dot(phi_u,hor_o),Er_phi);
		//const float2 Eipar=sca_complex_prod(dot(beta_o_u,ver_o),Er_beta) + sca_complex_prod(dot(phi_u,ver_o),Er_phi);
		//const float2 Einorm=sca_complex_prod(dot(beta_o_u,beta_p),Er_beta) + sca_complex_prod(dot(phi_u,beta_p),Er_phi);
		//const float2 Eipar=sca_complex_prod(dot(beta_o_u,phi_p),Er_beta) + sca_complex_prod(dot(phi_u,phi_p),Er_phi);
		//float2 E=Einorm+Eipar;

		//The above formulation is equivalent to this below 
		//Geometric part due to polarization at the receiver
		//Get polarization for receiver for this ray rxRay is already in the direction receiver to DP 
		//const float3 Er_pol=getLinearPolarizationForRaySimple(pol,rxRay); //It is normalized
		const float3 Er_pol=getLinearPolarizationForRaySimple(receiverPolarizationBuffer[rxIndex],rxRay); //It is normalized

			//rtPrintf("%u\t%u\t%u\trxRay=(%f,%f,%f) Er_pol=(%f,%f,%f) Ev=(%f,%f,%f) A=%f )\n",launchIndex.x, launchIndex.y,launchIndex.z,rxRay.x,rxRay.y,rxRay.z,Er_pol.x,Er_pol.y,Er_pol.z,Ev.x,Ev.y,Ev.z,A);

		const float Er_beta_v=dot(beta_o_u,Er_pol); 
		const float Er_phi_v=dot(phi_u,Er_pol); 
		Er_beta=sca_complex_prod(Er_beta_v,Er_beta);
		Er_phi=sca_complex_prod(Er_phi_v,Er_phi);
		//This is actually the induced voltage on the antenna. From it we can compute the received power
		float2 E=Er_beta+Er_phi;
		if (useAntennaGain) {

			E=sca_complex_prod(g,E);
		}
		//float4 sangles=angles*57.2968f;

			//rtPrintf("%u\t%u\t%u  E=(%f,%f) E_b=(%f,%f) E_phi=(%f,%f) L=%f dif=%6e  \n",launchIndex.x,launchIndex.y,launchIndex.z, E.x,E.y,Er_beta.x,Er_beta.y,Er_phi.x, Er_phi.y, L, (sangles.x-sangles.y));
			//rtPrintf("G\t |E|=%6e index=(%u,%u,%u) %f \n",length(E),  index.x,index.y,index.z,angles.z*57.2968f);

			//rtPrintf("%u\t%u\t%u  sangles(beta, beta',phi, phi')=(%f,%f,%f,%f) L=%f dif=%6e  \n",launchIndex.x,launchIndex.y,launchIndex.z, sangles.x,sangles.y,sangles.z,sangles.w, L, (sangles.x-sangles.y));
		//RDNHit aHit;
		//aHit.EEx=make_float4(E.x,E.y, 0.0f,0.0f);
		aHit.EEx=make_float4(E.x,E.y,Ex.x,Ex.y); //Use 1 on EEx.x as flag for real hit
		aHit.EyEz=make_float4(Ey.x,Ey.y,Ez.x,Ez.y);
		//Additional output
		float unfoldedPathLength = s+s_prime;
		aHit.doaD = make_float4(rxRay.x, rxRay.y,rxRay.z, unfoldedPathLength);
		aHit.doDu = make_float4(txRay.x, txRay.y,txRay.z, s);
	
		difBuffer[launchIndex]=aHit;
		

}


//Diffraction Launch program
RT_PROGRAM void computeSingleDiffraction() {

	//1D launch [number_or_real_hist]
		//rtPrintf("%u\t%u\t%u Launch \n",launchIndex.x,launchIndex.y,launchIndex.z);
	const uint difBufferIndex=launchIndex;
	//Initialize buffer to make sure it does not carry values from previous launches
	RDNHit initHit;
	initHit.EEx=make_float4(0.0f,0.0f, 0.0f,0.0f);
	
	initHit.EyEz=make_float4(0.0f,0.0f,0.0f,0.0f);
	difBuffer[difBufferIndex]=initHit;
	const uint3 rhi=losBuffer[difBufferIndex];

	const uint edgeIndex=	rhi.x;
	const uint rxIndex = rhi.y;
	const uint txIndex = rhi.z;

	Transmitter tx = txBuffer[txIndex];

	const float3 origin = make_float3(tx.origin_p);
	const float4 sphere = receiverPositionsBuffer[rxIndex];
	//Check if ray is hitting his own tx (transmitter are also receivers usually) A transmitter cannot receive while it is transmitting, unless other channel is used.
	if (static_cast<int>(sphere.w)==tx.externalId) {
		return;	
	}	
	const float3 destination =make_float3(sphere.x,sphere.y,sphere.z);
	Edge e = edgeBuffer[edgeIndex];
	if (!isDiffractingEdge(origin,e)) {
		//rtPrintf("%u\t%u\t%u Not diffracting edge %d\n",launchIndex.x,launchIndex.y,launchIndex.z,e.id);
		return;
	}

	//Compute diffraction point (DP) between transmitter, receiver and edge
	float3 dp; 
	if (computeDiffractionPoint(origin,destination,e,dp)){
		//rtPrintf("%u\t%u\t%u e=%u dp=(%f,%f,%f) \n",launchIndex.x,launchIndex.y,launchIndex.z,e.id,dp.x,dp.y,dp.z);
		float3 originToDP=dp-origin;
		float dist_originToDp=length(originToDP);
		float3 txRayDirection = originToDP/dist_originToDp;
		float3 destinationToDP=dp-destination;
		//rtPrintf("%u\t%u\t%u e=%u destinationToDP=(%f,%f,%f) tx not blocked \n",launchIndex.x,launchIndex.y,launchIndex.z,e.id, destinationToDP.x,destinationToDP.y,destinationToDP.z);
		float dist_destinationToDp=length(destinationToDP);
		float3 rxRayDir=destinationToDP/dist_destinationToDp;
		float3 rxRay=rxRayDir;
		//We are not going to trac, but keep it for backward compatibility
		VisibilityPayload visibilityRayPayload;
		visibilityRayPayload.polarization_k = tx.polarization_k; 
		visibilityRayPayload.result.x=OPAL_DIFFRACTION_LOS;
		visibilityRayPayload.faces=e.faces; //Faces of the edge, where the ray can hit to compute the diffraction. If it hits any other face, there is no LoS
		visibilityRayPayload.result.y=0;
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
		//rtPrintf("%u\t%u\t%u e=%u dp=(%f,%f,%f) rx not blocked \n",launchIndex.x,launchIndex.y,launchIndex.z,e.id, dp.x,dp.y,dp.z);
		float gain;
		if (useAntennaGain) {
			//Tx gain here. Rx gain in included in computeElectricFieldAtReceiver
			//rtPrintf("tx useAntennaGain\n");
			rtBufferId<float,2> bid=tx.gainId;
			const Matrix<4,4> tp=tx.transformToPolarization;
			gain=getAntennaGain(txRayDirection, bid, tp);
			//printf("%u\t%u\t%u gain=%f\n",launchIndex.x,launchIndex.y,launchIndex.z,gain) ; 
		} else {
			gain=1.0f;
		}
		computeElectricFieldAtReceiver(angles, spar, txRayDirection, tx.polarization_k,n_iplane,n_dplane,e.pn.w, rxRayDir, R_0, R_n, gain);
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
		printf("Diffraction computation user exception %d at (%d, %d)\n", code - RT_EXCEPTION_USER, launchIndex);
	}
	else
	{
		printf("Diffraction computation Exception code 0x%X at (%d, %d)\n", code, launchIndex);
	}

}

