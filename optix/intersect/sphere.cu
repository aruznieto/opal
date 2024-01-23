/***************************************************************/
//
//Copyright (c) 2019 Esteban Egea-Lopez http://girtel.upct.es/~simulations
//
/**************************************************************/


//Modified from NVIDIA OptiX samples

/*
* Copyright (c) 2016, NVIDIA CORPORATION. All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions
* are met:
*  * Redistributions of source code must retain the above copyright
*    notice, this list of conditions and the following disclaimer.
*  * Redistributions in binary form must reproduce the above copyright
*    notice, this list of conditions and the following disclaimer in the
*    documentation and/or other materials provided with the distribution.
*  * Neither the name of NVIDIA CORPORATION nor the names of its
*    contributors may be used to endorse or promote products derived
*    from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
* PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
* CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
* EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
* PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
* OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "../../Common.h"
#include <optixu/optixu_math_namespace.h>
#include <optixu/optixu_aabb_namespace.h>
using namespace optix;





//Sphere variables
rtDeclareVariable(float4, sphere, , );
rtDeclareVariable(SphereHit, hit_attr, attribute hit_attr, );
rtDeclareVariable(optix::Ray, ray, rtCurrentRay, );

template<bool use_robust_method>
static __device__
void intersect_sphere(void)
{
	float3 center = make_float3(sphere);
	float3 O = ray.origin - center;
	float3 D = ray.direction;
	float radius = sphere.w;

	float b = dot(O, D);
	float c = dot(O, O) - radius*radius;
	float disc = b*b - c;
	if (disc > 0.0f) {
		float sdisc = sqrtf(disc);
		float root1 = (-b - sdisc);

		bool do_refine = false;

		float root11 = 0.0f;

		if (use_robust_method && fabsf(root1) > 10.f * radius) {
			do_refine = true;
		}

		if (do_refine) {
			// refine root1
			float3 O1 = O + root1 * ray.direction;
			b = dot(O1, D);
			c = dot(O1, O1) - radius*radius;
			disc = b*b - c;

			if (disc > 0.0f) {
				sdisc = sqrtf(disc);
				root11 = (-b - sdisc);
			}
		}

		bool check_second = true;
		float t = root1 + root11;
		if (rtPotentialIntersection(t)) {
			//shading_normal = geometric_normal = (O + (root1 + root11)*D) / radius;
			SphereHit h;
//			h.t = t;
//			h.geom_normal = (O + (root1 + root11)*D) / radius;
		float3 gn = (O + (root1 + root11)*D) / radius;		
	//packed version
			h.geom_normal_t=make_float4(gn.x,gn.y,gn.z,t);

			hit_attr = h;
			//Only one material. Change here if more materials used
			if (rtReportIntersection(0)) {
				check_second = false;
			}
		}
		if (check_second) {
			float root2 = (-b + sdisc) + (do_refine ? root1 : 0);
			if (rtPotentialIntersection(root2)) {
				SphereHit h;
				//h.t = root2;
				//h.geom_normal = (O + root2*D) / radius;
				float3 gn2 = (O + root2*D) / radius;
				//Packed version
				h.geom_normal_t=make_float4(gn2.x,gn2.y,gn2.z,root2);
				hit_attr = h;
				//shading_normal = geometric_normal = (O + root2*D) / radius;
			
				//Only one material. Change here if more materials used
				rtReportIntersection(0);
			}
		}
	}
}


RT_PROGRAM void intersectSphere(int primIdx)
{
	intersect_sphere<false>();
}


RT_PROGRAM void robust_intersectSphere(int primIdx)
{
	intersect_sphere<true>();
}
RT_PROGRAM void rtgem_intersectSphere(int primIdx) {
	float3 center = make_float3(sphere);
	float3 f = ray.origin - center;
	float3 d = ray.direction; 
	float radius = sphere.w;
	//float a=dot(d,d); Removed because we assume the ray direction is already normalized

	float b=-1.0f*dot(f,d);
	float3 aux=f+(b*d);
	float rsq=radius*radius;
	float disc=rsq - dot(aux,aux);
	if (disc > 0.0f) { //A tangential hit is not considered a hit
		float c=dot(f,f)-rsq;
		float s=((b > 0) ? 1.0f : -1.0f); 
		float q=b+(s*sqrt(disc));
		float t=c/q;
		bool check_second = true;
		if (rtPotentialIntersection(t)) {
			//rtPrintf("primIdx=%d,s=%f,q=%f,t=%f\n",primIdx,s,q,t);
			//shading_normal = geometric_normal = (O + (root1 + root11)*D) / radius;
			SphereHit h;
//			h.t = t;
//			h.geom_normal = (O + (root1 + root11)*D) / radius;
			float3 gn = (f + (t*d)) / radius;		
			//packed version
			h.geom_normal_t=make_float4(gn.x,gn.y,gn.z,t);

			hit_attr = h;
			//Only one material. Change here if more materials used
			if (rtReportIntersection(0)) {
				check_second = false;
			}
		}
		if (check_second) {
			t=q;	
			if (rtPotentialIntersection(t)) {
				SphereHit h;
				//rtPrintf("primIdx=%d,2s=%f,q=%f,t=%f\n",primIdx,s,q,t);
				float3 gn2 = (f + (t*d)) / radius;		
				//Packed version
				h.geom_normal_t=make_float4(gn2.x,gn2.y,gn2.z,t);
				hit_attr = h;
				//Only one material. Change here if more materials used
				rtReportIntersection(0);
			}
		}
	}

}


RT_PROGRAM void boundsSphere(int, float result[6])
{
	const float3 cen = make_float3(sphere);
	const float3 rad = make_float3(sphere.w);

	optix::Aabb* aabb = (optix::Aabb*)result;

	if (rad.x > 0.0f && !isinf(rad.x)) {
		aabb->m_min = cen - rad;
		aabb->m_max = cen + rad;
	}
	else {
		aabb->invalidate();
	}
}

