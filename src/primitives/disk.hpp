#pragma once
#include "..\core\surface.hpp"

namespace HydrogenCG
{
	/*!
		\brief	A disk defined through its origin and 2 edges (you can define a squished elliptical disk using this)
	*/
	struct Disk : public Surface
	{
		vec3 o;		//! "origin" vertex = v0
		vec3 e1;	//! v1-v0
		vec3 e2;	//! v2-v0

		Disk() {}
		explicit Disk(const vec3& origin, const vec3& e1, const vec3& e2) : o(origin), e1(e1), e2(e2) {}

		/*!
			\brief		Plane intersection
		*/
		Intersection intersect(const Ray& ray, float minT = 0.0f, float maxT = std::numeric_limits<float>::infinity()) const final
		{

			// see plane for an explanation
			vec3 tuv = inverse(fromColumns(-ray.d, e1, e2)) * (ray.o - o);
			float t = tuv.x;
			float u = tuv.y;
			float v = tuv.z;

			// a disk is defined by filtering out all points||(u,v)||>=1
			return minT < t && t < maxT && u*u+v*v<1 ? Intersection(t, normalize(cross(e1, e2)), tuv.yz) : Intersection();
		}

		/*!
			\brief		Computes the disk's normal

			\param[in]	p A point in space - redundant, is there in order to be consistent with the interface

			\return		The normal of the plane
		*/
		vec3 normal(const vec3& /*p*/) const
		{
			return normalize(cross(e1, e2));
		}
	};
}