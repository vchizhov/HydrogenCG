#pragma once
#include "..\core\surface.hpp"
#include "..\core\mat.hpp"

namespace HydrogenCG
{
	/*!
		\brief	A triangle defined by a vertex and 2 edges
	*/
	struct Triangle : public Surface
	{
		vec3 o;		//! "origin" vertex = v0
		vec3 e1;	//! v1-v0
		vec3 e2;	//! v2-v0

		Triangle() {}
		explicit Triangle(const vec3& v0, const vec3& v1, const vec3& v2) : o(v0), e1(v1-v0),e2(v2-v0) {}

		/*!
			\brief		Triangle intersection
		*/
		Intersection intersect(const Ray& ray, float minT = 0.0f, float maxT = std::numeric_limits<float>::infinity()) const final
		{
			// see plane for an explanation
			vec3 tuv = inverse(fromColumns(-ray.d, e1, e2)) * (ray.o-o);
			float t = tuv.x;
			float u = tuv.y;
			float v = tuv.z;

			// u,v,(1-u-v) are barycentric coordinates, a triangle is defined by considering only the region 0<u,0<v, u+v<1
			return minT < t && t < maxT && 0 < u && 0 < v && u + v < 1 ? Intersection(t, normalize(cross(e1,e2)), tuv.yz) : Intersection();
			
		}

		/*!
			\brief		Computes the triangle normal

			\param[in]	p A point in space - redundant, is there in order to be consistent with the interface

			\return		The normal of the triangle
		*/
		vec3 normal(const vec3& /*p*/) const
		{
			return normalize(cross(e1, e2));
		}
	};
}