#pragma once
#include "..\core\surface_instance.hpp"

namespace HydrogenCG
{
	/*!
		\brief	A parametrical plane equation.
	*/
	struct Plane : public Surface
	{
		vec3 o;		//!< some point on the plane
		// note: the two basis vectors should not be collinear (parallel)
		vec3 e1;	//!< The first basis vector spanning the plane
		vec3 e2;	//!< The second basis vector spanning the plane

		Plane() {}
		explicit Plane(const vec3& origin, const vec3& e1, const vec3& e2, const vec3& col = vec3(1)) : o(origin), e1(e1), e2(e2) {}

		/*!
			\brief		Plane intersection
		*/
		Intersection intersect(const Ray& ray, float minT = 0.0f, float maxT = std::numeric_limits<float>::infinity()) const final
		{
			/*
				Derivation:
				Any point on a parametrically defined plane can be written as:
				p(a,b) = o + u * e1 + v * e2; a,b in (-inf,inf)
				Any point on a ray can be written as:
				p(t) = r.o + t * r.d; t in [0,inf)

				Equating these two yields:
				o + u * e1 + v * e2 = r.o + t * r.d
				t * (-r.d) + u * e1 + v * e2 = r.o - o
				which is a linear system of equations that can be written in matrix form as:
				M * (t,u,v) = r.o - o
				where M has columns: -r.d, e1, e2
				the formal solution of the system is:
				(t,u,v) = inverse(M) * (r.o-o)

				Finally we check that t falls within the predefined bounds (minT,maxT)
			*/
			vec3 tuv = inverse(fromColumns(-ray.d, e1, e2)) * (ray.o - o);
			return minT < tuv.x && tuv.x < maxT ? Intersection(tuv.x, normalize(cross(e1,e2)), tuv.yz) : Intersection();
		}

		/*!
			\brief Returns the (closest) distance from a point p to the plane
		*/
		float distance(const vec3& p) const final
		{
			/*
				Derivation:
				Let q be orthogonal projection of p onto the plane, then 
				q = p - k*n, and <q-o,n> = 0, thus:
				<p-k*n-o,n> = 0
				k*<n,n> = <p-o,n>
				k = <p-o,n>
			*/
			return dot(p-o,normalize(cross(e1,e2)));
		}

		static const std::shared_ptr<const Plane> yPlane;

	};
	const std::shared_ptr<const Plane> Plane::yPlane = std::make_shared<Plane>(vec3(0), vec3(1,0,0), vec3(0,0,1));
	std::shared_ptr<SurfaceInstanceBSDF> plane(const vec3& origin, const vec3& e1, const vec3& e2, std::shared_ptr<BSDF> bsdf = BSDF::white)
	{
		return std::make_shared<SurfaceInstanceBSDF>(std::make_shared<Plane>(origin, e1, e2),bsdf);
	}
}