#pragma once
#include "vec.h"
#include "ray.h"
#include <limits>

namespace HydrogenCG
{
	/*!
		\brief	A base class for all surfaces that we will provide intersect and normal functions for.
	*/
	struct Surface
	{
		vec3 col;	//!< albedo of the surface
		Surface(const vec3& col = vec3(1)) : col(col) {}
		/*!
			\brief		Returns the closest surface intersection within the (minT,maxT) range

			\return		Returns maxT+1 in case no intersection occurs, otherwise returns
						the distance from the ray origin to the closest intersection point
						with the surface in the direction of the ray.
		*/
		virtual float intersect(const Ray& ray, float minT = 0.0f, float maxT = std::numeric_limits<float>::max()) const = 0;
		/*!
			\brief		Returns the normal (orthogonal unit vector) of a surface at a point p on it

			\return		The result can be different if the point is not on
						(close enough) to the surface.
		*/
		virtual vec3 normal(const vec3& p) const = 0;
	};

	/*!
		\brief	A mathematical description of a sphere.
	*/
	struct Sphere : public Surface
	{
		vec3 o;		//!< the origin/center of the sphere
		float r;	//!< The radius of the sphere

		Sphere() {}
		Sphere(const vec3& origin, float radius, const vec3& col = vec3(1)) : o(origin), r(radius), Surface(col) {}

		/*!
			\brief		Sphere intersection
		*/
		float intersect(const Ray& ray, float minT = 0.0f, float maxT = std::numeric_limits<float>::max()) const final
		{
			/*
				Sphere intersection derivation, 
				||x|| is used to denote the euclidean norm of vector x 
				<a,b> is the dot(inner) product of vectors a and b
				r^2 = r*r

				||ray(t) - pos|| == r <->
				||ray(t) - pos||^2 == r^2 <->
				<ray(t)-o,ray(t)-o> == r^2 <->
				<ray.o-o + t * ray.d, ray.o-o + t * ray.d> == r^2 <->
				<ray.d,ray.d> * t^2 - 2 * <ray.d, o - ray.o> * t + <o-ray.o,o-ray.o> - r^2 == 0
				A = <ray.d,ray.d>, but |ray.d|==1 if the direction is normalized -> A = 1
				B = <ray.d, o - ray.o>
				C = <o-ray.o,o-ray.o> - r^2

				D = B^2 - C
				D<0  -> no intersection
				D==0 -> grazing intersection
				D>0  -> 2 intersections

				One can ignore grazing intersections since they are actually numerical error.

				D>0 -> sqrtD = sqrt(D)
				t_1 = B - D
				t_2 = B + D

				If maxT>t_1>minT -> intersection at ray(t_1)
				Else If maxT>t_2>min_t -> intersection at ray(t_2)
				Else -> no intersection
			*/

			// vector from the ray origin to the sphere origin
			vec3 ro = o - ray.o;
			// quadratic polynomial terms (a==1)
			float b = dot(ray.d, ro);
			float c = dot(ro, ro) - r*r;
			// discriminant
			float discriminant = b * b - c;

			// If the discriminant is 0 or negative -> no (real) intersection
			if (discriminant <= 0.0f) return maxT+1.0f;

			float sqrtD = sqrt(discriminant);
			float t1 = b - sqrtD; // closer intersection (where the ray enters the sphere)
			// if it is within the defined ray segment by (minT,maxT) return it as the closest intersection
			if (t1 > minT && t1 < maxT) return t1;

			float t2 = b + sqrtD; // further intersection (where the ray exits the sphere)
			// if it is within the defined ray segment by (minT,maxT) return it as the closest intersection
			if (t2 > minT && t2 < maxT) return t2;

			// otherwise return no intersection
			return maxT + 1.0f;
		}

		/*!
			\brief		Computes the sphere normal for a point p

			\param[in]	p A point in space

			\return		Returns the outwards facing unit normal at point p 
						if p is on the surface of the sphere, otherwise 
						returns a scaled version of the normal
		*/
		vec3 normal(const vec3& p) const final
		{
			return (p - o) / r;
		}
	};
}