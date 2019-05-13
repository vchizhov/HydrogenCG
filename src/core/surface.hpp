#pragma once
#include "vec.hpp"
#include "mat.hpp"
#include "ray.hpp"
#include <limits>

namespace HydrogenCG
{
	struct Surface;
	/*!
		\brief		Return type for the scene intersection routine

		Provides the distance of the intersection from the ray origin
		(provided that the ray direction is normalized), as well as a
		pointer to the intersected surface. A cast to bool is defined
		for convenience which returns true if an intersection occurred.
	*/
	struct Intersection
	{
		float t;	//!< Distance to the intersection from the ray origin, defaults to infinity
		vec3 n;		//!< Normal at the intersection
		vec2 uv;	//!< Surface coords at the intersection
		Intersection() : t(std::numeric_limits<float>::infinity()) {}
		Intersection(float t, const vec3& n, const vec2& uv) : t(t), n(n), uv(uv) {}
		/*!
			\brief		Convenience bool cast

			\return		Returns true if an intersection occurred, false otherwise.
		*/
		inline explicit operator bool() const
		{
			// if the intersection is at any point closer than infinity, it's valid
			return t < std::numeric_limits<float>::infinity();
		}

		bool operator<=(const Intersection& rhs)
		{
			return t <= rhs.t;
		}
	};

	/*!
		\brief	A base class for all surfaces that we will provide intersect and normal functions for.
	*/
	struct Surface
	{

		virtual ~Surface() {}


		explicit Surface() {}
		/*!
			\brief		Returns the closest surface intersection within the (minT,maxT) range

			\return		Returns maxT+1 in case no intersection occurs, otherwise returns
						the distance from the ray origin to the closest intersection point
						with the surface in the direction of the ray.
		*/
		virtual Intersection intersect(const Ray& ray, float minT = 0.0f, float maxT = std::numeric_limits<float>::infinity()) const = 0;
		/*!
			\brief		Returns the normal (orthogonal unit vector) of a surface at a point p on it

			\return		The result can be different if the point is not on
						(close enough) to the surface.
		*/

		virtual float distance(const vec3& p) const = 0;
	};


}