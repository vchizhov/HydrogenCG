#pragma once
#include "sphere.h"
#include "light.h"
#include <vector>

namespace HydrogenCG
{
	/*
		\brief		Return type for the scene intersection routine

		Provides the distance of the intersection from the ray origin
		(provided that the ray direction is normalized), as well as a 
		pointer to the intersected surface. A cast to bool is defined 
		for convenience which returns true if an intersection occurred.
	*/
	struct Intersection
	{
		float t;		//!< Distance to the intersection from the ray origin, defaults to -1
		Surface* s;		//!< Pointer to the intersected surface, defaults to nullptr
		Intersection() : s(nullptr), t(-1.0f) {}

		/*!
			\brief		Convenience bool cast

			\return		Returns true if an intersection occurred, false otherwise.
		*/
		inline operator bool() const
		{
			return s!=nullptr;
		}

	};

	/*
		\brief		A structure that holds a collection of surfaces and lights

		Provides intersection methods with the collection.
	*/
	struct Scene
	{
		std::vector<Surface*> surfaces;		//!< The collection of all surfaces in the scene
		std::vector<LightPoint> lights;		//!< The collection of all pointlights in the scene

		~Scene() 
		{
			for (int i = 0; i < surfaces.size(); ++i) delete surfaces[i]; 
		}

		// convenience
		Surface* operator[](uint32_t i) { return surfaces[i]; }
		const Surface* operator[](uint32_t i) const { return surfaces[i]; }

		/*!
			\brief			An intersection routine with the collection

			\param[in]		ray	The ray to be tested for an intersection with
			\param[in]		minT The lower bound above which an intersection is accepted
			\param[in]		maxT The upper bound below which an intersection is accepted

			\return			A structure filled with the distance to the intersection and
							a pointer to the intersected surface, the distance returns in 
							the case that no intersection occurred is maxT
		*/
		Intersection intersect(const Ray& ray, float minT = 0.0f, float maxT = std::numeric_limits<float>::max()) const
		{
			Intersection intersection;

			// initialize distance to maximum draw distance
			// this will be used as the closest intersection so far
			intersection.t = maxT;

			// iterate over all surfaces and find the intersection with the closest one
			// what this algorithm essentially does is finding the minimum of an "array"
			// in this case the "array" is implicit and is defined by the sequence of values 
			// produced by calling intersect on each surface
			for (auto s : surfaces)
			{
				// we feed in the currently closest intersection as the upper bound,
				// thus if the intersection in the current iteration is not closer 
				// than intersection.t (the closest thus far), we would get no intersection
				// Since this is the case, the invariant that intersection.t contains the closest
				// intersection so far will be maintained
				float t = s->intersect(ray, minT, intersection.t);

				// intersect returns maxT+1 if no intersection occurrs
				// so we use the check t<maxT as an indicator for an intersection
				if (t < intersection.t)
				{
					intersection.s = s;
					intersection.t = t;
				}
			}

			return intersection;
		}

		// convenience/brevity
		inline Intersection operator()(const Ray& ray, float minT = 0.0f, float maxT = std::numeric_limits<float>::max()) const
		{
			return intersect(ray, minT, maxT);
		}

		/*!
			\brief			Returns true at the first found intersection, false otherwise

			Useful for testing for visibility between two points, 
			set maxT to the distance between those minus epsilon in that case.
		*/
		bool intersectAny(const Ray& ray, float minT = 0.0f, float maxT = std::numeric_limits<float>::max()) const
		{
			for (auto s : surfaces)
			{
				if (s->intersect(ray, minT, maxT)<maxT) return true;
			}
			return false;
		}
	};
}