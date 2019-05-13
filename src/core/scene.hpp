#pragma once
#include "surface.hpp"
#include "light.hpp"
#include <vector>
#include "surface_instance.hpp"

namespace HydrogenCG
{

	
	/*!
		\brief		A structure that holds a collection of surfaces and lights

		Provides intersection methods with the collection.
	*/
	struct Scene
	{
		SurfaceInstanceCollection collection;					//!< A collection of surface instances
		std::vector<std::unique_ptr<Light>> lights;				//!< A collection of all lights in the scene

		// convenience/brevity
		inline InstanceIntersection operator()(const Ray& ray, float minT = 0.0f, float maxT = std::numeric_limits<float>::infinity()) const
		{
			return collection.intersect(ray, minT, maxT);
		}

		bool intersectAny(const Ray& ray, float minT = 0.0f, float maxT = std::numeric_limits<float>::infinity()) const
		{
			return collection.intersectAny(ray, minT, maxT);
		}

		float distance(const vec3& p) const
		{
			return collection.distance(p);
		}
	};
}