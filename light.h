#pragma once
#include "vec.h"

namespace HydrogenCG
{
	/*!
		\brief A basic point light structure.
	*/
	struct LightPoint
	{
		vec3 o;		//!< Position in 3d space of the point light
		vec3 i;		//!< Intensity = Flux / 4PI, the "color" and strength of the light

		LightPoint(const vec3& origin = vec3(0), const vec3& intensity = vec3(1)) : o(origin), i(intensity) {}
	};
}