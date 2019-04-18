#pragma once
#include "vec.h"

namespace HydrogenCG
{

	struct LightPoint
	{
		vec3 pos;
		vec3 intensity;

		LightPoint(const vec3& pos = vec3(0), const vec3& intensity = vec3(1)) : pos(pos), intensity(intensity) {}
	};
}