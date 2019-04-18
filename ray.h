#pragma once
#include "vec.h"


namespace HydrogenCG
{
	struct Ray
	{
		vec3 o;
		vec3 d;

		Ray() {}
		Ray(const vec3& o, const vec3& d) : o(o), d(d) {}

		// I do not like this syntax, as it looks funky when used.
		// I propose you rename this function to "at"
		vec3 operator()(float t) const { return o + t * d; }
	};
}
