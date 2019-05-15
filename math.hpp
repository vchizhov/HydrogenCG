#pragma once
/*
	@author: Vassillen Chizhov, 2019
	Some convenience math functions
*/

#include <math.h>
#undef INFINITY

namespace HydrogenCG
{
	const float INFINITY = std::numeric_limits<float>::infinity();

	float min(float lhs, float rhs)
	{
		return lhs <= rhs ? lhs : rhs;
	}

	float max(float lhs, float rhs)
	{
		return lhs >= rhs ? lhs : rhs;
	}

	float clamp(float arg, float minVal, float maxVal)
	{
		return max(min(arg, maxVal), minVal);
	}
}
