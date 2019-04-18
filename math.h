#pragma once
/*
	@author: Vassillen Chizhov, 2019
	Raycaster, part 0: Ray tracing a sphere
*/

namespace HydrogenCG
{
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
