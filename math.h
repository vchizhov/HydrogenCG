#pragma once

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