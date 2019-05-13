#pragma once

namespace HydrogenCG
{
	template<typename T>
	T min(T lhs, T rhs)
	{
		return lhs <= rhs ? lhs : rhs;
	}

	template<typename T>
	T max(T lhs, T rhs)
	{
		return lhs >= rhs ? lhs : rhs;
	}

	template<typename T>
	T clamp(T arg, T minVal, T maxVal)
	{
		return max(min(arg, maxVal), minVal);
	}

	template<typename T>
	inline T sign(T scalar)
	{
		return static_cast<T>(T(0) < scalar) - static_cast<T>(T(0) > scalar);
	}

	// https://en.wikipedia.org/wiki/Smoothstep
	// we do not clamp
	float smoothstepCubic(float x)
	{
		return x*x*(3.0f - 2.0f*x);
	}

	// Ken Perlin, "Improving Noise"
	// we do not clamp
	float smoothstepQuintic(float x)
	{
		return x*x*x*(x*(6.0f*x-15.0f)+10.0f);
	}

	const float EPSILON = 0.0001f;
	const float pi = static_cast<float>(3.141592653589793238L);
}