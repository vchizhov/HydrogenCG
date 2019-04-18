#pragma once
#include "vec.h"
#include "ray.h"

namespace HydrogenCG
{
	struct Camera
	{
		vec3 e0, e1, e2;
		vec3 pos;

		Camera(const vec3& pos = vec3(0), const vec3& right = vec3(1.0f,0.0f,0.0f), 
			const vec3& up = vec3(0.0f, 1.0f,0.0f), const vec3& forward = vec3(0.0f, 0.0f, 1.0f)) 
			: e0(right), e1(up), e2(forward) {}

		Ray operator()(float u, float v) const
		{
			return Ray(pos, normalize(u * e0 + v * e1 + e2));
		}
	};
}