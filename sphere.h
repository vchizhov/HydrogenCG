#pragma once
#include "vec.h"
#include "ray.h"
#include <limits>

namespace HydrogenCG
{
	struct Surface
	{
		vec3 col;
		Surface(const vec3& col = vec3(1)) : col(col) {}
		virtual float intersect(const Ray& ray, float minT = 0.0f, float maxT = std::numeric_limits<float>::max()) const = 0;
		virtual vec3 normal(const vec3& pointInSpace) const = 0;
	};

	struct Sphere : public Surface
	{
		vec3 pos;
		float radius;

		Sphere() {}
		Sphere(const vec3& pos, float radius, const vec3& col = vec3(1)) : pos(pos), radius(radius), Surface(col) {}

		// Evaluates whether the ray intersects the sphere withint the ray segment (minT, maxT)
		// returns minT - 1.0f if there is no intersection, so that a check for intersection can be performed as intersect(...)<minT ?
		// otherwise returns the distance from the ray origin to the closest intersection
		float intersect(const Ray& ray, float minT = 0.0f, float maxT = std::numeric_limits<float>::max()) const final
		{
			/*
				|ray(t) - pos| == r <->
				|ray(t) - pos|^2 == r^2 <->
				<ray(t)-pos,ray(t)-pos> == r^2 <->
				<ray.o-pos + t * ray.d, ray.o-pos + t * ray.d> == r^2 <->
				<ray.d,ray.d> * t^2 - 2 * <ray.d, pos - ray.o> * t + <pos-ray.o,pos-ray.o> - r^2 == 0
				A = <ray.d,ray.d>, but |ray.d|==1 if the direction is normalized -> A = 1
				B = <ray.d, pos - ray.o>
				C = <pos-ray.o,pos-ray.o> - r^2

				D = B^2 - C
				D<0  -> no intersection
				D==0 -> grazing intersection
				D>0  -> 2 intersections

				One can ignore grazing intersections since they are actually numerical error.

				D>0 -> sqrtD = sqrt(D)
				t_1 = B - D
				t_2 = B + D

				If maxT>t_1>minT -> intersection at ray(t_1)
				Else If maxT>t_2>min_t -> intersection at ray(t_2)
				Else -> no intersection
			*/
			vec3 oPos = pos - ray.o;
			float b = dot(ray.d, oPos);
			float c = dot(oPos, oPos) - radius*radius;
			float d = b * b - c;

			// If the discriminant is 0 or negative -> no (actual) intersection
			if (d <= 0.0f) return minT - 1.0f;

			float sqrtD = sqrt(d);
			float t1 = b - sqrtD; // closer intersection
			// is it within the defined ray segment by (minT,maxT) ?
			if (t1 > minT && t1 < maxT) return t1;

			float t2 = b + sqrtD; // further intersection
			// is it withint the defined ray segment by (minT,maxT) ?
			if (t2 > minT && t2 < maxT) return t2;

			return minT - 1.0f;
		}

		// Returns the normal for any point on the sphere, and scaled normals for points not on the sphere
		vec3 normal(const vec3& pointInSpace) const final
		{
			return (pointInSpace - pos) / radius;
		}
	};
}