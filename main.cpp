/*
	@author: Vassillen Chizhov, 2019
	Raycaster, part 0: Ray tracing a sphere
*/

#include <limits>
#include "vec.h"
#include "image.h"

using namespace HydrogenCG;
struct Ray
{
	vec3 o;
	vec3 d;

	Ray() {}
	Ray(const vec3& o, const vec3& d) : o(o), d(d) {}

	vec3 operator()(float t) const { return o + t * d; }
};


struct Camera
{
	vec3 e0, e1, e2;
	vec3 pos;

	Camera(const vec3& pos = vec3(0), const vec3& right = vec3(1.0f, 0.0f, 0.0f),
		const vec3& up = vec3(0.0f, 1.0f, 0.0f), const vec3& forward = vec3(0.0f, 0.0f, 1.0f))
		: e0(right), e1(up), e2(forward) {}

	Ray operator()(float u, float v) const
	{
		return Ray(pos, normalize(u * e0 + v * e1 + e2));
	}
};

struct Sphere
{
	vec3 pos;
	float radius;
	vec3 col;

	Sphere() {}
	Sphere(const vec3& pos, float radius, const vec3& col = vec3(1)) : pos(pos), radius(radius), col(col) {}

	// Evaluates whether the ray intersects the sphere withint the ray segment (minT, maxT)
	// returns minT - 1.0f if there is no intersection, so that a check for intersection can be performed as intersect(...)<minT ?
	// otherwise returns the distance from the ray origin to the closest intersection
	float intersect(const Ray& ray, float minT = 0.0f, float maxT = std::numeric_limits<float>::max()) const
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
		float c = dot(oPos, oPos) - radius * radius;
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
};

int main()
{
	Image image;
	image.init(640,480);

	Camera camera;

	Sphere sphere(vec3(0, 0, 3), 1);

	float aspectRatio = (float)image.w() / (float)image.h();

	// iterate over all image pixels in row-major order
	// for each image pixel shoot a ray and check for an intersection, set white if there's an intersection, and black if there's none
	for (u32 y = 0; y < image.h(); ++y)
	{
		for (u32 x = 0; x < image.w(); ++x)
		{
			// map [0,width]x[0,height] to [-aspectRatio,aspectRatio] x [1,-1]
			// multiply by the aspect ratio to stretch/squeeze the virtual film size to match the screen's aspect ratio
			float u = aspectRatio * (2.0f * ((float)x + 0.5f) / (float)image.w() - 1.0f);
			float v = -2.0f * ((float)y + 0.5f) / (float)image.h() + 1.0f;

			Ray ray = camera(u, v);

			image(x, y) = vec3(float(0.0f<sphere.intersect(ray)));
		}
	}
	image.savePPM("out.ppm");
	return 0;
}