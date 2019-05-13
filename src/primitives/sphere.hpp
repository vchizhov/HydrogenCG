#pragma once
#include "..\core\surface.hpp"
#include "..\core\surface_instance.hpp"

namespace HydrogenCG
{
	/*!
		\brief	A sphere described through the canonical equation ||p-o|| = r
	*/
	struct Sphere : public Surface
	{
		vec3 o;		//!< the origin/center of the sphere
		float r;	//!< the radius of the sphere

		Sphere() {}
		explicit Sphere(const vec3& origin, float radius) : o(origin), r(radius) {}

		/*!
			\brief		Sphere intersection
		*/
		Intersection intersect(const Ray& ray, float minT = 0.0f, float maxT = std::numeric_limits<float>::infinity()) const final
		{
			/*
				Sphere intersection derivation:
				||x|| is used to denote the euclidean norm of vector x
				<a,b> is the dot(inner) product of vectors a and b
				r^2 = r*r

				||ray(t) - pos|| == r <->
				||ray(t) - pos||^2 == r^2 <->
				<ray(t)-o,ray(t)-o> == r^2 <->
				<ray.o-o + t * ray.d, ray.o-o + t * ray.d> == r^2 <->
				<ray.d,ray.d> * t^2 - 2 * <ray.d, o - ray.o> * t + <o-ray.o,o-ray.o> - r^2 == 0
				A = <ray.d,ray.d>, but |ray.d|==1 if the direction is normalized -> A = 1
				B = <ray.d, o - ray.o>
				C = <o-ray.o,o-ray.o> - r^2

				D = B^2 - C
				D<0  -> no intersection
				D==0 -> grazing intersection
				D>0  -> 2 intersections

				One can ignore grazing intersections since they are actually numerical error.

				D>0 -> sqrtD = sqrt(D)
				t_1 = B - D
				t_2 = B + D

				We want to return the closer intersection in the range (minT,maxT)
			*/

			// vector from the ray origin to the sphere origin
			vec3 ro = o - ray.o;
			// quadratic polynomial terms (a==1)
			float b = dot(ray.d, ro);
			float c = dot(ro, ro) - r * r;
			// discriminant
			float discriminant = b * b - c;

			// If the discriminant is 0 or negative -> no (real) intersection
			if (discriminant <= 0.0f) return Intersection();

			float sqrtD = sqrt(discriminant);
			float tN = b - sqrtD; // closer intersection (where the ray enters the sphere)
			// if it is within the defined ray segment by (minT,maxT) return it as the closest intersection
			if (tN > minT && tN < maxT)
			{
				vec3 normal;
				vec2 uv;
				normalAndUV(ray(tN), normal, uv);
				return Intersection(tN, normal, uv);
			}

			float tF = b + sqrtD; // further intersection (where the ray exits the sphere)
			// if it is within the defined ray segment by (minT,maxT) return it as the closest intersection
			if (tF > minT && tF < maxT)
			{
				vec3 normal;
				vec2 uv;
				normalAndUV(ray(tF), normal, uv);
				return Intersection(tF, normal, uv);
			}

			// otherwise return no intersection
			return Intersection();
		}


		/*!
			\brief		Computes the sphere normal for a point p and its uv coordinates

			\param[in]	p A point in space

			\return		Returns the outwards facing unit normal at point p
						if p is on the surface of the sphere, otherwise
						returns a scaled version of the normal
		*/
		void normalAndUV(const vec3& p, vec3& normal, vec2& uv) const
		{
			normal = (p - o) / r;
			uv = vec2(atan2f(normal.y, normal.x), acosf(normal.z));
		}

		float distance(const vec3& p) const final
		{
			// the (signed) distance from a point to a sphere is simply the distance
			// from the point to the sphere center minus the radius of the sphere
			return length(p - o) - r;
		}

		static const std::shared_ptr<const Sphere> unitSphere;
	};
	const std::shared_ptr<const Sphere> Sphere::unitSphere = std::make_shared<const Sphere>(vec3(0.0f), 1.0f);
	std::shared_ptr<SurfaceInstanceBSDF> sphere()
	{
		return std::make_shared<SurfaceInstanceBSDF>(Sphere::unitSphere);
	}

	std::shared_ptr<SurfaceInstanceBSDF> sphere(const vec3& origin, float radius, const vec3& col = vec3(1))
	{
		return std::make_shared<SurfaceInstanceBSDF>(std::make_shared<Sphere>(origin, radius), bsdfDiffuse(textureConstant(col)));
	}
}