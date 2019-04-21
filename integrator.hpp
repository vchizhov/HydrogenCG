#pragma once
#include "scene.hpp"
#include "image.hpp"
#include "camera.hpp"

namespace HydrogenCG
{
	/*!
		\brief	Base class for all integrators.

		The purpose of this class and inheriting classes are to render a scene from
		the perspective of a specific camera into an image.
		It also defines a function that returns the light energy (radiance) arriving at a ray's
		origin from the direction of the ray.
	*/
	class Integrator
	{
	public:
		virtual ~Integrator() {}

		/*!
			\brief				Renders the scene from the camera's perspective into the image.

			\param[out]			image Image object to render into.
			\param[in]			camera Camera from which to generate the rays for rendering the scene
			\param[in]			scene The scene to be rendered

			\return				A ray passing through the point corresponding to (u,v) on the virtual film
		*/
		virtual void render(Image& image, const Camera& camera, const Scene& scene) const
		{
			float aspectRatio = (float)image.w() / (float)image.h();
			for (uint32_t y = 0; y < image.h(); ++y)
			{
				for (uint32_t x = 0; x < image.w(); ++x)
				{
					// map [0,width]x[0,height] to [-aspectRatio,aspectRatio] x [1,-1]
					// multiply by the aspect ratio to non-uniformly stretch/squeeze the virtual film size to match the screen's aspect ratio
					float u = aspectRatio * (2.0f * ((float)x + 0.5f) / (float)image.w() - 1.0f);
					float v = -2.0f * ((float)y + 0.5f) / (float)image.h() + 1.0f;

					// generate ray corresponding to the normalized screen coordinates
					Ray ray = camera(u, v);

					// evaluate radiance arriving along the ray
					image(x, y) = radiance(ray, scene);
				}
			}
		}

		/*!
			\brief	Computes the radiance arriving from the scene along the ray direction.

			In this case a simple binary color returning white for intersections, and black for no intersections.
		*/
		virtual vec3 radiance(const Ray& ray, const Scene& scene) const
		{
			return vec3(static_cast<bool>(scene(ray)));
		}
	};

	/*!
		\brief	Renders reciprocal depth
	*/
	class IntegratorDepth : public Integrator
	{
	public:
		vec3 radiance(const Ray& ray, const Scene& scene) const final
		{
			return vec3(1 / scene(ray).t);
		}
	};

	/*!
		\brief	Renders normals
	*/
	class IntegratorNormal : public Integrator
	{
	public:
		vec3 radiance(const Ray& ray, const Scene& scene) const final
		{
			Intersection intersection = scene(ray);
			vec3 col = vec3(0); // black background
			if (intersection)
			{
				vec3 intersectionPoint = ray(intersection.t);
				// maps normals from [-1,1]^3 to [0,1]^3 (xyz -> rgb)
				col = 0.5f*intersection.s->normal(intersectionPoint) + vec3(0.5f);
			}
			return col;
		}
	};

	/*!
		\brief	Renders flat colors
	*/
	class IntegratorColor : public Integrator
	{
	public:
		vec3 radiance(const Ray& ray, const Scene& scene) const final
		{
			Intersection intersection = scene(ray);
			vec3 col = vec3(0); // black background
			if (intersection)
			{
				col = intersection.s->col;
			}
			return col;
		}
	};

	/*!
		\brief	Renders the scene using direct (local) illumination from point lights.

		Not unlike how rasterized graphics work, it ignores occlusion between the light source
		and the point to be illuminated.
	*/
	class IntegratorDiffuseLocal : public Integrator
	{
	public:
		vec3 radiance(const Ray& ray, const Scene& scene) const final
		{

			vec3 col = vec3(0); // black background
			Intersection intersection = scene(ray);
			if (intersection)
			{
				vec3 intersectionPoint = ray(intersection.t);
				vec3 normal = intersection.s->normal(intersectionPoint);
				vec3 intersectionColor = intersection.s->col;
				// iterate over all lights and add their contribution
				for (auto light : scene.lights)
				{
					// vector from the point being shaded to the light
					vec3 isectL = light.o - intersectionPoint;

					// squared distance between the point being shaded and the light
					// see: https://en.wikipedia.org/wiki/Inverse-square_law
					float rSquared = lengthSquared(isectL);
					float r = sqrtf(rSquared);
					isectL /= r; // normalize vector to light

					// account for illumination on the other side of the surface
					float cosThetaRay = dot(normal, ray.d);
					vec3 normalPointingTowardsRay = cosThetaRay <= 0.0f ? normal : -normal;

					// Lambert's law, the max is there to disallow negative energy accumulation
					// see: https://en.wikipedia.org/wiki/Lambert%27s_cosine_law
					float cosTheta = max(0.0f, dot(isectL, normalPointingTowardsRay));

					// diffuse material evaluation
					col += light.i * intersectionColor * cosTheta / rSquared;
				}
			}
			return col;
		}
	};

	/*!
	\brief	Renders transparent colors
*/
	class IntegratorTransparency : public Integrator
	{
	public:
		vec3 radiance(const Ray& ray, const Scene& scene) const final
		{

			Ray r = ray;
			vec3 col = vec3(1); // white background

			// allow no more than 10 intersections
			for (int i = 0; i < 10; ++i)
			{
				Intersection intersection = scene(r);
				if (!intersection) break;
				// use the color as a multiplicative filter (it absorbs the color spectrum which it is missing)
				col *= intersection.s->col;
				vec3 normal = intersection.s->normal(r(intersection.t)); // compute the normal to use in the offset to avoid self-intersection

				// flip the normal to account for a ray intersecting the surface from the other side
				float cosTheta = dot(normal, r.d);
				vec3 normalPointingTowardsRay = cosTheta <= 0.0f ? normal : -normal;

				r.o = r(intersection.t) - EPSILON * normalPointingTowardsRay; // regenerate ray starting from the last intersection
			}
			return col;
		}
	};

	/*!
		\brief	Renders the scene using direct illumination from point lights with shadows.

		Unlike rasterization graphics it accounts for the visibility function accurately and
		thus produces shadows.
	*/
	class IntegratorDiffuse : public Integrator
	{
	public:
		vec3 radiance(const Ray& ray, const Scene& scene) const final
		{

			vec3 col = vec3(0); // black background
			Intersection intersection = scene(ray);
			if (intersection)
			{
				vec3 intersectionPoint = ray(intersection.t);
				vec3 normal = intersection.s->normal(intersectionPoint);
				vec3 intersectionColor = intersection.s->col;
				// iterate over all lights and add their contribution
				for (auto light : scene.lights)
				{
					// vector from the intersection point to the light
					vec3 isectL = light.o - intersectionPoint;

					// account for illumination on the other side of the surface
					float cosThetaRay = dot(normal, ray.d);
					vec3 normalPointingTowardsRay = cosThetaRay <= 0.0f ? normal : -normal;

					// squared distance between the point being shaded and the light
					// see: https://en.wikipedia.org/wiki/Inverse-square_law
					float rSquared = lengthSquared(isectL);
					float r = sqrtf(rSquared);
					isectL /= r; // normalize shading vector

					// Lambert's law
					// see: https://en.wikipedia.org/wiki/Lambert%27s_cosine_law
					float cosTheta = dot(isectL, normalPointingTowardsRay);

					// if the contribution will be 0 anyways, skip the shadow ray casting below
					if (cosTheta > 0.0f)
					{
						// avoid shadow ray self intersection by an appropriate offset along the normal in the direction of the light
						float cosThetaL = dot(isectL, normal);
						vec3 offsetRayOrigin = intersectionPoint + EPSILON * (cosThetaL >= 0 ? normal : -normal);

						// add light contribution only if there's no surface between the point being shaded and the light
						if (!scene.intersectAny(Ray(offsetRayOrigin, isectL), 0, r - EPSILON))
						{
							// diffuse material evaluation
							// no need for max because of the if check
							col += light.i * intersectionColor * cosTheta / rSquared;
						}
					}
				}
			}
			return col;
		}
	};

}