#pragma once
#include "scene.h"
#include "image.h"
#include "camera.h"

namespace HydrogenCG
{
	class Integrator
	{
	public:
		virtual void render(Image& image, const Camera& camera, const Scene& scene) const
		{
			float aspectRatio = (float)image.w() / (float)image.h();
			for (uint32_t y = 0; y < image.h(); ++y)
			{
				for (uint32_t x = 0; x < image.w(); ++x)
				{
					// map [0,width]x[0,height] to [-aspectRatio,aspectRatio] x [1,-1]
					// multiply by the aspect ratio to non-uniformly stretch/squeeze the virtual film size to match the screen's aspect ratio
					float u = aspectRatio * (2.0f * ((float)x+0.5f) / (float)image.w() - 1.0f);
					float v = -2.0f * ((float)y + 0.5f) / (float)image.h() + 1.0f;

					Ray ray = camera(u, v);

					// evaluate radiance
					image(x, y) = li(ray, scene);
				}
			}
		}

		virtual vec3 li(const Ray& ray, const Scene& scene) const
		{
			IntersectInfo intersect = scene.intersect(ray);
			return vec3(intersect);
		}
	};

	class IntegratorDepth : public Integrator
	{
	public:
		virtual vec3 li(const Ray& ray, const Scene& scene) const override
		{
			float farPlane = 0.0f;
			IntersectInfo intersect = scene.intersect(ray);
			return vec3(1 / intersect.t);
		}
	};

	class IntegratorNormal : public Integrator
	{
	public:
		virtual vec3 li(const Ray& ray, const Scene& scene) const override
		{
			float farPlane = 0.0f;
			IntersectInfo intersect = scene.intersect(ray);
			vec3 col = vec3(0);
			if (intersect)
			{
				const Surface* s = scene(intersect.idx);
				vec3 intersectionPoint = ray(intersect.t);
				col = s->normal(intersectionPoint);
				col = 0.5f*col + vec3(0.5f);
			}
			return col;
		}
	};

	class IntegratorColor : public Integrator
	{
	public:
		virtual vec3 li(const Ray& ray, const Scene& scene) const override
		{
			float farPlane = 0.0f;
			IntersectInfo intersect = scene.intersect(ray);
			vec3 col = vec3(0);
			if (intersect)
			{
				col = scene(intersect.idx)->col;
			}
			return col;
		}
	};

	// transparent
	class IntegratorTransparent : public Integrator
	{
	public:
		virtual vec3 li(const Ray& ray, const Scene& scene) const override
		{
			
			Ray r = ray;
			vec3 col = vec3(1.0f);
			for (int i = 0; i < 10; ++i)
			{
				IntersectInfo intersect = scene.intersect(r);
				if (!intersect) break;
				col *= scene(intersect.idx)->col;
				r.o = r(intersect.t);
			}
			return col;
		}
	};

	// point light
	class IntegratorDiffuseLocal : public Integrator
	{
	public:
		virtual vec3 li(const Ray& ray, const Scene& scene) const override
		{

			vec3 col = vec3(0.0f);
			IntersectInfo intersect = scene.intersect(ray);
			if (intersect)
			{
				vec3 intersectionPoint = ray(intersect.t);
				vec3 normal = scene.surfaces[intersect.idx]->normal(intersectionPoint);
				vec3 intersectionColor = scene.surfaces[intersect.idx]->col;
				for (auto light : scene.lights)
				{
					float rSquared = lengthSquared(light.pos - intersectionPoint);
					float cosTheta = max(0.0f, dot(-ray.d, normal));
					col += light.intensity * intersectionColor * cosTheta / rSquared;
				}
			}
			return col;
		}
	};

}