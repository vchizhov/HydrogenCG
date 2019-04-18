#pragma once
#include "sphere.h"
#include "light.h"
#include <vector>

namespace HydrogenCG
{
	struct IntersectInfo
	{
		float t;
		int idx;
		IntersectInfo() : idx(-1), t(-1.0f) {}

		operator bool() const
		{
			return idx > -1;
		}

		/*operator int() const
		{
			return idx;
		}*/
	};

	struct Scene
	{
		std::vector<Surface*> surfaces;
		std::vector<LightPoint> lights;

		~Scene() 
		{
			for (int i = 0; i < surfaces.size(); ++i) delete surfaces[i]; 
		}

		Surface* operator()(uint32_t i) { return surfaces[i]; }
		const Surface* operator()(uint32_t i) const { return surfaces[i]; }

		IntersectInfo intersect(const Ray& ray, float minT = 0.0f, float maxT = std::numeric_limits<float>::max()) const
		{
			IntersectInfo info;
			// initialize distance to draw distance
			info.t = maxT;
			// iterate over all objects
			for (int i = 0; i < surfaces.size(); ++i)
			{
				float t = surfaces[i]->intersect(ray, minT, info.t);
				if (t > minT)
				{
					info.idx = (int)i;
					info.t = t;
				}
			}

			return info;
		}
	};
}