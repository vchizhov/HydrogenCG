#pragma once
#include "vec.hpp"
#include "bsdf.hpp"

namespace HydrogenCG
{
	struct LightSample
	{
		vec3 o; //!< origin
		vec3 r;	//!< radiance
		

		LightSample() {}
		LightSample(vec3 o, vec3 r) : o(o), r(r) {}
	};

	struct Light
	{
		std::shared_ptr<Texture> i; 

		Light(std::shared_ptr<Texture> i = Texture::white) : i(i) {}
		Light(const vec3& color) : i(textureConstant(color)) {}
		/*!
			Given the poisition of the point to be shaded p, and uv coords
			returns the sampled point o on the light source, and the radiance r
			emitted from o in direction p-o
		*/
		virtual LightSample sample(const vec3& p, const vec2& uv = vec2(0)) const = 0;
	};

	struct LightDirectional : public Light
	{
		vec3 d; //!< direction - needs to be unit length
		
		LightDirectional() {}
		explicit LightDirectional(const vec3& d, std::shared_ptr<Texture> i = Texture::white) : Light(i),  d(normalize(d)) {}
		explicit LightDirectional(const vec3& d, const vec3& color) : Light(color), d(normalize(d)) {}

		virtual LightSample sample(const vec3& p, const vec2& uv) const
		{
			return LightSample(p+d, i->sample(uv));
		}
	};
	/*!
		\brief A basic point light structure.
	*/
	struct LightPoint : public Light
	{
		vec3 o;		//!< Position in 3d space of the point light

		LightPoint() {}
		explicit LightPoint(const vec3& origin, std::shared_ptr<Texture> i = Texture::white) : Light(i), o(origin) {}
		explicit LightPoint(const vec3& origin, const vec3& color) : Light(color), o(origin) {}

		virtual LightSample sample(const vec3& p, const vec2& /*uv*/) const
		{
			vec3 d = normalize(p - o);
			return LightSample(o, i->sample(vec2(atan2f(d.y,d.x), acosf(d.z))));
		}
	};

	/*struct LightArea : public Light
	{
		std::shared_ptr<SurfaceInstance> surf;
	};*/
}