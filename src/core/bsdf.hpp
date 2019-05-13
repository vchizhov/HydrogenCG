#pragma once
#include "texture.hpp"
namespace HydrogenCG
{

	/*!
		A struct used as an argument for bsdf shading
	*/
	struct SurfacePoint
	{
		SurfacePoint() {}
		SurfacePoint(const vec3& normal, const vec2& uv) : normal(normal), uv(uv) {}
		vec3 normal;	//!< normal of the intersection
		vec2 uv;		//!< uv coordinates of the intersection
	};
	struct BSDFDiffuse;

	/*
		Bidirectional scattering distribution function (BSDF) structure
		(a function describing the physical properties of a material in the context of light interaction)

		Provides a way to evaluate a bsdf given a surface point and an incident and outgoing direction
	*/
	struct BSDF
	{
		virtual vec3 eval(const SurfacePoint& p, const vec3& out, const vec3& in) const = 0;

		static const std::shared_ptr<BSDF> white;

		static std::shared_ptr<BSDFDiffuse> diffuse(std::shared_ptr<Texture> albedo)
		{
			return std::make_unique<BSDFDiffuse>(albedo);
		}
	};

	struct BSDFDiffuse : public BSDF
	{
		std::shared_ptr<Texture> albedo;	//!< No component should be greater than 1 for energy conservation to hold
		BSDFDiffuse(std::shared_ptr<Texture> albedo) : albedo(albedo) {}
		vec3 eval(const SurfacePoint& p, const vec3& /*out*/, const vec3& /*in*/) const final
		{
			return albedo->sample(p.uv) / pi;
		}
	};

	std::shared_ptr<BSDF> bsdfDiffuse(std::shared_ptr<Texture> albedo)
	{
		return std::make_unique<BSDFDiffuse>(albedo);
	}

	std::shared_ptr<BSDFDiffuse> bsdfDiffuse(const vec3& color)
	{
		return std::make_unique<BSDFDiffuse>(textureConstant(color));
	}

	const std::shared_ptr<BSDF> BSDF::white = std::make_unique<BSDFDiffuse>(Texture::white);
}