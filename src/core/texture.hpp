#pragma once
#include "vec.hpp"
#include <memory>
#include "image.hpp"
namespace HydrogenCG
{

	/*!
		A base class for textures
	*/
	struct Texture
	{
		virtual vec3 sample(const vec2& uv) const = 0;

		static const std::shared_ptr<Texture> white;
	};

	/*!
		A constant (only one color everywhere) texture
	*/
	struct TextureConstant : public Texture
	{
		vec3 color;	//!< the color of the texture
		TextureConstant(const vec3& color) : color(color) {}
		vec3 sample(const vec2& /*uv*/) const final { return color; }
	};

	std::shared_ptr<Texture> textureConstant(const vec3& color)
	{
		return std::make_unique<TextureConstant>(color);
	}

	const std::shared_ptr<Texture> Texture::white = std::make_unique<TextureConstant>(vec3(1.0f));

	/*!
		A checkerboard texture
	*/
	struct TextureCheckerboard : public Texture
	{
		vec3 color1, color2;	//!< the two colours for the checkerboard
		float uFreq, vFreq;		//!< the frequency (the larger the frequency the smaller the regions) of the pattern
		float uOffset, vOffset;	//!< the initial offset of the pattern
		TextureCheckerboard(const vec3& color1, const vec3& color2,
			float uFreq = 1.0f, float vFreq = 1.0f, float uOffset = 0.0f, float vOffset = 0.0f)
			: color1(color1), color2(color2), uFreq(uFreq), vFreq(vFreq), uOffset(uOffset), vOffset(vOffset) {}
		vec3 sample(const vec2& uv) const final
		{
			int x = static_cast<int>(floorf(uv.x*uFreq + uOffset));
			int y = static_cast<int>(floorf(uv.y*vFreq + vOffset));
			// if both are even or both are odd -> color1
			return ((x % 2 == 0) == (y % 2 == 0)) ? color1 : color2;
		}
	};

	std::shared_ptr<Texture> textureCheckerboard(const vec3& color1, const vec3& color2, float uFreq = 1.0f, float vFreq = 1.0f, float uOffset = 0.0f, float vOffset = 0.0f)
	{
		return std::make_unique<TextureCheckerboard>(color1, color2, uFreq, vFreq, uOffset, vOffset);
	}

	/*!
		An image texture aiming to provide similar sampling options to OpenGL textures
	*/
	struct TextureImage : public Texture
	{
		enum TEXTURE_SAMPLING
		{
			TEXTURE_SAMPLING_NEAREST,
			TEXTURE_SAMPLING_BILINEAR,
			TEXTURE_SAMPLING_BICUBIC // left as an exercise for interested readers
		};

		enum TEXTURE_WRAPPING
		{
			TEXTURE_WRAPPING_CLAMP,
			TEXTURE_WRAPPING_REPEAT,
			TEXTURE_WRAPPING_MIRROR,
			TEXTURE_WRAPPING_BORDER
		};

		enum TEXTURE_INTERPOLATION
		{
			TEXTURE_INTERPOLATION_LINEAR,
			TEXTURE_INTERPOLATION_CUBIC,
			TEXTURE_INTERPOLATION_QUINTIC
		};

		std::shared_ptr<Image> image;
		vec3 borderColor;
		TEXTURE_SAMPLING sampling;
		TEXTURE_WRAPPING wrapping;
		TEXTURE_INTERPOLATION interpolation;
		TextureImage(std::shared_ptr<Image> image,
			TEXTURE_SAMPLING sampling = TEXTURE_SAMPLING_NEAREST,
			TEXTURE_WRAPPING wrapping = TEXTURE_WRAPPING_CLAMP,
			TEXTURE_INTERPOLATION interpolation = TEXTURE_INTERPOLATION_LINEAR,
			const vec3& borderColor = vec3(0))
			: image(image), sampling(sampling), wrapping(wrapping),
			interpolation(interpolation), borderColor(borderColor) {}

		// expects the whole and fraction parts of (u*w,y*h)
		vec3 samplePixel(float wx, float wy, float fx, float fy, TEXTURE_WRAPPING wrapMode) const
		{
			float x, y;
			if (wrapMode == TEXTURE_WRAPPING_CLAMP)
			{
				// from pixel centered coordinates we go to pixel corner coordinates
				// that's why we directly use the fwhole part, every interval is floored,
				// that way each pixel gets an equal interval: 0: [0,1), 1:[1,2), ..., w-1: [w-1,w)
				// so the [0,w] range is split evenly, similarly for [0,h]
				x = wx;
				y = wy;

				// clamp from below to disallow negative pixel coordinates
				x = max(0.0f, x);
				y = max(0.0f, y);
			}
			else if (wrapMode == TEXTURE_WRAPPING_REPEAT)
			{
				// tile the image with [0,1) intervals
				x = fx;
				y = fy;

				// map [0,1)x[0,1) -> [0,w)x[0,h)
				x = floorf(x*image->w());
				y = floorf(y*image->h());
			}
			else if (wrapMode == TEXTURE_WRAPPING_MIRROR)
			{
				// if the whole part is even use non-mirrored image,
				// if the whole part is odd mirror the image
				x = (static_cast<u32>(fx) % 2 == 0) ? fx : 1.0f - fx;
				y = (static_cast<u32>(fy) % 2 == 0) ? fy : 1.0f - fy;

				// map [0,1]x[0,1] -> [0,w]x[0,h]
				x = floorf(x*image->w());
				y = floorf(y*image->h());
			}
			else //if (wrapping == TEXTURE_WRAPPING_BORDER)
			{
				// the same as clamp, the only difference is that we return the border color if outside of bounds
				x = wx;
				y = wy;
				// clamp to border color if outside of [0,w)x[0,h)
				if (x < 0.0f || y < 0.0f || x >= image->w() || y >= image->h())
					return borderColor;
			}
			// clamp from above just in case
			u32 px = min(static_cast<u32>(x), image->w() - 1);
			u32 py = min(static_cast<u32>(y), image->h() - 1);
			return image->get(px, py);
		}

		vec3 sample(const vec2& uv) const final
		{
			// Note [a,b] and [a,b) should be considered equivalent here, since a single point has length zero. 
			// In practice this is not the case, since precision is finite, but we cannot do any 
			// better anyways (you will always have to bias the remaining bit either upwards or downwards)

			// (u,v) are normalized pixel centered coordinates -> map to image coordinates: [0,1]x[0,1] -> [0,w]x[0,h]
			float x = uv.x * image->w();
			float y = uv.y * image->h();

			// compute whole part
			float wx = floorf(x);
			float wy = floorf(y);
			// compute the fractional part, in [0,1)x[0,1)
			float fx = x - wx;
			float fy = y - wy;

			vec3 color;
			if (sampling == TEXTURE_SAMPLING_NEAREST) // uses 1 samples
			{
				color = samplePixel(wx, wy, fx, fy, wrapping);
			}
			else if (sampling == TEXTURE_SAMPLING_BILINEAR) // uses 4 samples
			{
				if (interpolation == TEXTURE_INTERPOLATION_CUBIC)
				{
					fx = smoothstepCubic(fx);
					fy = smoothstepCubic(fy);
				}
				else if (interpolation == TEXTURE_INTERPOLATION_QUINTIC)
				{
					fx = smoothstepQuintic(fx);
					fy = smoothstepQuintic(fy);
				}
				//else if (interpolation == TEXTURE_INTERPOLATION_LINEAR)
				//{
				//	// do nothing
				//}

				vec3 color00 = samplePixel(wx, wy, 0, 0, wrapping);
				vec3 color10 = samplePixel(wx + 1.0f, wy, 0, 0, wrapping);
				vec3 color01 = samplePixel(wx, wy + 1.0f, 0, 0, wrapping);
				vec3 color11 = samplePixel(wx + 1.0f, wy + 1.0f, 0, 0, wrapping);
				color = blerp(color00, color10, color01, color11, fx, fy);
			}
			else if (sampling == TEXTURE_SAMPLING_BICUBIC) // left as an exercise for the interested reader
			{
				// refer to: https://en.wikipedia.org/wiki/Bicubic_interpolation
				color = vec3(0);
			}
			return color;
		}
	};
}