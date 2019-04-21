#pragma once
#include "vec.hpp"
#include "ray.hpp"

namespace HydrogenCG
{
	/*!
		\brief A very simple pinhole camera class.
	*/
	struct Camera
	{
		vec3 o;			//!< The position of the the aperture of the camera
		vec3 e0;		//!< The first basis vector for the camera, "right" in our case
		vec3 e1;		//!< The second basis vector for the camera, "up" in our case
		vec3 e2;		//!< The third basis vector for the camera, "forward" in our case


		/*!
			\brief			Initializes the camera

			The vectors right, up, and forward should not be coplanar,
			usually they form an orthonormal basis (orthogonal + unit length).
		*/
		explicit Camera(const vec3& origin = vec3(0), const vec3& right = vec3(1.0f, 0.0f, 0.0f),
			const vec3& up = vec3(0.0f, 1.0f, 0.0f), const vec3& forward = vec3(0.0f, 0.0f, 1.0f))
			: o(origin), e0(right), e1(up), e2(forward) {}


		/*!
			\brief			Given normalized screen coordinates returns a ray passing through that point on the virtual film

			\param[in]		u The "x" coordinate on the film, is measured on the e0 basis vector
			\param[in]		v The "y" coordinate on the film, is measured on the e1 basis vector

			\return			A ray passing through the point corresponding to (u,v) on the virtual film
		*/
		inline Ray generateRay(float u, float v) const
		{
			return Ray(o, normalize(u * e0 + v * e1 + e2));
		}

		// convenience/brevity
		inline Ray operator()(float u, float v) const
		{
			return generateRay(u, v);
		}
	};
}