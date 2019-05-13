#pragma once

#include "../core/surface.hpp"
#include "../core/surface_instance.hpp"

namespace HydrogenCG
{
	/*!
		\brief	A box described as through its minimum and maximum coordinate vertices.
	*/
	struct Box : public Surface
	{
		vec3 minV;	//!< the vertex with minimum coordinates
		vec3 maxV;	//!< the vertex with maximum coordinates

		Box() {}
		explicit Box(const vec3& minVertex, const vec3& maxVertex)
			: minV(minVertex), maxV(maxVertex) {}

		/*!
			\brief		Box intersection
		*/
		Intersection intersect(const Ray& ray, float minT = 0.0f, float maxT = std::numeric_limits<float>::infinity()) const final
		{
			/*
				Box intersection derivation:
				A box can be represented as the set intersection of 6 axis-aligned half-spaces (defined through the 6 faces of the box)

				Let us consider one of those, for example the plane given by: x = minV.x
				This is the plane equation of the left face of the box in Hesse normal form, to see 
				that it is indeed the case one can rewrite it as:
				x = minV.x
				1 * (x-minV.x) = 0
				1 * (x - minV.x) + 0 * (y-minV.y) + 0 * (z-minV.z) = 0, or equivalently using the definition of the inner-product:
				<(x,y,z) - minV,(1,0,0)> = 0, thus it represent a plane with normal (1,0,0)
				and a point on it with coordinates minV, which is precisely the leftmost plane defining the box.

				We can show the same thing for the other 5 planes.

				To find an intersection with an axis aligned plane we need only consider the coordinate for which 
				the normal is non-zero. In the case above we have:
				<o+t_x*d - minV,(1,0,0)> = 0
				o.x+t_x * d.x - minV.x = 0
				t_y = (minV.x - o.x)/d.x

				We can find the other 5 intersections in the same way.
				One can concisely summarize the intersections with 
				the 3 planes defined through minV, through the vector expression
				t0 = (minV-o)/d
				Where all operations are to be read as coordinatewise and t0.x gives 
				the distance along the ray to the intersection with the left plane,
				t0.y with the bottom plane, t0.z with the "back" plane of the box

				A similar expression can be written for the right, top, and front planes
				t1 = (maxV - o)/d

				This produces 3 intervals [t0.x,t1.x], [t0.y,t1.y], [t0.z,t1.z], since it is possible that
				t0.x>t1.x for example, we can swap the t values so that t0 < t1:
				tmin = min(t0,t1)
				tmax = max(t0,t1)
				Where min and max act coordinatewise

				Finally we know that the box is defined as the intersection of all 6 halfspaces, thus we need to 
				find the set intersection of the intervals [tmin.x,tmax.x], [tmin.y,tmax.y], [tmin.z,tmax.z]
				this can be achieved in a concise manner by taking the maximum component of the minimums
				and the minimum component of the maximums:
				tN = maxComp(tmin);
				tF = minComp(tmax);
				
				if the interval [tN,tF] is empty (tF<=tN) then there the ray misses the box,
				otherwise we pick the closest intersection withing the bounds [minT,maxT]
			*/
			vec3 invRaydir = 1.0f / ray.d;
			// intersection with the far and near plane for each coordinate
			vec3 t0 = (minV - ray.o) * invRaydir;
			vec3 t1 = (maxV - ray.o) * invRaydir;
			// swap intersection interval ends such that t0<=t1
			vec3 tmin = min(t0, t1), tmax = max(t0, t1);
			// find the closest actual intersection (it needs to be the closest distance withing the set intersection of all intervals)
			float tN = maxComp(tmin);
			float tF = minComp(tmax);
			// if the set interval intersection is empty -> no collision
			if (tN >= tF) return Intersection();
			// otherwise verify that it is in the required bounds
			if (minT < tN && tN < maxT)
			{
				vec3 normal;
				vec2 uv;
				normalAndUV(ray(tN), normal, uv);
				return Intersection(tN, normal, uv);
			}
			// similarly for the far intersection (in case we are inside the box)
			if (minT < tF && tF < maxT)
			{
				vec3 normal;
				vec2 uv;
				normalAndUV(ray(tF), normal, uv);
				return Intersection(tF, normal, uv);
			}
			return Intersection();
		}

		/*!
			\brief	Returns the closest distance from a point p to the box
		*/
		float distance(const vec3& p) const final
		{
			// https://iquilezles.org/www/articles/distfunctions/distfunctions.htm
			vec3 boxOrigin = 0.5f*(maxV + minV);
			vec3 boxScale = 0.5f*(maxV - minV);

			
			// handle only the positive octant case and reduce all other cases to it
			// by taking the absolute value of the coordinates transformed into the local 
			// coordinate system of the box (as if it is centered at (0,0,0))
			vec3 centeredAbsPoint = abs(p - boxOrigin);

			// compute the distance between the point and the cube's planes along each coordinate separately
			// we can do this since the cube is axis aligned, so we can just compute the distance along each coordinate
			vec3 d = centeredAbsPoint - boxScale;

			// if not all coordinates of d are negative, then the point is outside of the box, since
			// it does not lie withing the set intersection all half-spaces defining the cube
			// simply compute the euclidean distance to the cube
			float distanceOutside = length(max(d, 0.0f));

			// if all coordinates of d are negative, then the point is inside of the box, and 
			// the closest distance is the least negative one - that is, the greatest coordinate:
			float distanceInside = min(maxComp(d), 0.0f);

			// the point cannot be both inside and outside (unless it's at the boundary, but then the distance is 0)
			// thus the two cases can never happen at the same time, so simply return the sum (since at least one of the 
			// distances will always be 0)
			return distanceOutside + distanceInside;
		}

		/*!
			\brief		Computes the sphere normal for a point p

			\param[in]	p A point in space

			\return		Returns the outwards facing unit normal at point p
						if p is on the surface of the sphere, otherwise
						returns a scaled version of the normal
		*/
		void normalAndUV(const vec3& p, vec3& normal, vec2& uv) const
		{
			// similar to the distance computation
			vec3 boxOrigin = 0.5f*(maxV + minV);
			vec3 boxScale = 0.5f*(maxV - minV);
			vec3 d = p - boxOrigin;
			vec3 absD = abs(d)-boxScale;
			// compute uvw coordinatewise
			vec3 uvw = (p - minV) / (maxV - minV);
			if (absD.x >= absD.y && absD.x >= absD.z)
			{
				normal = vec3(sign(d.x), 0, 0);
				uv = vec2(uvw.y, uvw.z);
			}
			else if (absD.y > absD.x && absD.y >= absD.z)
			{
				normal = vec3(0, sign(d.y), 0);
				uv = vec2(uvw.x, uvw.z);
			}
			else //if (absD.z > absD.x && absD.z > absD.y)
			{
				normal = vec3(0, 0, sign(d.z));
				uv = vec2(uvw.x, uvw.y);
			}
		}

		static const std::shared_ptr<const Box> unitBox;
	};
	// convenience/brevity
	const std::shared_ptr<const Box> Box::unitBox = std::make_shared<const Box>(vec3(-1), vec3(1));

	std::shared_ptr<SurfaceInstanceBSDF> box()
	{
		return std::make_shared<SurfaceInstanceBSDF>(Box::unitBox);
	}

	std::shared_ptr<SurfaceInstanceBSDF> box(const vec3& minVertex, const vec3& maxVertex, const vec3& col = vec3(1))
	{
		return std::make_shared<SurfaceInstanceBSDF>(std::make_shared<Box>(minVertex, maxVertex), bsdfDiffuse(textureConstant(col)));
	}
}