#pragma once
#include "surface.hpp"
#include "mat.hpp"
#include <memory>
#include "bsdf.hpp"
namespace HydrogenCG
{
	struct InstanceIntersection : public Intersection
	{
		BSDF* bsdf;

		InstanceIntersection() {}
		InstanceIntersection(const Intersection& i, BSDF* b) : Intersection(i), bsdf(b) {}
	};
	 

	class SRTTransformation
	{
	protected:
		vec3 negativePos;	//!< The position coordinates negated (it's the inverse of the translation)
		mat3 transposedRot;	//!< The rotation matrix transposed (it's the inverse of the rotation)
		vec3 invScale;		//!< One divided by the scale vector

	public:
		SRTTransformation() : negativePos(0), transposedRot(mat3::identity), invScale(1) {}
		SRTTransformation(const vec3& pos, const mat3& rot = mat3::identity, const vec3& scale = vec3(1)) 
			: negativePos(-pos), transposedRot(transpose(rot)), invScale(1.0f/scale) {}

		SRTTransformation& translate(const vec3& translation)
		{
			negativePos -= translation;
			return *this;
		}

		// resets the position
		SRTTransformation& position(const vec3& position)
		{
			negativePos = -position;
			return *this;
		}
		// uniform scaling
		SRTTransformation& scale(float scale)
		{
			invScale = vec3(1.0f / scale);
			return *this;
		}
		// possibly non-uniform scale/reflection
		SRTTransformation& scale(const vec3& scale)
		{
			invScale = 1.0f / scale;
			return *this;
		}
		SRTTransformation& turn(const vec3& angles)
		{
			transposedRot *= transpose(rotationFromXYZ(angles));
			return *this;
		}
		// resets the rotation matrix
		SRTTransformation& rotate(const vec3& angles)
		{
			transposedRot = transpose(rotationFromXYZ(angles));
			return *this;
		}

		vec3 transformPoint(const vec3& p) const
		{
			// (TRS)^-1 * v = S^-1 * R^-1 * T^-1 * v = 1/S * tranpose(R) * (v-t)
			return invScale * (transposedRot*(p + negativePos));
		}

		vec3 transformVector(const vec3& v) const
		{
			// (RS)^-1 * v = S^-1 * R^-1 * v = 1/S * tranpose(R) * v
			return invScale * (transposedRot * v);
		}

		vec3 transformNormal(const vec3& n) const
		{
			// normalize(transpose((RS)^-1)*n) = normalize(R * 1/S * n) = scaleDet * R * (n/S)
			return transpose(transposedRot) * (invScale * n);
		}

		Ray transform(const Ray& ray) const
		{
			return Ray(transformPoint(ray.o), transformVector(ray.d));
		}
	};

	// note the determinant can be precomputed for static geometry
	struct SurfaceInstance : public SRTTransformation
	{
		SurfaceInstance() {}
		virtual ~SurfaceInstance() {}

		float transformRayAndBounds(Ray& ray, float& minT, float& maxT) const
		{
			// T * R * S * obj -> (T * R * S)^-1 * ray = S^-1 * R^-1 * T^-1 * ray = 1/S * tranpose(R) * (ray-t)
			ray = transform(ray);
			float invDirLength = 1.0f / length(ray.d);
			ray.d *= invDirLength;
			// if we make a larger ray smaller, distances become larger and vice versa
			minT /= invDirLength;
			maxT /= invDirLength;

			return invDirLength;
		}

		InstanceIntersection intersect(const Ray& ray, float minT = 0.0f, float maxT = std::numeric_limits<float>::infinity()) const
		{
			Ray transformedRay = ray;
			float rescaledMinT = minT;
			float rescaledMaxT = maxT;
			float invDirLength = transformRayAndBounds(transformedRay, rescaledMinT, rescaledMaxT);
			InstanceIntersection isect = intersectOriginal(transformedRay, rescaledMinT, rescaledMaxT);
			if (isect)
			{
				// normal transformation requires the tranpose of the inverse, in our case we have R*S
				// transpose((R*S)^-1) = tranpose(S^-1 * transpose(R)) = R * 1/S
				// finally we need to normalize the ray that we get like this, but we know the factor by which it was
				// lengthened/shortened, which is DirLength, so we just need to multiply by invDirLength
				// for uniform scaling scaleDet and invScale cancel out
				isect.n = invDirLength * transformNormal(isect.n);
				// rescale the distance traveled in the stretched/squeezed space
				isect.t *= invDirLength;
			}
			return isect;
		}

		virtual InstanceIntersection intersectOriginal(const Ray& /*ray*/, float /*minT*/ = 0.0f, float /*maxT*/ = std::numeric_limits<float>::infinity()) const = 0;

		bool intersectAny(const Ray& ray, float minT = 0.0f, float maxT = std::numeric_limits<float>::infinity()) const
		{
			Ray transformedRay = ray;
			float rescaledMinT = minT;
			float rescaledMaxT = maxT;
			transformRayAndBounds(transformedRay, rescaledMinT, rescaledMaxT);
			return intersectAnyOriginal(transformedRay, rescaledMinT, rescaledMaxT);
		}

		virtual bool intersectAnyOriginal(const Ray& ray, float minT = 0.0f, float maxT = std::numeric_limits<float>::infinity()) const
		{
			return static_cast<bool>(intersectOriginal(ray, minT, maxT));
		}



		virtual std::shared_ptr<SurfaceInstance> doClone() const = 0;

		// simulate covariance with smart pointers
		std::shared_ptr<SurfaceInstance> Clone() const
		{
			return doClone();
		}

		float distance(const vec3& p) const
		{
			return distanceOriginal(transformPoint(p));
		}

		virtual float distanceOriginal(const vec3& p) const = 0;
	};

	struct SurfaceInstanceBSDF : public SurfaceInstance
	{
		std::shared_ptr<const Surface> archetype;	//!< surface to instance - for example a 1 mil tri count model
		std::shared_ptr<BSDF> bsdf;					//!< the material used for the instance
		SurfaceInstanceBSDF() {}
		SurfaceInstanceBSDF(const std::shared_ptr<const Surface> surf, const std::shared_ptr<BSDF> bsdf = BSDF::white) : archetype(surf), bsdf(bsdf) {}
		SurfaceInstanceBSDF(const std::shared_ptr<const Surface> surf, const vec3& color) : archetype(surf), bsdf(bsdfDiffuse(color)) {}

		InstanceIntersection intersectOriginal(const Ray& ray, float minT = 0.0f, float maxT = std::numeric_limits<float>::infinity()) const final
		{
			return InstanceIntersection(archetype->intersect(ray, minT, maxT), bsdf.get());
		}



		std::shared_ptr<SurfaceInstance> doClone() const final
		{
			return std::make_shared<SurfaceInstanceBSDF>(*this);
		}

		// simulate covariance with smart pointers
		std::shared_ptr<SurfaceInstanceBSDF> Clone() const
		{
			return std::static_pointer_cast<SurfaceInstanceBSDF>(doClone());
		}

		float distanceOriginal(const vec3& p) const final
		{
			return archetype->distance(p);
		}
	};

	struct SurfaceInstanceCollection : public SurfaceInstance
	{
		std::vector<std::shared_ptr<SurfaceInstance>> surfaces;		//!< a collection of surfaces
		SurfaceInstanceCollection() {}

		// convenience
		std::shared_ptr<SurfaceInstance>& operator[](uint32_t i) { return surfaces[i]; }
		const std::shared_ptr<SurfaceInstance>& operator[](uint32_t i) const { return surfaces[i]; }

		// All of those take ownership of the instance
		// for l-values
		template<typename SurfaceInstanceDerived>
		SurfaceInstance& add(std::shared_ptr<SurfaceInstanceDerived>& s)
		{
			assert(s && "SurfaceInstanceCollection::add():: The pointer is null.");
			surfaces.push_back(std::static_pointer_cast<SurfaceInstance>(s));
			return *this;
		}
		// for r-values
		template<typename SurfaceInstanceDerived>
		SurfaceInstance& add(std::shared_ptr<SurfaceInstanceDerived>&& s)
		{
			assert(s && "SurfaceInstanceCollection::add():: The pointer is null.");
			surfaces.push_back(std::static_pointer_cast<SurfaceInstance>(s));
			return *this;
		}
		// for l-values
		SurfaceInstance& operator+=(std::shared_ptr<SurfaceInstance>& s)
		{
			assert(s && "SurfaceInstanceCollection::add():: The pointer is null.");
			surfaces.push_back(std::move(s));
			return *this;
		}
		// for r-values
		SurfaceInstance& operator+=(std::shared_ptr<SurfaceInstance>&& s)
		{
			assert(s && "SurfaceInstanceCollection::add():: The pointer is null.");
			surfaces.push_back(std::move(s));
			return *this;
		}

		InstanceIntersection intersectOriginal(const Ray& ray, float minT = 0.0f, float maxT = std::numeric_limits<float>::infinity()) const final
		{
			InstanceIntersection intersection;

			// initialize distance to maximum draw distance
			// this will be used as the closest intersection so far
			float closestT = maxT;

			// iterate over all surfaces and find the intersection with the closest one
			// what this algorithm essentially does is finding the minimum of an "array"
			// in this case the "array" is implicit and is defined by the sequence of values 
			// produced by calling intersect on each surface
			for (const auto& s : surfaces)
			{
				// we feed in the currently closest intersection as the upper bound,
				// thus if the intersection in the current iteration is not closer 
				// than intersection.t (the closest thus far), we would get no intersection
				// Since this is the case, the invariant that intersection.t contains the closest
				// intersection so far will be maintained
				InstanceIntersection tempI = s->intersect(ray, minT, closestT);

				// intersect returns maxT+1 if no intersection occurrs
				// so we use the check t<maxT as an indicator for an intersection
				if (tempI)
				{
					intersection = tempI;
					closestT = tempI.t;
				}
			}
			return intersection;
		}

		/*!
			\brief			Returns true at the first found intersection, false otherwise

			Useful for testing for visibility between two points,
			set maxT to the distance between those minus epsilon in that case.
		*/
		bool intersectAnyOriginal(const Ray& ray, float minT = 0.0f, float maxT = std::numeric_limits<float>::infinity()) const final
		{
			for (const auto& s : surfaces)
			{
				if (s->intersect(ray, minT, maxT)) return true;
			}
			return false;
		}
	
		std::shared_ptr<SurfaceInstance> doClone() const final
		{
			auto inst = std::make_shared<SurfaceInstanceCollection>();
			inst->negativePos = negativePos;
			inst->transposedRot = transposedRot;
			inst->invScale = invScale;
			inst->surfaces.resize(surfaces.size());
			for (size_t i=0;i<inst->surfaces.size();++i)
			{
				inst->surfaces[i] = surfaces[i]->Clone();
			}

			return inst;
		}

		// simulate covariance for smart pointers
		std::shared_ptr<SurfaceInstanceCollection> Clone() const
		{
			return std::static_pointer_cast<SurfaceInstanceCollection>(doClone());
		}
	
		float distanceOriginal(const vec3& p) const final
		{
			float dist = std::numeric_limits<float>::infinity();
			for (auto& e : surfaces) dist = min(dist, e->distance(p));
			return dist;
		}
	};

	std::shared_ptr<SurfaceInstanceCollection> collection()
	{
		return std::make_unique<SurfaceInstanceCollection>();
	}

	template<typename A>
	std::shared_ptr<A> operator+(const std::shared_ptr<A>& surf, const vec3& translation)
	{
		auto res = surf->Clone();
		res.get()->translate(translation);
		return res;
	}

	template<typename A>
	std::shared_ptr<A> operator*(const vec3& scale, std::shared_ptr<A> surf)
	{
		surf.get()->scale(scale);
		return surf;
	}

	template<typename A>
	std::shared_ptr<A> operator*(float scale, std::shared_ptr<A> surf)
	{
		surf.get()->scale(scale);
		return surf;
	}

	template<typename A>
	std::shared_ptr<A> turn(std::shared_ptr<A> surf, const vec3& angles)
	{
		surf.get()->turn(angles);
		return surf;
	}

	template<typename A>
	std::shared_ptr<A> rotate(std::shared_ptr<A> surf, const vec3& angles)
	{
		surf.get()->rotate(angles);
		return surf;
	}


	/*std::unique_ptr<SurfaceInstance> instance(std::shared_ptr<Surface> surf)
	{
		return std::make_unique<SurfaceInstance>(surf);
	}*/


	//std::unique_ptr<Surface> operator+(std::unique_ptr<Surface> lhs, const vec3& rhs)
	//{
	//	mat4& m = lhs.get()->transform;
	//	m.setColumn(3, m.getColumn(3) - vec4(rhs, 0));
	//	return lhs;
	//}

	//std::unique_ptr<Surface> operator+(const vec3& lhs, std::unique_ptr<Surface> rhs)
	//{
	//	mat4& m = rhs.get()->transform;
	//	m.setColumn(3, m.getColumn(3) - vec4(lhs, 0));
	//	return rhs;
	//}

	//std::unique_ptr<Surface> operator*(const mat4& lhs, std::unique_ptr<Surface> rhs)
	//{
	//	mat4& m = rhs.get()->transform;
	//	m = mult(inverse(lhs), m);
	//	return rhs;
	//}

	//std::unique_ptr<Surface> operator*(std::unique_ptr<Surface> lhs, const mat4& rhs)
	//{
	//	mat4& m = lhs.get()->transform;
	//	m = mult(m, inverse(rhs));
	//	return lhs;
	//}
}