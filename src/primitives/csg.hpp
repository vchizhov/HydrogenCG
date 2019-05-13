#pragma once

#include <vector>
#include "..\core\surface.hpp"
#include "..\core\surface_instance.hpp"

namespace HydrogenCG
{
	/*

		We will use * for set intersection: A intersection B = A * B
		We will use unary - for set complement: complement(A) = -A
		For union we will use +: A union B = A + B
		For set difference we will use binary -: difference(A,B) = A - B
		For symmetric difference we will use ^: symmetricDifference(A,B) = A^B

		Complement and intersection are enough for functional completeness.

		For instance the following set operations can be defined through complement and intersection only:
		A-B = A * -B
		A+B = -(-A)+(-(-B)) = -((-A)*(-B)), where the second equality holds due to De Morgan's law
		A^B = (A-B) + (B-A) = -(-(A*-B)*(-(B*-A)))
		
	*/


	/*
		\brief Takes an implicitly defined volume through a surface and a distance function as an argument and returns its complement
	*/
	class SurfaceInstanceCSGComplement : public SurfaceInstance
	{
	public:
		std::shared_ptr<SurfaceInstance> a;

		SurfaceInstanceCSGComplement(std::shared_ptr<SurfaceInstance> a) : a(a) {}

		InstanceIntersection intersectOriginal(const Ray& ray, float minT = 0.0f, float maxT = std::numeric_limits<float>::infinity()) const final
		{
			InstanceIntersection isect = a->intersect(ray, minT, maxT);
			// flips the normal on intersection, the complement "flips" the surface
			// the inside of a closed surface becomes the outside
			if (isect)
				return InstanceIntersection(Intersection(isect.t,-isect.n,isect.uv),isect.bsdf);
			else
				return InstanceIntersection();
		}

		std::shared_ptr<SurfaceInstance> doClone() const final
		{
			auto inst = std::make_shared<SurfaceInstanceCSGComplement>(a->Clone());
			inst->negativePos = negativePos;
			inst->transposedRot = transposedRot;
			inst->invScale = invScale;

			return inst;
		}

		std::shared_ptr<SurfaceInstanceCSGComplement> Clone() const
		{
			return std::static_pointer_cast<SurfaceInstanceCSGComplement>(doClone());
		}

		// note - for CSG a set membership function can be used that isn't the exact distance
		float distanceOriginal(const vec3& p) const final
		{
			// since the inside becomes the outside we simply return the negative distance
			return -a->distance(p);
		}
	};



	/*
		\brief Takes two implicitly defined volumes as an argument and returns the volume resulting from their intersection
	*/
	class SurfaceInstanceCSGIntersection : public SurfaceInstance
	{
	public:
		std::shared_ptr<SurfaceInstance> a;
		std::shared_ptr<SurfaceInstance> b;

		SurfaceInstanceCSGIntersection(std::shared_ptr<SurfaceInstance> a, std::shared_ptr<SurfaceInstance> b) : a(a), b(b) {}
		//CSGIntersection(const CSGIntersection& arg) : SurfaceInstance(arg), a(std::make_unique<SurfaceInstance>(*a.get())), b(std::make_unique<SurfaceInstance>(*b.get())) {}

		InstanceIntersection intersectOriginal(const Ray& ray, float minT = 0.0f, float maxT = std::numeric_limits<float>::infinity()) const final
		{
			// check whether the ray origin is inside any of the volumes
			bool insideA = a->distance(ray.o) <= 0.0f;
			bool insideB = b->distance(ray.o) <= 0.0f;

			InstanceIntersection iA = a->intersect(ray, minT, maxT);

			// if the ray origins is not inside the first volume and does not intersect its boundary -> the intersection set is empty
			// since it never enters the first volume
			if (!insideA && !iA) return InstanceIntersection();

			InstanceIntersection iB = b->intersect(ray, minT, maxT);

			// if the ray origins is not inside the second volume and does not intersect its boundary -> the intersection is empty
			// since it never enters the second volume
			if (!insideB && !iB) return InstanceIntersection();

			// even if the ray origin is inside both volumes, if it does not intersect their boundaries -> there's nothing to see
			// so return no intersection
			if (!iA && !iB) return InstanceIntersection();

			// loop until you find a hit, or miss everything
			// in general the normal offset should be used and not iA.t+EPSILON similarly to what we do in the transparency integrator, 
			// but we keep it that way for simplicity
			int iter = 0;
			const int maxIter = 100;
			while (iter<maxIter)
			{
				++iter;
				if (!insideA && !insideB) // the current point on the ray is outside of both volumes
				{
					if (!iA || !iB) // misses at least one of the volumes, and is outside both volumes -> no intersection
					{
						return InstanceIntersection();
					}
					else // intersects both volumes
					{
						// advance the closer intersection point, we denote the closer volume C and the farther volume F
						InstanceIntersection& iCloser = iA <= iB ? iA : iB;
						InstanceIntersection& iFarther = iA <= iB ? iB : iA;
						const SurfaceInstance* closer = iA <= iB ? a.get() : b.get();
						bool& insideCloser = iA <= iB ? insideA : insideB;
						bool& insideFarther = iA <= iB ? insideB : insideA;

						// the distance at which it's supposed to exit volume C
						iCloser = closer->intersect(ray, iCloser.t + EPSILON, maxT);
						// if there's no intersection, it does not exit C ever (C is an unbounded volume)
						// we already know there's an intersection where it enters F, and it is in C at the same time-> return that intersection
						// If there is an intersection where it exits volume C, but that intersection is further than the intersection at which we
						// enter F, then F is the intersection we're looking for since it's in both volumes (thus inside their set intersection)
						if (iFarther<=iCloser)
						{
							return iFarther;
						}

						// here we know that the ray has exited C for sure, and that it has not entered F yet, so we are once again
						// in the case where it's outside both volumes, we need to advance the closer intersection again to have some 
						// new data for the next iteration of the loop
						insideFarther = false;
						insideCloser = false;
						iCloser = closer->intersect(ray, iCloser.t + EPSILON, maxT);
					}
				}
				else if (insideA && insideB) // is inside both volumes
				{
					if (!iA && !iB) // misses both volumes -> even though it's inside both, there's nothing to see
					{
						return InstanceIntersection();
					}
					else // intersects at least one of the volumes -> return the closer intersection
					{
						return iA <= iB ? iA : iB;
					}
				}
				// The two cases below are symmetric and may be merged in one case, but are left as it is for clarity
				else if (insideA && !insideB) // is inside volume A, but not inside volume B
				{
					if (!iB) // if we do not intersect B, we never enter it -> no intersection
					{
						return InstanceIntersection();
					}
					else // there are two cases: either we intersect B before leaving A -> return it, or we leave A before B, and we need to perform more checks
					{
						if (iB <= iA) // this also handles the case where there is no intersection for A
							return iB;
						else // here we are sure we have intersected A, and that its intersection is closer
						{
							// advance A since we need to be inside A and B at the same time
							iA = a->intersect(ray, iA.t+EPSILON, maxT);
							// we were inside A, exited it at the previous value of iA.t, and now we possibly enter it again if there is an intersection
							// if we do not enter it again -> the intersection is empty
							if (!iA) return InstanceIntersection();
							// we have exited both, and intersected both -> it goes in the case at the top
							insideA = false;
							insideB = false;
						}
					}
				}
				else //if (!insideA && insideB) // a symmetric case to the above
				{
					if (!iA) // if we do not intersect B, we never enter it -> no intersection
					{
						return InstanceIntersection();
					}
					else // there are two cases: either we intersect A before leaving B -> return it, or we leave B before A, and we need to perform more checks
					{
						if (iA <= iB) // this also handles the case where there is no intersection for A
							return iA;
						else // here we are sure we have intersected B, and that its intersection is closer
						{
							// advance B since we need to be inside B and A at the same time
							iB = b->intersect(ray, iB.t + EPSILON, maxT);
							// we were inside B, exited it at the previous value of iB.t, and now we possibly enter it again (if ther was an intersection
							// if we do not enter it again -> the intersection is empty
							if (!insideB) return InstanceIntersection();
							// check whether we are still inside B at the new position
							// we have exited both, and intersected both -> it goes in the case at the top
							insideA = false;
							insideB = false;

						}
					}
				}
			}
			return InstanceIntersection();
		}

		std::shared_ptr<SurfaceInstance> doClone() const final
		{
			auto inst = std::make_shared<SurfaceInstanceCSGIntersection>(a->Clone(), b->Clone());
			inst->negativePos = negativePos;
			inst->transposedRot = transposedRot;
			inst->invScale = invScale;

			return inst;
		}

		std::shared_ptr<SurfaceInstanceCSGIntersection> Clone() const
		{
			return std::static_pointer_cast<SurfaceInstanceCSGIntersection>(doClone());
		}

		float distanceOriginal(const vec3& p) const final
		{
			return max(a->distance(p), b->distance(p));
		}
	};

	// set intersection operator overload for lvalues
	template<class A, class B>
	std::shared_ptr<SurfaceInstanceCSGIntersection> operator*(std::shared_ptr<A>& lhs, std::shared_ptr<B>& rhs)
	{
		return std::make_shared<SurfaceInstanceCSGIntersection>(lhs->Clone(), rhs->Clone());
	}

	// set intersection operator overload for rvalues
	template<class A, class B>
	std::shared_ptr<SurfaceInstanceCSGIntersection> operator*(std::shared_ptr<A>&& lhs, std::shared_ptr<B>&& rhs)
	{
		return std::make_shared<SurfaceInstanceCSGIntersection>(lhs, rhs);
	}

	// set complement operator overload for lvalues
	template<class A>
	std::shared_ptr<SurfaceInstanceCSGComplement> operator-(std::shared_ptr<A>& arg)
	{
		return std::make_shared<SurfaceInstanceCSGComplement>(arg->Clone());
	}

	// set complement operator overload for rvalues
	template<class A>
	std::shared_ptr<SurfaceInstanceCSGComplement> operator-(std::shared_ptr<A>&& arg)
	{
		return std::make_shared<SurfaceInstanceCSGComplement>(arg);
	}

	// set union operator overload for lvalues
	template<class A, class B>
	std::shared_ptr<SurfaceInstanceCSGIntersection> operator+(std::shared_ptr<A>& lhs, std::shared_ptr<B>& rhs)
	{
		// using De Morgan's law:
		// -(-A*-B) = A+B
		return -((-lhs->Clone()) * (-rhs->Clone()));
	}

	// set union operator overload for rvalues
	template<class A, class B>
	std::shared_ptr<SurfaceInstanceCSGIntersection> operator+(std::shared_ptr<A>&& lhs, std::shared_ptr<B>&& rhs)
	{
		// using De Morgan's law:
		// -(-A*-B) = A+B
		return -((-lhs) * (-rhs));
	}

	// set difference operator overload for lvalues
	template<class A, class B>
	std::shared_ptr<SurfaceInstanceCSGIntersection> operator-(std::shared_ptr<A>& lhs, std::shared_ptr<B>& rhs)
	{
		// A*-B = A-B
		return lhs->Clone() * (-rhs->Clone());
	}

	// set difference operator overload for rvalues
	template<class A, class B>
	std::shared_ptr<SurfaceInstanceCSGIntersection> operator-(std::shared_ptr<A>&& lhs, std::shared_ptr<B>&& rhs)
	{
		// A*-B = A-B
		return lhs * (-rhs);
	}

	// symmetric difference operator overload for rvalues
	template<class A, class B>
	std::shared_ptr<SurfaceInstance> operator^(std::shared_ptr<A>& lhs, std::shared_ptr<B>& rhs)
	{
		// (B-A) + (A-B) = B*-A + A*-B
		auto lhsCloned = lhs->Clone();
		auto rhsCloned = rhs->Clone();
		return (rhsCloned - lhsCloned) + (lhsCloned - rhsCloned);
	}

	// symmetric difference operator overload for lvalues
	template<class A, class B>
	std::shared_ptr<SurfaceInstance> operator^(std::shared_ptr<A>&& lhs, std::shared_ptr<B>&& rhs)
	{
		// (B-A) + (A-B) = B*-A + A*-B
		return (rhs - lhs) + (lhs - rhs);
	}

}