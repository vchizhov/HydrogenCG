#pragma once
#include "..\core\surface_instance.hpp"

namespace HydrogenCG
{


	/*!
		\brief	A mathematical description of a quadric
	*/
	struct Quadric : public Surface
	{
		vec3 o;		//!< the origin/center of the quadric
		mat4 m;		//!< the matrix defining the coefficients of the quadric, needs to be symmetric for sane results

		Quadric() {}
		explicit Quadric(const vec3& origin, const mat4& m) : o(origin), m(m) {}

		/*!
			\brief		Quadric intersection
		*/
		Intersection intersect(const Ray& ray, float minT = 0.0f, float maxT = std::numeric_limits<float>::infinity()) const final
		{
			/*
				Quadric intersection derivation:
				a_11 * x^2 + a_22 * y^2 + a_33 * z^2 +
				(a_12 + a_21) * x*y + (a_13 + a_31) * x*z + (a_23+a_32) * y*z +
				(a_14 + a_41) * x + (a_24 + a_42) * y + (a_34 + a_43) * z + a_44 = 0

					Rewrite as a vector-matrix-vector equation and plug in the ray:
													<=>

					 <(x,y,z,1), A * (x,y,z,1)> = 0 <=> <(r(t)-c,1),A * (r(t)-c,1)> = 0

									 Multiply both side by 1 = (-1)*(-1)
													<=>

								(-1)*(-1) * <(r(t)-c,1),A * (r(t)-c,1)> = 0

						Flip the sign using the commutativity of scalar-matrix product:
													<=>

				 <(c-r(t),-1),A * (c-r(t),-1)> = 0  <=>  <(c-o-t*d,-1),A * (c-o-t*d,-1)> = 0

			 Use the additivity additivity of vectors and vector-scalar product properties to separate t:
													<=>

								<[(c-o,-1)-t * (d,0)],A * [(c-o,-1)-t * (d,0)]> = 0	

					Use the distributivity of the matrix (inner) product to separate the terms with t:
													<=>

					<[(c-o,-1)-t * (d,0)],A * (c-o,-1)> + <[(c-o,-1)-t * (d,0)],A * -t * (d,0)> = 0
												
					Use the distributivity of the matrix (inner) product to separate the terms with t:
													<=>

			<(c-o,-1),A * (c-o,-1)> - <t * (d,0),A * (c-o,-1)> - <(c-o,-1),A * t * (d,0)> + <-t * (d,0),A * -t * (d,0)> = 0

			Use the property that the adjoint of a real matrix A is A^T (A transposed), and scalar distributivity to take out t:
													<=>

			<(c-o,-1),A * (c-o,-1)> - <(d,0),A * (c-o,-1)> * t - <A^T * (c-o,-1), (d,0)> * t + <(d,0),A * (d,0)> * t^2 = 0

				Use the symmetry of the inner product <u,v> = u^T * v = ((u^T * v)^T)^T = (v^T * u)^T = (<v,u>)^T = <v,u>:
													<=>

			<(c-o,-1),A * (c-o,-1)> - <(d,0),A * (c-o,-1)> * t - <(d,0), A^T * (c-o,-1)> * t + <(d,0),A * (d,0)> * t^2 = 0

			    Require that A is square symmetric (A^T=A), this is not an actual restriction in our case since the terms (a_12+a_21) 
			might as well be written as (a_12+a_21) = (a_12+a_12) = 2 * a_12 without changing the expressvity of the initial formulation:
													<=>

				<(c-o,-1),A * (c-o,-1)> - <(d,0),A * (c-o,-1)> * t - <(d,0), A * (c-o,-1)> * t + <(d,0),A * (d,0)> * t^2 = 0

													<=>

					<(d,0),A * (d,0)> * t^2 - 2*<(d,0),A * (c-o,-1)> * t + <(c-o,-1),A * (c-o,-1)> = 0
												A * t^2 - 2*B * t + C = 0
					A = <(d,0),A * (d,0)>, B = <(d,0),A * (c-o,-1)>, C = <(c-o,-1),A * (c-o,-1)>, D = B^2 - A*C

													<=>
									t_1 = (B-sqrt(D))/A, t_2 = (B+sqrt(D))/A

			*/

			// (c-o,-1)
			vec4 ro = vec4(o - ray.o,-1.0f);
			// (d,0)
			vec4 d = vec4(ray.d, 0);

			// A * (c-o,-1)
			vec4 mro = m * ro;

			// A = <(d,0), A * (d,0)>
			float a = dot(d, m*d);
			// B = <(d,0),A * (c-o,-1)>
			float b = dot(d, mro);
			// C = <(c-o,-1),A * (c-o,-1)>
			float c = dot(ro, mro);
			// D = B^2 - A*C
			float discriminant = b * b - a * c;

			// If the discriminant is 0 or negative -> no (real) intersection
			if (discriminant <= 0.0f) return Intersection();

			float sqrtD = sqrt(discriminant);

			// handle negative values of the A term - specifically to get the closer/further intersection
			float invA;
			if (a < 0.0f)
			{
				b = -b;
				invA = -1.0f / a;
			}
			else
			{
				invA = 1.0f / a;
			}
			
			// t_1 = (B-sqrt(D))/A
			float tN = (b - sqrtD) * invA; // closer intersection (where the ray enters the sphere)

			// if it is within the defined ray segment by (minT,maxT) return it as the closest intersection
			if (tN > minT && tN < maxT)
			{
				vec3 normal;
				vec2 uv;
				normalAndUV(ray(tN), normal, uv);
				return Intersection(tN, normal, uv);
			}
			// if we didn't get the closer intersection, it may be that we are "inside" the quadric,
			// so we still need to intersect the further intersection

			// t_2 = (B+sqrt(D))/A
			float tF = (b + sqrtD)*invA; // further intersection (where the ray exits the sphere)

			// if it is within the defined ray segment by (minT,maxT) return it as the closest intersection
			if (tF > minT && tF < maxT)
			{
				vec3 normal;
				vec2 uv;
				normalAndUV(ray(tF), normal, uv);
				return Intersection(tF, normal, uv);
			}

			// otherwise return no intersection
			return Intersection();
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
			// We know that the normal of an implicit surface f(p) = 0, is given as:
			// normalize(d/dp(f)(p)), where df/dp is to be understood as the gradient.
			// In our case: f(p) = (p-o,1)^T * A * (p-o,1) = x^T * A * x
			// We have: d(x^T * A * x)/dx = (A^T+A) * x, where x is a vector
			// since A=A^T -> n = normalize(2 * A * x)
			normal = normalize((m * vec4((p - o), 1)).xyz);
			// I did not derive a general expression for the uv coordinates -> just return the x,y intersection coordinates
			uv = vec2(p.x,p.y);
		}

		/*!
			Returns some "approximation of the distance" between a point p and a quadric - 
			not an accurate representation of the distance (works for inside/outside checks though)
		*/
		float distance(const vec3& p) const final
		{
			return dot(vec4(p, 1), m*vec4(p, 1));
		}

		// a_11 * x^2 + a_22 * y^2 + a_33 * z^2 - 1 = 0 -> a_11=a_22=a_33=1, a_44=-r^2; set  r=1
		static const std::shared_ptr<const Quadric> unitSphere;
		static std::shared_ptr<SurfaceInstanceBSDF> sphere()
		{
			return std::make_shared<SurfaceInstanceBSDF>(Quadric::unitSphere);
		}
		static std::shared_ptr<SurfaceInstanceBSDF> sphere(const vec3& origin, float radius = 1.0f, const vec3& col = vec3(1))
		{
			mat4 sphereMat = mat4::identity;
			sphereMat[3][3] = -radius * radius;
			return std::make_unique<SurfaceInstanceBSDF>(std::make_unique<Quadric>(origin, sphereMat), bsdfDiffuse(col));
		}

		// a_11 * x^2 + a_33 * z^2 - 1 = 0 -> a_11 = 1/(r_1)^2, a_33 = 1/(r_2)^2; set r_1=r_2=1
		static const std::shared_ptr<const Quadric> unitYCylinder;

		// a_11 * x^2 - a_22 * y^2 + a_33 * z^2= 0 -> 1/(r_1)^2, a_33 = 1/(r_2)^2, a_22 = yScaling^2; set r_1=r_2=yScaling=1
		static const std::shared_ptr<const Quadric> unitYCone;

		// a_11 * x^2 - (a_24 + a_42) * y + a_33 * z^2 = 0 -> a_11 = 1/(r_1)^2, a_33 = 1/(r_2)^2, a_24=a_42=0.5*yScaling; set r_1=r_2=1 = yScaling
		static const std::shared_ptr<const Quadric> unitYParaboloid;

		// a_11 * x^2 - a_22 * y^2 + a_33 * z^2 - 1 = 0 -> 1/(r_1)^2, a_33 = 1/(r_2)^2, a_22 = 1/yScaling^2; set r_1=r_2=yScaling=1
		static const std::shared_ptr<const Quadric> unitYHyperboloid; // single sheet

		// a_11 * x^2 - a_22 * y^2 + a_33 * z^2 + 1 = 0 -> 1/(r_1)^2, a_33 = 1/(r_2)^2, a_22 = 1/yScaling^2; set r_1=r_2=yScaling=1
		//static const std::shared_ptr<const Surface> unitYHyperboloid2; // two sheets

		// a_11 * x^2 - (a_24+a_42) * y - a_33 * z^2 = 0 -> 1/(r_1)^2, a_33 = 1/(r_2)^2, a_24=a_42 = 0.5*yScaling; set r_1=r_2=yScaling=1
		static const std::shared_ptr<const Quadric> unitYHyperbolicParaboloid;

		// a_11 * x^2 - a_33 * z^2 - 1 = 0 -> 1/(r_1)^2, a_33 = 1/(r_2)^2; set r_1=r_2=1
		static const std::shared_ptr<const Quadric> unitYHyperbolicCylinder;

		// a_11 * x^2 + (a_34+a_43)*z = 0 -> 1/(r_1)^2, a_32=a_43 = zScaling; set r_1=1, a_32=a_43 = 1
		static const std::shared_ptr<const Quadric> unitYParabolicCylinder;
	};
	const std::shared_ptr<const Quadric> Quadric::unitSphere = std::make_unique<const Quadric>(vec3(0.0f), mat4(vec4(1,0,0,0),
		vec4(0,1,0,0),vec4(0,0,1,0), vec4(0,0,0,-1)));
	

	

	const std::shared_ptr<const Quadric> Quadric::unitYCylinder = std::make_shared<const Quadric>(vec3(0.0f), mat4(vec4(1, 0, 0, 0),
		vec4(0, 0, 0, 0), vec4(0, 0, 1, 0), vec4(0, 0, 0, -1)));
	std::shared_ptr<SurfaceInstanceBSDF> quadricCylinder()
	{
		return std::make_shared<SurfaceInstanceBSDF>(Quadric::unitYCylinder);
	}

	const std::shared_ptr<const Quadric> Quadric::unitYCone = std::make_shared<const Quadric>(vec3(0.0f), mat4(vec4(1, 0, 0, 0),
		vec4(0, -1, 0, 0), vec4(0, 0, 1, 0), vec4(0, 0, 0, 0)));
	std::shared_ptr<SurfaceInstanceBSDF> quadricCone()
	{
		return std::make_shared<SurfaceInstanceBSDF>(Quadric::unitYCone);
	}

	const std::shared_ptr<const Quadric> Quadric::unitYParaboloid = std::make_shared<const Quadric>(vec3(0.0f), mat4(vec4(1, 0, 0, 0),
		vec4(0, 0, 0, -0.5f), vec4(0, 0, 1, 0), vec4(0, -0.5f, 0, 0)));
	std::shared_ptr<SurfaceInstanceBSDF> quadricParaboloid()
	{
		return std::make_shared<SurfaceInstanceBSDF>(Quadric::unitYParaboloid);
	}

	const std::shared_ptr<const Quadric> Quadric::unitYHyperboloid = std::make_shared<const Quadric>(vec3(0.0f), mat4(vec4(1, 0, 0, 0),
		vec4(0, -1, 0, 0), vec4(0, 0, 1, 0), vec4(0, 0, 0, -1)));
	std::shared_ptr<SurfaceInstanceBSDF> quadricHyperboloid()
	{
		return std::make_shared<SurfaceInstanceBSDF>(Quadric::unitYHyperboloid);
	}

	const std::shared_ptr<const Quadric> Quadric::unitYHyperbolicParaboloid = std::make_shared<const Quadric>(vec3(0.0f), mat4(vec4(1, 0, 0, 0),
		vec4(0, 0, 0, -0.5f), vec4(0, 0, -1, 0), vec4(0, -0.5f, 0, 0)));
	std::shared_ptr<SurfaceInstanceBSDF> quadricHyperbolicParaboloid()
	{
		return std::make_shared<SurfaceInstanceBSDF>(Quadric::unitYHyperbolicParaboloid);
	}

	const std::shared_ptr<const Quadric> Quadric::unitYHyperbolicCylinder = std::make_shared<const Quadric>(vec3(0.0f), mat4(vec4(1, 0, 0, 0),
		vec4(0, 0, 0, 0), vec4(0, 0, -1, 0), vec4(0, 0, 0, -1)));
	std::shared_ptr<SurfaceInstanceBSDF> quadricHyperbolicCylinder()
	{
		return std::make_shared<SurfaceInstanceBSDF>(Quadric::unitYHyperbolicCylinder);
	}


	//self intersection artifacts with shadow rays
	const std::shared_ptr<const Quadric> Quadric::unitYParabolicCylinder = std::make_shared<const Quadric>(vec3(0.0f), mat4(vec4(1, 0, 0, 0),
		vec4(0, 0, 0, 0), vec4(0, 0, 0, 1), vec4(0, 0, 1, 0)));
	std::shared_ptr<SurfaceInstanceBSDF> quadricParabolicCylinder()
	{
		return std::make_shared<SurfaceInstanceBSDF>(Quadric::unitYParabolicCylinder);
	}
}