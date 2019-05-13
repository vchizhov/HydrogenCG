#pragma once
/*
	@author: Vassillen Chizhov, 2019

	A very simple mathematical vector library
*/

#include <assert.h>
#include <math.h>
#include "math.hpp"
#include "typedef.hpp"

namespace HydrogenCG
{

	struct vec3;
	struct vec4;

	/*
		\brief	A 3-dimensional float-valued vector structure
	*/
	struct vec2
	{
		// this union is not standard conformant, but supported by major compilers
		// allows to alias the same variables through xy, rg, st, e[0],e[1],e[2]
		union
		{
			float e[2];		// e stands for element here
			struct
			{
				float x, y;	//!< Coordinates (components) of the vector
			};
			

			// convenience
			struct
			{
				float r, g;
			};
			struct
			{
				float s, t;
			};
			
		};

		vec2() {}
		vec2(const vec2& arg) : x(arg.x), y(arg.y) {}

		//! Creates a vector with all components equal to 'scalar'
		explicit vec2(float scalar) : x(scalar), y(scalar) {}
		vec2(float x, float y) : x(x), y(y) {}

		vec2& operator=(const vec2& arg)
		{
			x = arg.x;
			y = arg.y;
			return *this;
		}

		float& at(u32 i) { assert(i < 2 && "HydrogenCG::vec2::at():: Array index out of bounds"); return e[i]; }
		const float& at(u32 i) const { assert(i < 2 && "HydrogenCG::vec2::at():: Array index out of bounds"); return e[i]; }
		float& operator()(u32 i) { return e[i]; }
		const float& operator()(u32 i) const { return e[i]; }
		float& operator[](u32 i) { return e[i]; }
		const float& operator[](u32 i) const { return e[i]; }


		// Componentwise / Hadamard operations

		vec2& operator+=(const vec2& arg)
		{
			x += arg.x;
			y += arg.y;
			return *this;
		}

		vec2& operator-=(const vec2& arg)
		{
			x -= arg.x;
			y -= arg.y;
			return *this;
		}

		vec2& operator*=(const vec2& arg)
		{
			x *= arg.x;
			y *= arg.y;
			return *this;
		}

		vec2& operator/=(const vec2& arg)
		{
			x /= arg.x;
			y /= arg.y;
			return *this;
		}

		// Scalar operations
		vec2& operator*=(float scalar)
		{
			x *= scalar;
			y *= scalar;
			return *this;
		}

		vec2& operator/=(float scalar)
		{
			float invScalar = 1.0f / scalar; // note: precision is not as good as dividing each component
			x *= invScalar;
			y *= invScalar;
			return *this;
		}
	};

	inline vec2 operator+(const vec2& arg)
	{
		return arg;
	}

	inline vec2 operator-(const vec2& arg)
	{
		return vec2(-arg.x, -arg.y);
	}

	// Componentwise / Hadamard operations

	inline vec2 operator+(const vec2& lhs, const vec2& rhs)
	{
		return vec2(lhs.x + rhs.x, lhs.y + rhs.y);
	}

	inline vec2 operator-(const vec2& lhs, const vec2& rhs)
	{
		return vec2(lhs.x - rhs.x, lhs.y - rhs.y);
	}

	inline vec2 operator*(const vec2& lhs, const vec2& rhs)
	{
		return vec2(lhs.x * rhs.x, lhs.y * rhs.y);
	}

	inline vec2 operator/(const vec2& lhs, const vec2& rhs)
	{
		return vec2(lhs.x / rhs.x, lhs.y / rhs.y);
	}

	// Scalar operations

	inline vec2 operator*(const vec2& lhs, float rhs)
	{
		return vec2(lhs.x * rhs, lhs.y * rhs);
	}

	inline vec2 operator/(const vec2& lhs, float rhs)
	{
		float invRhs = 1.0f / rhs; // note: precision is not as good as dividing each component
		return vec2(lhs.x * invRhs, lhs.y * invRhs);
	}

	inline vec2 operator*(float lhs, const vec2& rhs)
	{
		return vec2(lhs * rhs.x, lhs * rhs.y);
	}

	inline vec2 operator/(float lhs, const vec2& rhs)
	{
		return vec2(lhs / rhs.x, lhs / rhs.y);
	}

	//! Returns the dot/inner/scalar product of two vectors
	inline float dot(const vec2& lhs, const vec2& rhs)
	{
		return lhs.x * rhs.x + lhs.y * rhs.y;
	}

	//! Returns the squared euclidean norm of a vector
	inline float lengthSquared(const vec2& arg)
	{
		return dot(arg, arg);
	}

	//! Returns the euclidean norm of a vector
	inline float length(const vec2& arg)
	{
		return sqrtf(lengthSquared(arg));
	}

	//! Returns a vector with the same direction as arg, but unit length (undefined for 0)
	inline vec2 normalize(const vec2& arg)
	{
		return arg / length(arg);
	}

	//! Returns the signed parallelogram area of two 2D vectors (||lhs|| * ||rhs|| * sin(lhs,rhs))
	inline float cross(const vec2& lhs, const vec2& rhs)
	{
		return lhs.x * rhs.y - lhs.y * rhs.x;
	}

	inline vec2 clamp(const vec2& arg, float minVal, float maxVal)
	{
		return vec2(clamp(arg.x, minVal, maxVal), clamp(arg.y, minVal, maxVal));
	}

	// convenience componentwise operations
	inline vec2 pow(const vec2& arg, float scalar)
	{
		return vec2(powf(arg.x, scalar), powf(arg.y, scalar));
	}

	inline vec2 max(const vec2& lhs, const vec2& rhs)
	{
		return vec2(max(lhs.x, rhs.x), max(lhs.y, rhs.y));
	}

	inline vec2 max(const vec2& arg, float scalar)
	{
		return vec2(max(arg.x, scalar), max(arg.y, scalar));
	}

	inline float maxComp(const vec2& arg)
	{
		return max(arg.x, arg.y);
	}


	inline vec2 min(const vec2& lhs, const vec2& rhs)
	{
		return vec2(min(lhs.x, rhs.x), min(lhs.y, rhs.y));
	}

	inline vec2 min(const vec2& arg, float scalar)
	{
		return vec2(min(arg.x, scalar), min(arg.y, scalar));
	}

	inline float minComp(const vec2& arg)
	{
		return min(arg.x, arg.y);
	}

	inline vec2 round(const vec2& arg)
	{
		return vec2(roundf(arg.x), roundf(arg.y));
	}

	inline vec2 floor(const vec2& arg)
	{
		return vec2(floorf(arg.x), floorf(arg.y));
	}

	inline vec2 ceil(const vec2& arg)
	{
		return vec2(ceilf(arg.x), ceilf(arg.y));
	}

	/*
		\brief	A 3-dimensional float-valued vector structure
	*/
	struct vec3
	{
		// this union is not standard conformant, but supported by major compilers
		// allows to alias the same variables through xyz, rgb, stp, e[0],e[1],e[2]
		union
		{
			float e[3];		// e stands for element here
			struct
			{
				float x, y, z;	//!< Coordinates (components) of the vector
			};

			// convenience
			struct
			{
				vec2 xy;
				float z;
			};
			struct
			{
				float x;
				vec2 yz;
			};
			struct
			{
				float r, g, b;
			};
			struct
			{
				vec2 rg;
				float b;
			};
			struct
			{
				float r;
				vec2 gb;
			};
			struct
			{
				float s, t, p;
			};
			struct
			{
				vec2 st;
				float p;
			};
			struct
			{
				float s;
				vec2 tp;
			};
		
		};

		vec3() {}
		//! Creates a vector with all components equal to 'scalar'
		explicit vec3(float scalar) : x(scalar), y(scalar), z(scalar) {}
		vec3(float x, float y, float z) : x(x), y(y), z(z) {}
		vec3(const vec3& arg) : x(arg.x), y(arg.y), z(arg.z) {}

		// convenience
		vec3(const vec2& xy, float z) : x(xy.x), y(xy.y), z(z) {}
		vec3(float x, const vec2& yz) : x(x), y(yz.x), z(yz.y) {}

		vec3& operator=(const vec3& arg)
		{
			x = arg.x;
			y = arg.y;
			z = arg.z;
			return *this;
		}

		float& at(u32 i) { assert(i < 3 && "HydrogenCG::vec3::at():: Array index out of bounds"); return e[i]; }
		const float& at(u32 i) const { assert(i < 3 && "HydrogenCG::vec3::at():: Array index out of bounds"); return e[i]; }
		float& operator()(u32 i) { return e[i]; }
		const float& operator()(u32 i) const { return e[i]; }
		float& operator[](u32 i) { return e[i]; }
		const float& operator[](u32 i) const { return e[i]; }


		// Coordinatewise / Hadamard operations

		vec3& operator+=(const vec3& arg)
		{
			x += arg.x;
			y += arg.y;
			z += arg.z;
			return *this;
		}

		vec3& operator-=(const vec3& arg)
		{
			x -= arg.x;
			y -= arg.y;
			z -= arg.z;
			return *this;
		}

		vec3& operator*=(const vec3& arg)
		{
			x *= arg.x;
			y *= arg.y;
			z *= arg.z;
			return *this;
		}

		vec3& operator/=(const vec3& arg)
		{
			x /= arg.x;
			y /= arg.y;
			z /= arg.z;
			return *this;
		}

		// Scalar operations
		vec3& operator*=(float scalar)
		{
			x *= scalar;
			y *= scalar;
			z *= scalar;
			return *this;
		}

		vec3& operator/=(float scalar)
		{
			float invScalar = 1.0f / scalar; // note: precision is not as good as dividing each component
			x *= invScalar;
			y *= invScalar;
			z *= invScalar;
			return *this;
		}
	};

	inline vec3 operator+(const vec3& arg)
	{
		return arg;
	}

	inline vec3 operator-(const vec3& arg)
	{
		return vec3(-arg.x, -arg.y, -arg.z);
	}

	// Coordinatewise / Hadamard operations

	inline vec3 operator+(const vec3& lhs, const vec3& rhs)
	{
		return vec3(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z);
	}

	inline vec3 operator-(const vec3& lhs, const vec3& rhs)
	{
		return vec3(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
	}

	inline vec3 operator*(const vec3& lhs, const vec3& rhs)
	{
		return vec3(lhs.x * rhs.x, lhs.y * rhs.y, lhs.z * rhs.z);
	}

	inline vec3 operator/(const vec3& lhs, const vec3& rhs)
	{
		return vec3(lhs.x / rhs.x, lhs.y / rhs.y, lhs.z / rhs.z);
	}

	// Scalar operations

	inline vec3 operator*(const vec3& lhs, float rhs)
	{
		return vec3(lhs.x * rhs, lhs.y * rhs, lhs.z * rhs);
	}

	inline vec3 operator/(const vec3& lhs, float rhs)
	{
		float invRhs = 1.0f / rhs; // note: precision is not as good as dividing each component
		return vec3(lhs.x * invRhs, lhs.y * invRhs, lhs.z * invRhs);
	}

	inline vec3 operator*(float lhs, const vec3& rhs)
	{
		return vec3(lhs * rhs.x, lhs * rhs.y, lhs * rhs.z);
	}

	inline vec3 operator/(float lhs, const vec3& rhs)
	{
		return vec3(lhs / rhs.x, lhs / rhs.y, lhs / rhs.z);
	}

	//! Returns the dot/inner/scalar product of two vectors
	inline float dot(const vec3& lhs, const vec3& rhs)
	{
		return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
	}

	//! Returns the squared euclidean norm of a vector
	inline float lengthSquared(const vec3& arg)
	{
		return dot(arg, arg);
	}

	//! Returns the euclidean norm of a vector
	inline float length(const vec3& arg)
	{
		return sqrtf(lengthSquared(arg));
	}

	//! Returns a vector with the same direction as arg, but unit length (undefined for 0)
	inline vec3 normalize(const vec3& arg)
	{
		return arg / length(arg);
	}

	//! Returns the cross product of two vectors (a vector orthogonal to both, with length ||lhs|| * ||rhs|| * |sin(lhs,rhs)|)
	inline vec3 cross(const vec3& lhs, const vec3& rhs)
	{
		return vec3(lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z, lhs.x * rhs.y - lhs.y * rhs.x);
	}

	inline vec3 clamp(const vec3& arg, float minVal, float maxVal)
	{
		return vec3(clamp(arg.x, minVal, maxVal), clamp(arg.y, minVal, maxVal), clamp(arg.z, minVal, maxVal));
	}

	// convenience coordinatewise operations
	inline vec3 pow(const vec3& arg, float scalar)
	{
		return vec3(powf(arg.x, scalar), powf(arg.y, scalar), powf(arg.z, scalar));
	}

	inline vec3 max(const vec3& lhs, const vec3& rhs)
	{
		return vec3(max(lhs.x, rhs.x), max(lhs.y, rhs.y), max(lhs.z, rhs.z));
	}

	inline vec3 max(const vec3& arg, float scalar)
	{
		return vec3(max(arg.x, scalar), max(arg.y, scalar), max(arg.z, scalar));
	}

	inline float maxComp(const vec3& arg)
	{
		return max(max(arg.x, arg.y), arg.z);
	}


	inline vec3 min(const vec3& lhs, const vec3& rhs)
	{
		return vec3(min(lhs.x, rhs.x), min(lhs.y, rhs.y), min(lhs.z, rhs.z));
	}

	inline vec3 min(const vec3& arg, float scalar)
	{
		return vec3(min(arg.x, scalar), min(arg.y, scalar), min(arg.z, scalar));
	}

	inline float minComp(const vec3& arg)
	{
		return min(min(arg.x, arg.y), arg.z);
	}

	inline vec3 abs(const vec3& arg)
	{
		return vec3(fabs(arg.x), fabs(arg.y), fabs(arg.z));
	}

	inline vec3 round(const vec3& arg)
	{
		return vec3(roundf(arg.x), roundf(arg.y), roundf(arg.z));
	}

	inline vec3 floor(const vec3& arg)
	{
		return vec3(floorf(arg.x), floorf(arg.y), floorf(arg.z));
	}

	inline vec3 ceil(const vec3& arg)
	{
		return vec3(ceilf(arg.x), ceilf(arg.y), ceilf(arg.z));
	}

	inline vec3 mix(const vec3& a0, const vec3& a1, float t)
	{
		return (1.0f - t) * a0 + t * a1;
	}

	inline vec3 lerp(const vec3& a0, const vec3& a1, float t)
	{
		return mix(a0, a1, t);
	}

	inline vec3 blerp(const vec3& a00, const vec3& a10, 
		const vec3& a01, const vec3& a11, 
		float x, float y)
	{
		return lerp(lerp(a00,a10,x),lerp(a01,a11,x),y);
	}

	inline vec3 tlerp(const vec3& a000, const vec3& a100, 
		const vec3& a010, const vec3& a110,
		const vec3& a001, const vec3& a101,
		const vec3& a011, const vec3& a111,
		float x, float y, float z)
	{
		return lerp(blerp(a000, a100, a010, a110, x, y), blerp(a001, a101, a011, a111, x, y), z);
	}

	/*
		\brief	A 3-dimensional float-valued vector structure
	*/
	struct vec4
	{
		// this union is not standard conformant, but supported by major compilers
		// allows to alias the same variables through xyzw, rgba, stpq, e[0],e[1],e[2], e[4]
		union
		{
			float e[4];		// e stands for element here
			struct
			{
				float x, y, z, w;	//!< Coordinates (components) of the vector
			};

			// convenience
			struct
			{
				vec2 xy;
				vec2 zw;
			};
			struct
			{
				float x;
				vec2 yz;
				float w;
			};
			struct
			{
				vec3 xyz;
				float w;
			};
			struct
			{
				float x;
				vec3 yzw;
			};
			struct
			{
				float r, g, b, a;
			};
			struct
			{
				vec2 rg;
				vec2 ba;
			};
			struct
			{
				float r;
				vec2 gb;
				float a;
			};
			struct
			{
				vec3 rgb;
				float a;
			};
			struct
			{
				float r;
				vec3 fba;
			};
			struct
			{
				float s, t, p, q;
			};
			struct
			{
				vec2 st;
				vec2 pq;
			};
			struct
			{
				float s;
				vec2 tp;
				float q;
			};
			struct
			{
				vec3 stp;
				float q;
			};
			struct
			{
				float s;
				vec3 tpq;
			};
			
		};

		vec4() {}
		//! Creates a vector with all components equal to 'scalar'
		explicit vec4(float scalar) : x(scalar), y(scalar), z(scalar), w(scalar) {}
		vec4(float x, float y, float z, float w) : x(x), y(y), z(z), w(w) {}
		vec4(const vec4& arg) : x(arg.x), y(arg.y), z(arg.z), w(arg.w) {}

		// convenience
		vec4(const vec2& xy, const vec2& zw) : x(xy.x), y(xy.y), z(zw.x), w(zw.y) {}
		vec4(const vec2& xy, float z, float w) : x(xy.x), y(xy.y), z(z), w(w) {}
		vec4(float x, const vec2& yz, float w) : x(x), y(yz.x), z(yz.y), w(w) {}
		vec4(float x, float y, const vec2& zw) : x(x), y(y), z(zw.x), w(zw.y) {}
		vec4(const vec3& xyz, float w) : x(xyz.x), y(xyz.y), z(xyz.z), w(w) {}
		vec4(float x, const vec3& yzw) : x(x), y(yzw.x), z(yzw.y), w(yzw.z) {}
		

		vec4& operator=(const vec4& arg)
		{
			x = arg.x;
			y = arg.y;
			z = arg.z;
			w = arg.w;
			return *this;
		}

		float& at(u32 i) { assert(i < 4 && "HydrogenCG::vec4::at():: Array index out of bounds"); return e[i]; }
		const float& at(u32 i) const { assert(i < 4 && "HydrogenCG::vec4::at():: Array index out of bounds"); return e[i]; }
		float& operator()(u32 i) { return e[i]; }
		const float& operator()(u32 i) const { return e[i]; }
		float& operator[](u32 i) { return e[i]; }
		const float& operator[](u32 i) const { return e[i]; }


		// Componentwise / Hadamard operations

		vec4& operator+=(const vec4& arg)
		{
			x += arg.x;
			y += arg.y;
			z += arg.z;
			w += arg.w;
			return *this;
		}

		vec4& operator-=(const vec4& arg)
		{
			x -= arg.x;
			y -= arg.y;
			z -= arg.z;
			w -= arg.w;
			return *this;
		}

		vec4& operator*=(const vec4& arg)
		{
			x *= arg.x;
			y *= arg.y;
			z *= arg.z;
			w *= arg.w;
			return *this;
		}

		vec4& operator/=(const vec4& arg)
		{
			x /= arg.x;
			y /= arg.y;
			z /= arg.z;
			w /= arg.w;
			return *this;
		}

		// Scalar operations
		vec4& operator*=(float scalar)
		{
			x *= scalar;
			y *= scalar;
			z *= scalar;
			w *= scalar;
			return *this;
		}

		vec4& operator/=(float scalar)
		{
			float invScalar = 1.0f / scalar; // note: precision is not as good as dividing each component
			x *= invScalar;
			y *= invScalar;
			z *= invScalar;
			w *= invScalar;
			return *this;
		}
	};

	inline vec4 operator+(const vec4& arg)
	{
		return arg;
	}

	inline vec4 operator-(const vec4& arg)
	{
		return vec4(-arg.x, -arg.y, -arg.z, -arg.w);
	}

	// Componentwise / Hadamard operations

	inline vec4 operator+(const vec4& lhs, const vec4& rhs)
	{
		return vec4(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z, lhs.w + rhs.w);
	}

	inline vec4 operator-(const vec4& lhs, const vec4& rhs)
	{
		return vec4(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z, lhs.w - rhs.w);
	}

	inline vec4 operator*(const vec4& lhs, const vec4& rhs)
	{
		return vec4(lhs.x * rhs.x, lhs.y * rhs.y, lhs.z * rhs.z, lhs.w * rhs.w);
	}

	inline vec4 operator/(const vec4& lhs, const vec4& rhs)
	{
		return vec4(lhs.x / rhs.x, lhs.y / rhs.y, lhs.z / rhs.z, lhs.w / rhs.w);
	}

	// Scalar operations

	inline vec4 operator*(const vec4& lhs, float rhs)
	{
		return vec4(lhs.x * rhs, lhs.y * rhs, lhs.z * rhs, lhs.w * rhs);
	}

	inline vec4 operator/(const vec4& lhs, float rhs)
	{
		float invRhs = 1.0f / rhs; // note: precision is not as good as dividing each component
		return vec4(lhs.x * invRhs, lhs.y * invRhs, lhs.z * invRhs, lhs.w * invRhs);
	}

	inline vec4 operator*(float lhs, const vec4& rhs)
	{
		return vec4(lhs * rhs.x, lhs * rhs.y, lhs * rhs.z, lhs * rhs.w);
	}

	inline vec4 operator/(float lhs, const vec4& rhs)
	{
		return vec4(lhs / rhs.x, lhs / rhs.y, lhs / rhs.z, lhs / rhs.w);
	}

	//! Returns the dot/inner/scalar product of two vectors
	inline float dot(const vec4& lhs, const vec4& rhs)
	{
		return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z + lhs.w * rhs.w;
	}

	//! Returns the squared euclidean norm of a vector
	inline float lengthSquared(const vec4& arg)
	{
		return dot(arg, arg);
	}

	//! Returns the euclidean norm of a vector
	inline float length(const vec4& arg)
	{
		return sqrtf(lengthSquared(arg));
	}

	//! Returns a vector with the same direction as arg, but unit length (undefined for 0)
	inline vec4 normalize(const vec4& arg)
	{
		return arg / length(arg);
	}

	inline vec4 clamp(const vec4& arg, float minVal, float maxVal)
	{
		return vec4(clamp(arg.x, minVal, maxVal), clamp(arg.y, minVal, maxVal), clamp(arg.z, minVal, maxVal), clamp(arg.w, minVal, maxVal));
	}

	// convenience componentwise operations
	inline vec4 pow(const vec4& arg, float scalar)
	{
		return vec4(powf(arg.x, scalar), powf(arg.y, scalar), powf(arg.z, scalar), powf(arg.w, scalar));
	}

	inline vec4 max(const vec4& lhs, const vec4& rhs)
	{
		return vec4(max(lhs.x, rhs.x), max(lhs.y, rhs.y), max(lhs.z, rhs.z), max(lhs.w, rhs.w));
	}

	inline vec4 max(const vec4& arg, float scalar)
	{
		return vec4(max(arg.x, scalar), max(arg.y, scalar), max(arg.z, scalar), max(arg.w, scalar));
	}

	inline float maxComp(const vec4& arg)
	{
		return max(max(max(arg.x, arg.y), arg.z),arg.w);
	}


	inline vec4 min(const vec4& lhs, const vec4& rhs)
	{
		return vec4(min(lhs.x, rhs.x), min(lhs.y, rhs.y), min(lhs.z, rhs.z), min(lhs.w,rhs.w));
	}

	inline vec4 min(const vec4& arg, float scalar)
	{
		return vec4(min(arg.x, scalar), min(arg.y, scalar), min(arg.z, scalar), min(arg.w, scalar));
	}

	inline float minComp(const vec4& arg)
	{
		return min(min(min(arg.x, arg.y), arg.z), arg.w);
	}

	inline vec4 round(const vec4& arg)
	{
		return vec4(roundf(arg.x), roundf(arg.y), roundf(arg.z), roundf(arg.w));
	}

	inline vec4 floor(const vec4& arg)
	{
		return vec4(floorf(arg.x), floorf(arg.y), floorf(arg.z), floorf(arg.w));
	}

	inline vec4 ceil(const vec4& arg)
	{
		return vec4(ceilf(arg.x), ceilf(arg.y), ceilf(arg.z), ceilf(arg.w));
	}

}