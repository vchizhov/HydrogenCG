#pragma once
#include <assert.h>
#include <math.h>
#include "math.h"

namespace HydrogenCG
{
	struct vec3
	{
		// this union is not standard conformant, but supported by all major compilers
		// allows to use aliases for the variables
		union
		{
			struct
			{
				float x, y, z;
			};
			struct
			{
				float r, g, b;
			};
			struct
			{
				float s, t, p;
			};
			float e[3];
		};

		// Imo we should not have default constructors that do nothing.
		// 0 initialization should be the default, and we can use tags to
		// create constructors that leave data uninitialized.
		// This helps people from all skill levels to prevent bugs from
		// happening caused by uninitialized data.
		vec3() {}
		explicit vec3(float scalar) : x(scalar), y(scalar), z(scalar) {}
		vec3(float x, float y, float z) : x(x), y(y), z(z) {}
		vec3(const vec3& arg) : x(arg.x), y(arg.y), z(arg.z) {}

		vec3& operator=(const vec3& arg)
		{
			x = arg.x;
			y = arg.y;
			z = arg.z;
			return *this;
		}

		float& at(int i) { return e[i]; }
		const float& at(int i) const { return e[i]; }
		float& operator()(int i) { assert(i >= 0 && i < 3 && "HydrogenCG::vec3::operator():: Array index out of bounds"); return at(i); }
		const float& operator()(int i) const { assert(i >= 0 && i < 3 && "HydrogenCG::vec3::operator():: Array index out of bounds"); return at(i); }
		float& operator[](int i) { assert(i >= 0 && i < 3 && "HydrogenCG::vec3::operator():: Array index out of bounds"); return at(i); }
		const float& operator[](int i) const { assert(i >= 0 && i < 3 && "HydrogenCG::vec3::operator():: Array index out of bounds"); return at(i); }


		// Componentwise / Hadamard operations

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

	// Componentwise / Hadamard operations

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

	inline float dot(const vec3& lhs, const vec3& rhs)
	{
		return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
	}

	inline float lengthSquared(const vec3& arg)
	{
		return dot(arg, arg);
	}

	inline float length(const vec3& arg)
	{
		return sqrtf(lengthSquared(arg));
	}

	inline vec3 normalize(const vec3& arg)
	{
		return arg / length(arg);
	}

	inline vec3 cross(const vec3& lhs, const vec3& rhs)
	{
		return vec3(lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z, lhs.x * rhs.y - lhs.y * rhs.x);
	}

	vec3 clamp(const vec3& arg, float minVal, float maxVal)
	{
		return vec3(clamp(arg.x, minVal, maxVal), clamp(arg.y, minVal, maxVal), clamp(arg.z, minVal, maxVal));
	}

}
