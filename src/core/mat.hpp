#pragma once
#include "vec.hpp"


namespace HydrogenCG
{
	struct mat3
	{
		union
		{
			struct
			{
				vec3 row0, row1, row2;
			};
			vec3 r[3];
			float e[3][3];
		};

		mat3() {}
		explicit mat3(float scalar) : row0(scalar), row1(scalar), row2(scalar) {}
		mat3(float e00, float e01, float e02,
			float e10, float e11, float e12,
			float e20, float e21, float e22)
			: row0(e00,e01,e02), 
			row1(e10, e11, e12), 
			row2(e20, e21, e22)
		{}
		mat3(const vec3& row0, const vec3& row1, const vec3& row2)
			: row0(row0), row1(row1), row2(row2) {}
		mat3(const mat3& arg) : row0(arg.row0), row1(arg.row1), row2(arg.row2) {}

		mat3& operator=(const mat3& arg)
		{
			row0 = arg.row0;
			row1 = arg.row1;
			row2 = arg.row2;
			return *this;
		}

		mat3& fromColumns(const vec3& col0, const vec3& col1, const vec3& col2)
		{
			row0 = vec3(col0[0], col1[0], col2[0]);
			row1 = vec3(col0[1], col1[1], col2[1]);
			row2 = vec3(col0[2], col1[2], col2[2]);
			return *this;
		}

		vec3& at(u32 rowIdx) { assert(rowIdx < 3 && "HydrogenCG::mat3::at():: Row index out of bounds"); return r[rowIdx]; }
		const vec3& at(u32 rowIdx) const { assert(rowIdx < 3 && "HydrogenCG::mat3::at():: Row index out of bounds"); return r[rowIdx]; }
		vec3& operator[](u32 rowIdx) { return r[rowIdx]; }
		const vec3& operator[](u32 rowIdx) const { return r[rowIdx]; }
		vec3& operator()(u32 rowIdx) { return r[rowIdx]; }
		const vec3& operator()(u32 rowIdx) const { return r[rowIdx]; }

		float& at(u32 rowIdx, u32 colIdx) 
		{ 
			assert(rowIdx < 3 && "HydrogenCG::mat3::at():: Row index out of bounds."); 
			assert(colIdx < 3 && "HydrogenCG::mat3::at():: Column index out of bounds.");
			return e[rowIdx][colIdx]; 
		}
		const float& at(u32 rowIdx, u32 colIdx) const
		{
			assert(rowIdx < 3 && "HydrogenCG::mat3::at():: Row index out of bounds.");
			assert(colIdx < 3 && "HydrogenCG::mat3::at():: Column index out of bounds.");
			return e[rowIdx][colIdx];
		}
		float& operator()(u32 rowIdx, u32 colIdx) { return r[rowIdx][colIdx]; }
		const float& operator()(u32 rowIdx, u32 colIdx) const { return r[rowIdx][colIdx]; }


		const mat3& operator+=(const mat3& rhs)
		{
			row0 += rhs.row0;
			row1 += rhs.row1;
			row2 += rhs.row2;
			return *this;
		}

		const mat3& operator-=(const mat3& rhs)
		{
			row0 -= rhs.row0;
			row1 -= rhs.row1;
			row2 -= rhs.row2;
			return *this;
		}

		// Hadamard product
		const mat3& operator*=(const mat3& rhs)
		{
			row0 *= rhs.row0;
			row1 *= rhs.row1;
			row2 *= rhs.row2;
			return *this;
		}

		// Hadamard division
		const mat3& operator/=(const mat3& rhs)
		{
			row0 /= rhs.row0;
			row1 /= rhs.row1;
			row2 /= rhs.row2;
			return *this;
		}

		// Scalar operations
		const mat3& operator*=(float scalar)
		{
			row0 *= scalar;
			row1 *= scalar;
			row2 *= scalar;
			return *this;
		}

		const mat3& operator/=(float scalar)
		{
			float invScalar = 1 / scalar;
			row0 *= invScalar;
			row1 *= invScalar;
			row2 *= invScalar;
			return *this;
		}

		static const mat3 identity;
	};

	const mat3 mat3::identity = mat3(vec3(1, 0, 0), vec3(0, 1, 0), vec3(0, 0, 1));

	mat3 fromColumns(const vec3& col0, const vec3& col1, const vec3& col2)
	{
		return mat3(vec3(col0[0], col1[0], col2[0]),
			vec3(col0[1], col1[1], col2[1]),
			vec3(col0[2], col1[2], col2[2]));
	}

	inline mat3 operator+(const mat3& arg)
	{
		return arg;
	}

	inline mat3 operator-(const mat3& arg)
	{
		return mat3(-arg[0], -arg[1], -arg[2]);
	}

	// Componentwise / Hadamard operations

	inline mat3 operator+(const mat3& lhs, const mat3& rhs)
	{
		return mat3(lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2]);
	}

	inline mat3 operator-(const mat3& lhs, const mat3& rhs)
	{
		return mat3(lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2]);
	}

	inline mat3 operator*(const mat3& lhs, const mat3& rhs)
	{
		return mat3(lhs[0] * rhs[0], lhs[1] * rhs[1], lhs[2] * rhs[2]);
	}

	inline mat3 operator/(const mat3& lhs, const mat3& rhs)
	{
		return mat3(lhs[0] / rhs[0], lhs[1] / rhs[1], lhs[2] / rhs[2]);
	}

	// Scalar operations

	inline mat3 operator*(const mat3& lhs, float rhs)
	{
		return mat3(lhs[0] * rhs, lhs[1] * rhs, lhs[2] * rhs);
	}

	inline mat3 operator/(const mat3& lhs, float rhs)
	{
		float invRhs = 1.0f / rhs; // note: precision is not as good as dividing each component
		return mat3(lhs[0] * invRhs, lhs[1] * invRhs, lhs[2] * invRhs);
	}

	inline mat3 operator*(float lhs, const mat3& rhs)
	{
		return mat3(lhs * rhs[0], lhs * rhs[1], lhs * rhs[2]);
	}

	inline mat3 operator/(float lhs, const mat3& rhs)
	{
		return mat3(lhs / rhs[0], lhs / rhs[1], lhs / rhs[2]);
	}

	inline mat3 transpose(const mat3& arg)
	{
		return mat3(vec3(arg(0, 0), arg(1, 0), arg(2, 0)),
			vec3(arg(0, 1), arg(1, 1), arg(2, 1)),
			vec3(arg(0, 2), arg(1, 2), arg(2, 2)));
	}

	//! matrix-vector multiplication
	inline vec3 operator*(const mat3& lhs, const vec3& rhs)
	{
		return vec3(dot(lhs[0], rhs), dot(lhs[1], rhs), dot(lhs[2], rhs));
	}

	//! matrix-matrix multiplication
	inline mat3 mult(const mat3& lhs, const mat3& rhs)
	{
		mat3 rhsTransposed = transpose(rhs);
		// Here we use the equivalence:
		// (R^T * L^T)^T = (L^T)^T * (R^T)^T = L * R
		return mat3(rhsTransposed * lhs[0], rhsTransposed * lhs[1], rhsTransposed * lhs[2]);
	}

	//! matrix determinant
	inline float det(const mat3& arg)
	{
		// recursive minor decomposition
		return dot(arg[0], cross(arg[1], arg[2]));
	}

	// does not check whether the inverse actually exists
	inline mat3 inverse(const mat3& arg)
	{
		// refer to: http://mathworld.wolfram.com/MatrixInverse.html
		// compute minors' determinants:
		vec3 col0 = cross(arg[1], arg[2]);
		vec3 col1 = cross(arg[2], arg[0]);
		vec3 col2 = cross(arg[0], arg[1]);
		float det = dot(arg[0], col0);
		return fromColumns(col0, col1, col2) / det;
	}

	mat3 rotationFromXYZ(const vec3& angles)
	{
		float c1 = cosf(angles[0]);
		float c2 = cosf(angles[1]);
		float c3 = cosf(angles[2]);
		float s1 = sinf(angles[0]);
		float s2 = sinf(angles[1]);
		float s3 = sinf(angles[2]);
		float s12 = s1 * s2;
		return mat3(vec3(c2*c3, -c2 * s3, s2),
			vec3(c1*s3 + c3 * s12, c1*c3 - s12 * s3, -c2 * s1),
			vec3(s1*s3 - c1 * c3*s2, c3*s1 + c1 * s2*s3, c1*c2));
	}

	struct mat4
	{
		union
		{
			struct
			{
				vec4 row0, row1, row2, row3;
			};
			vec4 r[4];
			float e[4][4];
			float f[16]; // flat
		};

		mat4() {}
		explicit mat4(float scalar) : row0(scalar), row1(scalar), row2(scalar), row3(scalar) {}
		mat4(float e00, float e01, float e02, float e03,
			float e10, float e11, float e12, float e13,
			float e20, float e21, float e22, float e23,
			float e30, float e31, float e32, float e33)
			: row0(e00, e01, e02, e03),
			row1(e10, e11, e12, e13),
			row2(e20, e21, e22, e23),
			row3(e30,e31,e32,e33)
		{}
		mat4(const vec4& row0, const vec4& row1, const vec4& row2, const vec4& row3)
			: row0(row0), row1(row1), row2(row2), row3(row3) {}
		mat4(const mat4& arg) : row0(arg.row0), row1(arg.row1), row2(arg.row2), row3(arg.row3) {}

		mat4& operator=(const mat4& arg)
		{
			row0 = arg.row0;
			row1 = arg.row1;
			row2 = arg.row2;
			row3 = arg.row3;
			return *this;
		}

		vec4 getColumn(u32 colIdx) const
		{
			assert(colIdx < 4 && "HydrogenCG::mat4::setColumn():: Column index out of bounds");
			return vec4(e[0][colIdx], e[1][colIdx], e[2][colIdx], e[3][colIdx]);
		}

		mat4& setColumn(u32 colIdx, const vec4& col)
		{
			assert(colIdx < 4 && "HydrogenCG::mat4::setColumn():: Column index out of bounds");
			e[0][colIdx] = col[0];
			e[1][colIdx] = col[1];
			e[2][colIdx] = col[2];
			e[3][colIdx] = col[3];

			return *this;
		}

		mat4& fromColumns(const vec4& col0, const vec4& col1, const vec4& col2 ,const vec4& col3)
		{
			row0 = vec4(col0[0], col1[0], col2[0], col3[0]);
			row1 = vec4(col0[1], col1[1], col2[1], col3[1]);
			row2 = vec4(col0[2], col1[2], col2[2], col3[2]);
			row2 = vec4(col0[3], col1[3], col2[3], col3[3]);
			return *this;
		}

		vec4& at(u32 rowIdx) { assert(rowIdx < 4 && "HydrogenCG::mat4::at():: Row index out of bounds"); return r[rowIdx]; }
		const vec4& at(u32 rowIdx) const { assert(rowIdx < 4 && "HydrogenCG::mat4::at():: Row index out of bounds"); return r[rowIdx]; }
		vec4& operator[](u32 rowIdx) { return r[rowIdx]; }
		const vec4& operator[](u32 rowIdx) const { return r[rowIdx]; }
		vec4& operator()(u32 rowIdx) { return r[rowIdx]; }
		const vec4& operator()(u32 rowIdx) const { return r[rowIdx]; }

		float& at(u32 rowIdx, u32 colIdx)
		{
			assert(rowIdx < 4 && "HydrogenCG::mat4::at():: Row index out of bounds.");
			assert(colIdx < 4 && "HydrogenCG::mat4::at():: Column index out of bounds.");
			return e[rowIdx][colIdx];
		}
		const float& at(u32 rowIdx, u32 colIdx) const
		{
			assert(rowIdx < 4 && "HydrogenCG::mat4::at():: Row index out of bounds.");
			assert(colIdx < 4 && "HydrogenCG::mat4::at():: Column index out of bounds.");
			return e[rowIdx][colIdx];
		}
		float& operator()(u32 rowIdx, u32 colIdx) { return r[rowIdx][colIdx]; }
		const float& operator()(u32 rowIdx, u32 colIdx) const { return r[rowIdx][colIdx]; }


		const mat4& operator+=(const mat4& rhs)
		{
			row0 += rhs.row0;
			row1 += rhs.row1;
			row2 += rhs.row2;
			row3 += rhs.row3;
			return *this;
		}

		const mat4& operator-=(const mat4& rhs)
		{
			row0 -= rhs.row0;
			row1 -= rhs.row1;
			row2 -= rhs.row2;
			row3 -= rhs.row3;
			return *this;
		}

		// Hadamard product
		const mat4& operator*=(const mat4& rhs)
		{
			row0 *= rhs.row0;
			row1 *= rhs.row1;
			row2 *= rhs.row2;
			row3 *= rhs.row3;
			return *this;
		}

		// Hadamard division
		const mat4& operator/=(const mat4& rhs)
		{
			row0 /= rhs.row0;
			row1 /= rhs.row1;
			row2 /= rhs.row2;
			row3 /= rhs.row3;
			return *this;
		}

		// Scalar operations
		const mat4& operator*=(float scalar)
		{
			row0 *= scalar;
			row1 *= scalar;
			row2 *= scalar;
			row3 *= scalar;
			return *this;
		}

		const mat4& operator/=(float scalar)
		{
			float invScalar = 1 / scalar;
			row0 *= invScalar;
			row1 *= invScalar;
			row2 *= invScalar;
			row3 *= invScalar;
			return *this;
		}

		static const mat4 identity;
		static const mat4 flipHandedness;

		static mat4 scale(const vec3& scaleFactors)
		{
			return mat4(vec4(scaleFactors[0], 0, 0, 0),
				vec4(0, scaleFactors[1], 0, 0),
				vec4(0, 0, scaleFactors[2], 0),
				vec4(0, 0, 0, 1));
		}
	};

	const mat4 mat4::identity = mat4(vec4(1, 0, 0, 0), vec4(0, 1, 0, 0), vec4(0, 0, 1, 0), vec4(0, 0, 0, 1));
	const mat4 mat4::flipHandedness = mat4(vec4(1, 0, 0, 0), vec4(0, 1, 0, 0), vec4(0, 0, 1, 0), vec4(0, 0, 0, -1));
	/*mat4 rotationFromXYZ(const vec3& angles)
	{
		float c1 = cosf(angles[0]);
		float c2 = cosf(angles[1]);
		float c3 = cosf(angles[2]);
		float s1 = sinf(angles[0]);
		float s2 = sinf(angles[1]);
		float s3 = sinf(angles[2]);
		float s12 = s1 * s2;
		return mat4(vec4(c2*c3, -c2 * s3, s2, 0),
			vec4(c1*s3 + c3 * s12, c1*c3 - s12 * s3, -c2 * s1, 0),
			vec4(s1*s3 - c1 * c3*s2, c3*s1 + c1 * s2*s3, c1*c2, 0),
			vec4(0, 0, 0, 1));
	}*/

	

	mat4 fromColumns(const vec4& col0, const vec4& col1, const vec4& col2, const vec4& col3)
	{
		return mat4(vec4(col0[0], col1[0], col2[0], col3[0]),
			vec4(col0[1], col1[1], col2[1], col3[1]),
			vec4(col0[2], col1[2], col2[2], col3[2]),
			vec4(col0[3], col1[3], col2[3], col3[3]));
	}

	inline mat4 operator+(const mat4& arg)
	{
		return arg;
	}

	inline mat4 operator-(const mat4& arg)
	{
		return mat4(-arg[0], -arg[1], -arg[2], -arg[3]);
	}

	// Componentwise / Hadamard operations

	inline mat4 operator+(const mat4& lhs, const mat4& rhs)
	{
		return mat4(lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2], lhs[3] + rhs[3]);
	}

	inline mat4 operator-(const mat4& lhs, const mat4& rhs)
	{
		return mat4(lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2], lhs[3] - rhs[3]);
	}

	inline mat4 operator*(const mat4& lhs, const mat4& rhs)
	{
		return mat4(lhs[0] * rhs[0], lhs[1] * rhs[1], lhs[2] * rhs[2], lhs[3] * rhs[3]);
	}

	inline mat4 operator/(const mat4& lhs, const mat4& rhs)
	{
		return mat4(lhs[0] / rhs[0], lhs[1] / rhs[1], lhs[2] / rhs[2], lhs[3] / rhs[3]);
	}

	// Scalar operations

	inline mat4 operator*(const mat4& lhs, float rhs)
	{
		return mat4(lhs[0] * rhs, lhs[1] * rhs, lhs[2] * rhs, lhs[3] * rhs);
	}

	inline mat4 operator/(const mat4& lhs, float rhs)
	{
		float invRhs = 1.0f / rhs; // note: precision is not as good as dividing each component
		return mat4(lhs[0] * invRhs, lhs[1] * invRhs, lhs[2] * invRhs, lhs[3] * invRhs);
	}

	inline mat4 operator*(float lhs, const mat4& rhs)
	{
		return mat4(lhs * rhs[0], lhs * rhs[1], lhs * rhs[2], lhs * rhs[3]);
	}

	inline mat4 operator/(float lhs, const mat4& rhs)
	{
		return mat4(lhs / rhs[0], lhs / rhs[1], lhs / rhs[2], lhs / rhs[3]);
	}

	inline mat4 transpose(const mat4& arg)
	{
		return mat4(vec4(arg(0, 0), arg(1, 0), arg(2, 0), arg(3,0)),
			vec4(arg(0, 1), arg(1, 1), arg(2, 1), arg(3,1)),
			vec4(arg(0, 2), arg(1, 2), arg(2, 2), arg(3,2)),
			vec4(arg(0, 3), arg(1, 3), arg(2, 3), arg(3,3)));
	}

	//! matrix-vector multiplication
	inline vec4 operator*(const mat4& lhs, const vec4& rhs)
	{
		return vec4(dot(lhs[0], rhs), dot(lhs[1], rhs), dot(lhs[2], rhs), dot(lhs[3], rhs));
	}

	//! matrix-matrix multiplication
	inline mat4 mult(const mat4& lhs, const mat4& rhs)
	{
		mat4 rhsTransposed = transpose(rhs);
		// Here we use the equivalence:
		// (R^T * L^T)^T = (L^T)^T * (R^T)^T = L * R
		return mat4(rhsTransposed * lhs[0], rhsTransposed * lhs[1], rhsTransposed * lhs[2], rhsTransposed * lhs[3]);
	}

	// https://stackoverflow.com/questions/1148309/inverting-a-4x4-matrix
	// I know it's ugly
	inline float det(const mat4& arg)
	{

		float m0 = arg.f[5] * arg.f[10] * arg.f[15] -
			arg.f[5] * arg.f[11] * arg.f[14] -
			arg.f[9] * arg.f[6] * arg.f[15] +
			arg.f[9] * arg.f[7] * arg.f[14] +
			arg.f[13] * arg.f[6] * arg.f[11] -
			arg.f[13] * arg.f[7] * arg.f[10];

		float m1 = -arg.f[4] * arg.f[10] * arg.f[15] +
			arg.f[4] * arg.f[11] * arg.f[14] +
			arg.f[8] * arg.f[6] * arg.f[15] -
			arg.f[8] * arg.f[7] * arg.f[14] -
			arg.f[12] * arg.f[6] * arg.f[11] +
			arg.f[12] * arg.f[7] * arg.f[10];

		float m2 = arg.f[4] * arg.f[9] * arg.f[15] -
			arg.f[4] * arg.f[11] * arg.f[13] -
			arg.f[8] * arg.f[5] * arg.f[15] +
			arg.f[8] * arg.f[7] * arg.f[13] +
			arg.f[12] * arg.f[5] * arg.f[11] -
			arg.f[12] * arg.f[7] * arg.f[9];

		float m3 = -arg.f[4] * arg.f[9] * arg.f[14] +
			arg.f[4] * arg.f[10] * arg.f[13] +
			arg.f[8] * arg.f[5] * arg.f[14] -
			arg.f[8] * arg.f[6] * arg.f[13] -
			arg.f[12] * arg.f[5] * arg.f[10] +
			arg.f[12] * arg.f[6] * arg.f[9];

		return arg.f[0] * m0 + arg.f[1] * m1 + arg.f[2] * m2 + arg.f[3] * m3;
	}

	// https://stackoverflow.com/questions/1148309/inverting-a-4x4-matrix
	// does not check whether the inverse actually exists
	inline mat4 inverse(const mat4& arg)
	{
		mat4 inv;
		float det;

		inv.f[0] = arg.f[5] * arg.f[10] * arg.f[15] -
			arg.f[5] * arg.f[11] * arg.f[14] -
			arg.f[9] * arg.f[6] * arg.f[15] +
			arg.f[9] * arg.f[7] * arg.f[14] +
			arg.f[13] * arg.f[6] * arg.f[11] -
			arg.f[13] * arg.f[7] * arg.f[10];

		inv.f[4] = -arg.f[4] * arg.f[10] * arg.f[15] +
			arg.f[4] * arg.f[11] * arg.f[14] +
			arg.f[8] * arg.f[6] * arg.f[15] -
			arg.f[8] * arg.f[7] * arg.f[14] -
			arg.f[12] * arg.f[6] * arg.f[11] +
			arg.f[12] * arg.f[7] * arg.f[10];

		inv.f[8] = arg.f[4] * arg.f[9] * arg.f[15] -
			arg.f[4] * arg.f[11] * arg.f[13] -
			arg.f[8] * arg.f[5] * arg.f[15] +
			arg.f[8] * arg.f[7] * arg.f[13] +
			arg.f[12] * arg.f[5] * arg.f[11] -
			arg.f[12] * arg.f[7] * arg.f[9];

		inv.f[12] = -arg.f[4] * arg.f[9] * arg.f[14] +
			arg.f[4] * arg.f[10] * arg.f[13] +
			arg.f[8] * arg.f[5] * arg.f[14] -
			arg.f[8] * arg.f[6] * arg.f[13] -
			arg.f[12] * arg.f[5] * arg.f[10] +
			arg.f[12] * arg.f[6] * arg.f[9];

		inv.f[1] = -arg.f[1] * arg.f[10] * arg.f[15] +
			arg.f[1] * arg.f[11] * arg.f[14] +
			arg.f[9] * arg.f[2] * arg.f[15] -
			arg.f[9] * arg.f[3] * arg.f[14] -
			arg.f[13] * arg.f[2] * arg.f[11] +
			arg.f[13] * arg.f[3] * arg.f[10];

		inv.f[5] = arg.f[0] * arg.f[10] * arg.f[15] -
			arg.f[0] * arg.f[11] * arg.f[14] -
			arg.f[8] * arg.f[2] * arg.f[15] +
			arg.f[8] * arg.f[3] * arg.f[14] +
			arg.f[12] * arg.f[2] * arg.f[11] -
			arg.f[12] * arg.f[3] * arg.f[10];

		inv.f[9] = -arg.f[0] * arg.f[9] * arg.f[15] +
			arg.f[0] * arg.f[11] * arg.f[13] +
			arg.f[8] * arg.f[1] * arg.f[15] -
			arg.f[8] * arg.f[3] * arg.f[13] -
			arg.f[12] * arg.f[1] * arg.f[11] +
			arg.f[12] * arg.f[3] * arg.f[9];

		inv.f[13] = arg.f[0] * arg.f[9] * arg.f[14] -
			arg.f[0] * arg.f[10] * arg.f[13] -
			arg.f[8] * arg.f[1] * arg.f[14] +
			arg.f[8] * arg.f[2] * arg.f[13] +
			arg.f[12] * arg.f[1] * arg.f[10] -
			arg.f[12] * arg.f[2] * arg.f[9];

		inv.f[2] = arg.f[1] * arg.f[6] * arg.f[15] -
			arg.f[1] * arg.f[7] * arg.f[14] -
			arg.f[5] * arg.f[2] * arg.f[15] +
			arg.f[5] * arg.f[3] * arg.f[14] +
			arg.f[13] * arg.f[2] * arg.f[7] -
			arg.f[13] * arg.f[3] * arg.f[6];

		inv.f[6] = -arg.f[0] * arg.f[6] * arg.f[15] +
			arg.f[0] * arg.f[7] * arg.f[14] +
			arg.f[4] * arg.f[2] * arg.f[15] -
			arg.f[4] * arg.f[3] * arg.f[14] -
			arg.f[12] * arg.f[2] * arg.f[7] +
			arg.f[12] * arg.f[3] * arg.f[6];

		inv.f[10] = arg.f[0] * arg.f[5] * arg.f[15] -
			arg.f[0] * arg.f[7] * arg.f[13] -
			arg.f[4] * arg.f[1] * arg.f[15] +
			arg.f[4] * arg.f[3] * arg.f[13] +
			arg.f[12] * arg.f[1] * arg.f[7] -
			arg.f[12] * arg.f[3] * arg.f[5];

		inv.f[14] = -arg.f[0] * arg.f[5] * arg.f[14] +
			arg.f[0] * arg.f[6] * arg.f[13] +
			arg.f[4] * arg.f[1] * arg.f[14] -
			arg.f[4] * arg.f[2] * arg.f[13] -
			arg.f[12] * arg.f[1] * arg.f[6] +
			arg.f[12] * arg.f[2] * arg.f[5];

		inv.f[3] = -arg.f[1] * arg.f[6] * arg.f[11] +
			arg.f[1] * arg.f[7] * arg.f[10] +
			arg.f[5] * arg.f[2] * arg.f[11] -
			arg.f[5] * arg.f[3] * arg.f[10] -
			arg.f[9] * arg.f[2] * arg.f[7] +
			arg.f[9] * arg.f[3] * arg.f[6];

		inv.f[7] = arg.f[0] * arg.f[6] * arg.f[11] -
			arg.f[0] * arg.f[7] * arg.f[10] -
			arg.f[4] * arg.f[2] * arg.f[11] +
			arg.f[4] * arg.f[3] * arg.f[10] +
			arg.f[8] * arg.f[2] * arg.f[7] -
			arg.f[8] * arg.f[3] * arg.f[6];

		inv.f[11] = -arg.f[0] * arg.f[5] * arg.f[11] +
			arg.f[0] * arg.f[7] * arg.f[9] +
			arg.f[4] * arg.f[1] * arg.f[11] -
			arg.f[4] * arg.f[3] * arg.f[9] -
			arg.f[8] * arg.f[1] * arg.f[7] +
			arg.f[8] * arg.f[3] * arg.f[5];

		inv.f[15] = arg.f[0] * arg.f[5] * arg.f[10] -
			arg.f[0] * arg.f[6] * arg.f[9] -
			arg.f[4] * arg.f[1] * arg.f[10] +
			arg.f[4] * arg.f[2] * arg.f[9] +
			arg.f[8] * arg.f[1] * arg.f[6] -
			arg.f[8] * arg.f[2] * arg.f[5];

		det = arg.f[0] * inv.f[0] + arg.f[1] * inv.f[4] + arg.f[2] * inv.f[8] + arg.f[3] * inv.f[12];

		float invDet = 1.0f / det;

		return inv * invDet;
	}
}