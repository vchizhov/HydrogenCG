#pragma once


#include <cstdint>
#include <fstream>

#include "vec.h"
#include "math.h"
#include "typedef.h"


namespace HydrogenCG
{
	// row-major image - see definition of at()
	class Image
	{
	private:
		vec3* e;
		u32 width, height;
	public:
		Image() : e(nullptr), width(0), height(0) {}
		~Image() { free(); }

		void init(u32 w, u32 h)
		{
			free();
			e = new vec3[w*h];
			width = w;
			height = h;
		}

		void free()
		{
			delete[] e;
			width = 0;
			height = 0;
		}

		u32 w() const { return width; }
		u32 h() const { return height; }
		const vec3* data() const { return e; }

		vec3& at(u32 idx) { return e[idx]; }
		const vec3& at(u32 idx) const { return e[idx]; }

		vec3& operator()(u32 idx) { assert(idx < width*height"Array index out of bounds"); return at(idx); }
		const vec3& operator()(u32 idx) const { assert(idx < width*height"Array index out of bounds"); return at(idx); }

		vec3& at(u32 x, u32 y) { return e[x + width * y]; }
		const vec3& at(u32 x, u32 y) const { return e[x + width * y]; }
		vec3& operator()(u32 x, u32 y) { assert(x < width && y < height && "Array index out of bounds");  return at(x, y); }
		const vec3& operator()(u32 x, u32 y) const { assert(x < width && y < height && "Array index out of bounds");  return at(x, y); }

		// http://netpbm.sourceforge.net/doc/ppm.html
		bool savePPM(const char* filename)
		{
			std::ofstream file;
			file.open(filename);
			if (!file.good())
			{
				return false;
			}
			// header of a plain PPM file
			file << "P3\n";	// magic number identifying the file type
			file << w() << "\t" << h() << "\n"; // width and height of the image
			file << "255\n"; // the maximum color value, can be at most 2^16-1 = 65535

			// Iterate over all rows
			for (u32 y = 0; y < h(); ++y)
			{
				// Save each row
				for (u32 x = 0; x < w(); ++x)
				{
					vec3 col = clamp(at(x, y) * 256.0, 0.0, 255.0); // map from [0,1] to [0,255] and clamp
					u32 r = (u32)col.r;
					u32 g = (u32)col.g;
					u32 b = (u32)col.b;
					file << r << "\t" << g << "\t" << b << "\t\t";
				}
				file << "\t";
			}
			file.close();
			return true;
		}
	};
}