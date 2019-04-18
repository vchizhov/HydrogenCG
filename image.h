#pragma once
#include "vec.h"
#include "math.h"
#include <cstdint>
#include <fstream>


namespace HydrogenCG
{
	// row-major image - see definition of at()
	class Image
	{
	private:
		vec3* e;
		uint32_t width, height;
	public:
		Image() : e(nullptr), width(0), height(0) {}
		~Image() { free(); }

		void init(uint32_t w, uint32_t h)
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

		uint32_t w() const { return width; }
		uint32_t h() const { return height; }
		const vec3* data() const { return e; }

		vec3& at(uint32_t idx) { return e[idx]; }
		const vec3& at(uint32_t idx) const { return e[idx]; }

		vec3& operator()(uint32_t idx) { assert(idx < width*height && "Array index out of bounds"); return at(idx); }
		const vec3& operator()(uint32_t idx) const { assert(idx < width*height && "Array index out of bounds"); return at(idx); }

		vec3& at(uint32_t x, uint32_t y) { return e[x + width * y]; }
		const vec3& at(uint32_t x, uint32_t y) const { return e[x + width * y]; }
		vec3& operator()(uint32_t x, uint32_t y) { assert(x < width && y < height && "Array index out of bounds");  return at(x,y); }
		const vec3& operator()(uint32_t x, uint32_t y) const { assert(x < width && y < height && "Array index out of bounds");  return at(x,y); }

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
			for (uint32_t y = 0; y < h(); ++y)
			{
				// Save each row
				for (uint32_t x = 0; x < w(); ++x)
				{
					vec3 col = clamp(at(x, y) * 256.0, 0.0, 255.0); // map from [0,1] to [0,255] and clamp
					uint32_t r = (uint32_t)col.r;
					uint32_t g = (uint32_t)col.g;
					uint32_t b = (uint32_t)col.b;
					file << r << "\t" << g << "\t" << b << "\t\t";
				}
				file << "\t";
			}
			file.close();
			return true;
		}
	};
}