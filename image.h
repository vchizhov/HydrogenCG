#pragma once
#include "vec.h"
#include "math.h"
#include "typedef.h"
#include <cstdint>
#include <fstream>


namespace HydrogenCG
{
	/*!
		\brief		A class to represent a 3-component float valued image.

		Provides a method to save to PPM.
		We use row-major order for the array, so accessing elements 
		should be done along the x coordinate first preferably.
	*/
	class Image
	{
	private:
		vec3* e;		//!< The position of the the aperture of the camera
		u32 width;		//!< Number of columns
		u32 height;		//!< Number of rows

	public:
		Image() : e(nullptr), width(0), height(0) {}
		Image(u32 w, u32 h) { init(w, h); }
		Image(const Image& arg) { init(arg.w(), arg.h()); copy(arg); };
		~Image() { free(); }

		Image& operator=(const Image& arg)
		{
			copy(arg);
			return *this;
		}

		void init(u32 w, u32 h)
		{
			free();
			e = new vec3[w*h];
			width = w;
			height = h;
		}

		void copy(const Image& arg)
		{
			assert(w() == arg.w() && h() == arg.h() && "Image::copy():: Image dimensions do not match.");
			memcpy(e, arg.data(), arg.size() * sizeof(vec3));
		}

		void free()
		{
			delete[] e;
			width = 0;
			height = 0;
		}

		inline u32 w() const { return width; }
		inline u32 h() const { return height; }
		inline u32 size() const { return w() * h(); }
		const vec3* data() const { return e; }
		const float* dataPlain() const { return reinterpret_cast<float*>(e); }

		/*!
			Given the x (column) and y (row) index in the 2D array, returns the corresponding index
			in the 1D array, note that this defines the ordering of the array, in this case row-major
		*/
		inline u32 idx(u32 x, u32 y) const { return x + y * w(); }

		/*! 
			Given an index into the 1D array, returns the corresponding x index (column) in the 2D array,
			should be consistent with idx()
		*/
		inline u32 xIdx(u32 idx) const { return idx % w(); }

		/*! 
			Given an index into the 1D array, returns the corresponding y index (row) in the 2D array,
			should be consistent with idx()
		*/
		inline u32 yIdx(u32 idx) const { return idx / w(); }

		// accessors:

		// Convenience (rather than calling (*this)(idx) when defining other class methods)
		inline vec3& get(u32 idx) { return e[idx]; }
		inline const vec3& get(u32 idx) const { return e[idx]; }
		inline vec3& get(u32 x, u32 y) { return e[idx(x, y)]; }
		inline const vec3& get(u32 x, u32 y) const { return e[idx(x, y)]; }

		// accessors with bounds checks
		inline vec3& at(u32 idx) { assert(idx < size() && "Array index out of bounds."); return get(idx); }
		inline const vec3& at(u32 idx) const { assert(idx < size() && "Array index out of bounds."); return get(idx); }
		inline vec3& at(u32 x, u32 y) { assert(x < w() && y < h() && "Array index out of bounds."); return get(x,y); }
		inline const vec3& at(u32 x, u32 y) const { assert(x < w() && y < h() && "Array index out of bounds."); return get(x, y); }

		// accessors without bound checks
		inline vec3& operator()(u32 idx) { return get(idx); }
		inline const vec3& operator()(u32 idx) const { return get(idx); }
		inline vec3& operator()(u32 x, u32 y) { return get(x, y); }
		inline const vec3& operator()(u32 x, u32 y) const { return get(x, y); }

		inline vec3& operator[](u32 idx) { return get(idx); }
		inline const vec3& operator[](u32 idx) const { return get(idx); }
		
		
		/*!
			\brief		Saves the current image to a plain PPM file with [0,maxVal] range

			The default max value is set to 255, higher values may be used for higher 
			precision and less banding due to inadequate quantization

			see reference: http://netpbm.sourceforge.net/doc/ppm.html
		*/
		bool savePPM(const char* filename, float gamma = 1/2.2f, u16 maxVal = 255)
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
			file << maxVal << "\n"; // the maximum color value, can be at most 2^16-1 = 65535

			// Iterate over all rows
			for (u32 y = 0; y < h(); ++y)
			{
				// Save each row
				for (u32 x = 0; x < w(); ++x)
				{
					vec3 col = at(x, y);

					// throw away negative values
					col = max(col, 0.0f);
					// apply gamma
					col = pow(col, gamma);
					// map from [0,1] to [0,maxVal] and round
					col = round(float(maxVal)*col);
					// throw away values higher than maxVal
					col = min(col, float(maxVal));
					// convert to u32 and write to file
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