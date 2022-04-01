#pragma once

#include "Image_multichannel.h"
#include <vector>


#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>


class Pyramid {

	public:
		inline Pyramid(const size_t &depth, const Image_multichannel &A) {
			//build the pyramid
			std::cout << "build" << std::endl;
			m_pyramid.reserve(depth); //??
			std::cout << "resized" << std::endl;
			m_depth = depth;
			m_pyramid.push_back(A);
			std::cout << "push 1" << std::endl;


			for (int d = 1; d < m_depth; d++) {
				//build level d
				/*Image_multichannel A_d(width, height, A.num_channels());*/
				std::cout << "d = " << d << std::endl;
				Image_multichannel A_filtered = gaussian_filter(m_pyramid[d-1]);
				std::cout << "down" << std::endl;
				m_pyramid.push_back(A_filtered.downsample(2));
				std::cout << "pushback" << std::endl;
				//gaussian filter
				//downsample
			}
			std::cout << "built pyramid" << std::endl;
		}

		inline virtual ~Pyramid() {}

		inline std::vector<Image_multichannel> const pyramid() { return m_pyramid; }

		inline size_t depth() const { return m_depth; }

		inline const Image_multichannel& operator[] (size_t i) const { return m_pyramid[i]; }

		inline Image_multichannel& operator[] (size_t i) { return m_pyramid[i]; }

		
		//const Image_multichannel Pyramid::upsample(const Image_multichannel& A);
		const Image_multichannel gaussian_filter(const Image_multichannel &A);

	private:

		size_t downsampling;
		size_t m_depth;
		std::vector<Image_multichannel> m_pyramid;

};