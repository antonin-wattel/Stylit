// ----------------------------------------------
// Polytechnique - INF584 "Image Synthesis"
//
// Base code for practical assignments.
//
// Copyright (C) 2022 Tamy Boubekeur
// All rights reserved.
// ----------------------------------------------
#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>

#include <glm/glm.hpp>
#include <glm/ext.hpp>

class Image_multichannel {
public:
	inline Image_multichannel(size_t width = 64, size_t height = 64, size_t num_channels = 5 ) :
		m_width (width),
		m_height (height),
		m_num_channels(num_channels) {
		m_pixels.resize (width*height*m_num_channels, glm::vec3 (0.f, 0.f, 0.f));
	}

	inline virtual ~Image_multichannel() {}

	inline size_t width () const { return m_width; }

	inline size_t height () const { return m_height; }

	inline size_t num_channels() const { return m_num_channels; }

	inline const glm::vec3 & operator() (size_t channel, size_t x, size_t y) const { return m_pixels[m_num_channels*(y*m_width+x) + channel]; }
	inline glm::vec3 & operator() (size_t channel, size_t x, size_t y) { return m_pixels[m_num_channels*(y*m_width+x) + channel]; }

	inline const glm::vec3 & operator[] (size_t i) const { return m_pixels[i]; }

	inline glm::vec3 & operator[] (size_t i) { return m_pixels[i]; }

	inline const std::vector<glm::vec3> & pixels () const { return m_pixels; }

	/// Clear to 'color', black by default.
	inline void clear (const glm::vec3 & color = glm::vec3 (0.f, 0.f, 0.f)) {
		for (size_t channel = 0; channel < m_num_channels; channel++)
			for (size_t y = 0; y < m_height; y++)
				for (size_t x = 0; x < m_width; x++) 
					m_pixels[m_num_channels*(y*m_width+x)+channel] = color;
	}

	inline void savePPM (const std::string & filename, size_t channel) {
		std::ofstream out (filename.c_str ());
    	if (!out) {
        	std::cerr << "Cannot open file " << filename.c_str() << std::endl;
        	std::exit (1);
    	}
    	out << "P3" << std::endl
    		<< m_width << " " << m_height << std::endl
    		<< "255" << std::endl;
    	for (size_t y = 0; y < m_height; y++)
			for (size_t x = 0; x < m_width; x++) {
				out << std::min (255u, static_cast<unsigned int> (255.f * m_pixels[m_num_channels * (y * m_width + x) + channel][0])) << " "
					<< std::min (255u, static_cast<unsigned int> (255.f * m_pixels[m_num_channels * (y * m_width + x) + channel][1])) << " "
					<< std::min (255u, static_cast<unsigned int> (255.f * m_pixels[m_num_channels * (y * m_width + x) + channel][2])) << " ";
			}
			out << std::endl;
    	out.close ();
	}

	void save (const std::string & filename, size_t channel) const;
	void save_all_channels(const std::string& filename) const;

	glm::mat3 get_patch(glm::vec2 p) {
		//get neighborhood square patch at pixel p
		//returns a matrix of size ??
		float px = p.x;
		float py = p.y;

		//Eigen::Matrix(3, 3, m_num_channels*3)
		////indices
		//(px - 1, py + 1);
		//(px, py + 1);
		//(px + 1, py + 1);
		//(px - 1, py );
		//(px, py);
		//(px + 1, py);
		//(px - 1, py - 1);
		//(px , py - 1);
		//(px + 1, py - 1);

	}




private:
	size_t m_width;
	size_t m_height;
	size_t m_num_channels;

	std::vector<glm::vec3> m_pixels;
};