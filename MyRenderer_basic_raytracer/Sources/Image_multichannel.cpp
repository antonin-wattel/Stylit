#pragma once

//#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

#include "Image_multichannel.h"

#include "Image.h"
#include "Console.h"

void Image_multichannel::save (const std::string & filename, size_t channel) const {
	unsigned char * tmp = new unsigned char[3*m_width*m_height];
	stbi_write_png_compression_level = 1;
	for (unsigned int y = 0; y < m_height; y++)
#pragma omp parallel for
		for (int x = 0; x < m_width; x++)
			for (unsigned int k = 0; k < 3; k++)
				tmp[3*(y*m_width+x)+k] = 255*glm::clamp (this->operator()(channel, x,m_height-1-y)[k], 0.f, 1.f);
	stbi_write_png(filename.c_str(), m_width, m_height, 3, tmp, 3 * m_width);
	delete [] tmp;
}


void Image_multichannel::save_all_channels(const std::string& filename) const {
	for (int i = 0; i < m_num_channels; i++) {
		std::string filename_channel = filename + std::to_string(i) + ".png";
		this->save(filename_channel, i);
	}
}