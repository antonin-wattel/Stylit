#pragma once

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

#include "Image.h"
#include "Console.h"


//#define STB_IMAGE_IMPLEMENTATION
//#include "stb_image.h"
//
//Image::Image (const std::string& filename) {
//		std::cout<<"reading_image"<<std::endl;
//		int width, height, bpp;
//		unsigned char* rgb_image = stbi_load(filename.c_str(), &width, &height, &bpp, 3);
//
//		m_width = width;
//		m_height = height;
//		m_pixels.resize (width*height, glm::vec3 (0.f, 0.f, 0.f));
//		std::cout<<"width = "<< width;
//		std::cout<<"height = "<< height;
//		std::cout<<"bpp = "<< height;
//
//		float i = 0.1;
//		for (unsigned int y = 0; y < height; y++){
//			for (int x = 0; x < width; x++){
//				i+=0.00001;
//				for (unsigned int k = 0; k < 3; k++){
//					this->operator()(x,height-1-y)[k] = (float)rgb_image[3*(y*width+x)+k]/255.;
//					//std::cout<<(int)rgb_image[3*(y*width+x)+k]<<" "<<std::endl;
//				}
//			}
//		}
//		stbi_image_free(rgb_image);
//}

void Image::save (const std::string & filename) const {
	unsigned char * tmp = new unsigned char[3*m_width*m_height];
	stbi_write_png_compression_level = 1;
	for (unsigned int y = 0; y < m_height; y++)
#pragma omp parallel for
		for (int x = 0; x < m_width; x++)
			for (unsigned int k = 0; k < 3; k++)
				tmp[3*(y*m_width+x)+k] = 255*glm::clamp (this->operator()(x,m_height-1-y)[k], 0.f, 1.f);
	stbi_write_png(filename.c_str(), m_width, m_height, 3, tmp, 3 * m_width);
	delete [] tmp;
}
