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


Image_multichannel Image_multichannel::downsample(const int& f) const { //these should go in the Image_multichannel class
//adapted from
//https://www.geeksforgeeks.org/spatial-resolution-down-sampling-and-up-sampling-in-image-processing/

    std::cout << "downsampling " << std::endl;
    std::cout << "new width: " << this->width() / f << std::endl;
     std::cout << "new width: " << this->width() / f << std::endl;
    Image_multichannel res(this->width() / f, this->height() / f, this->num_channels());
    res.clear(glm::vec3(0., 0., 0.));

    for (int channel = 0; channel < res.num_channels(); channel++) {
        for (int x = 0; x < this->width(); x++) {
            for (int y = 0; y < this->height(); y++) {
                //careful with indices
                if ((0 <= (x / f) < res.width()) && (0 <= (y / f) < res.height()))
                    //std::cout << "oyeaaa !" << std::endl;
                    res.operator()(channel, x/f, y/f) = this->operator()(channel, x, y);
            }
        }
    }
    std::cout << "downsampled" << std::endl;

    return res;
}


Image_multichannel Image_multichannel::upsample(const int& f) const {
//adapted from
//https://www.geeksforgeeks.org/spatial-resolution-down-sampling-and-up-sampling-in-image-processing/

    Image_multichannel res(this->width() * f, this->height() * f, this->num_channels());
    res.clear(glm::vec3(0., 0., 0.));

    for (int channel = 0; channel < res.num_channels(); channel++) {
        for (int x = 0; x < res.width()-1; x+=f) {
            for (int y = 0; y < res.height()-1; y+=f) {
                //careful with indices
                //if ((0 <= (x / f) < this.width()) && (0 <= (y / f) < this.height()))
                    //std::cout << "oyeaaa !" << std::endl;
                    res.operator()(channel, x, y) = this->operator()(channel, x/f, y/f);
            }
        }
    }

    //replicating rows
    for (int channel = 0; channel < res.num_channels(); channel++) {
        for (int x = 1; x < (res.width() - (f-1)); x+=f) {
            for (int y = 0; y < (res.height() -(f-1)); y++) {
                for (int i = x; i < x + (f - 1); i++) 
                    res.operator()(channel, x, y) = res.operator()(channel, i-1, y);
            }
        }
    }

    for (int channel = 0; channel < res.num_channels(); channel++) {
        for (int x = 0; x < (res.width() -1); x ++) {
            for (int y = 1; y < (res.height() - 1); y+=f) {
                for (int j = y; j < y + (f - 1); j++)
                    res.operator()(channel, x, y) = res.operator()(channel, x, j-1);
            }
        }
    }

    return res;
} 
   
//TO DO : WRITE gaussian filter function