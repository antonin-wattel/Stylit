#pragma once

#include "Image.h"
#include <glm/glm.hpp>
#include <glm/ext.hpp>

#include "Image_multichannel.h"
#include "Image.h"

#include "Pyramid.h"

using namespace std;

//TO DO: MAKE A CLASS

void run_stylit(std::pair<Image_multichannel, Image_multichannel>&A, std::pair<Image_multichannel, Image_multichannel>& B);
float error(const pair<Image_multichannel, Image_multichannel> & A, const pair<Image_multichannel, Image_multichannel> & B, const glm::vec2 &p, const glm::vec2 &q, const float &mu);
glm::vec2 get_NNF(const glm::vec2 &q, const pair<Image_multichannel, Image_multichannel>& A, const pair<Image_multichannel, Image_multichannel>& B, const float & mu);
glm::vec3 average(const pair<Image_multichannel, Image_multichannel>& A, const std::vector<glm::vec2>& NNF, const glm::vec2 & q);

