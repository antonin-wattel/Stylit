#pragma once

#include "Image.h"
#include <glm/glm.hpp>
#include <glm/ext.hpp>

#include "Image_multichannel.h"
#include "Image.h"

using namespace std;

//TO DO: MAKE A CLASS

void run_stylit(std::pair<Image_multichannel, Image>&A, std::pair<Image_multichannel, Image>& B);
float error(const pair<Image_multichannel, Image> & A, const pair<Image_multichannel, Image> & B, const glm::vec2 &p, const glm::vec2 &q, const float &mu);
glm::vec2 get_NNF(const glm::vec2 &q, const pair<Image_multichannel, Image>& A, const pair<Image_multichannel, Image>& B, const float & mu);
glm::vec3 average(const pair<Image_multichannel, Image>& A, const std::vector<glm::vec2>& NNF, const glm::vec2 & q);

