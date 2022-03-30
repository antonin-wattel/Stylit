#pragma once

#include "Image.h"
#include <glm/glm.hpp>
#include <glm/ext.hpp>

#include "Image_multichannel.h"
#include "Image.h"

using namespace std;

//TO DO: MAKE A CLASS

void run_stylit(std::pair<Image_multichannel, Image> A, std::pair<Image_multichannel, Image> B);
float error(pair<Image_multichannel, Image> A, pair<Image_multichannel, Image> B, glm::vec2 p, glm::vec2 q, float mu);
glm::vec2 get_NNF(glm::vec2 q, pair<Image_multichannel, Image> A, pair<Image_multichannel, Image> B, float mu);
glm::vec3 average(pair<Image_multichannel, Image> A, std::vector<glm::vec2> NNF, glm::vec2 q);

