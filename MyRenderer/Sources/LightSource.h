// ----------------------------------------------
// Polytechnique - INF584 "Image Synthesis"
//
// Base code for practical assignments.
//
// Copyright (C) 2022 Tamy Boubekeur
// All rights reserved.
// ----------------------------------------------
#pragma once

#include <glm/glm.hpp>
#include <glm/ext.hpp>

//MOVE !
#include <random>

class LightSource {
public:
	LightSource(const glm::vec3& direction = glm::vec3 (0.f, 0.f, -1.f), const glm::vec3& color = glm::vec3 (1.f, 1.f, 1.f), float intensity = 1.f, bool area = true) : m_direction(direction), m_color(color), m_intensity(intensity), m_is_area (area) {}

	inline const glm::vec3& direction() const { return m_direction; }

	inline const glm::vec3& color() const { return m_color; }

	inline float intensity() const { return m_intensity; }

	//move to cpp
	inline glm::vec3& sample_point() const { 
		//functions inspired from
		//https://github.com/vascofazza/computer_graphics-Pathtrace/blob/master/src/pathtrace.cpp
		//TO DO: make it uniform using triangle areas
		int s = area->triangleIndices().size();
		int random_i = 0;
		glm::uvec3 random_triangle = area->triangleIndices()[random_i];
		return get_random_point_triangle(random_triangle);
	}

	int random_int(int min, int max){
   		return rand()%(max-min+1) + min;
	}

	//move to cpp
	inline glm::vec3& get_random_point_triangle(glm::uvec3 random_triangle) const
	{
	//https://math.stackexchange.com/questions/18686/uniform-random-point-in-triangle-in-3d
	glm::vec3 A = area->vertexPositions()[random_triangle[0]];
	glm::vec3 B = area->vertexPositions()[random_triangle[1]];
	glm::vec3 C = area->vertexPositions()[random_triangle[2]];
	//move
	std::random_device rd;
	std::default_random_engine eng(rd());
	std::uniform_real_distribution<> distr(0.f, 1.f);

	float r1 = (distr(eng));
	float r2 = (distr(eng));
	return (1-sqrt(r1))*A + (sqrt(r1)*(1-r2))*B + (r2*sqrt(r1))*C;
	}

	bool is_area() {
		return m_is_area;
	}

	//move to private
	std::shared_ptr<Mesh> area;

private:
	glm::vec3 m_direction;
	glm::vec3 m_color;
	float m_intensity;
	bool m_is_area;
	
	

};