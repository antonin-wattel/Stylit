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

	

	inline glm::vec3& sample_point() const { 
		//method that uniformly samples positions on the surface mesh

		//https://www.cs.cornell.edu/courses/cs6630/2015fa/pa2.html -> is how it should be done !
		//emits radiance uniformly in all directions from each triangle of the mesh

		//sampling by first choosing a triangle from the mesh according to its area

		//then sampling uniformly within that triangle

		//functions inspired from
		//https://github.com/vascofazza/computer_graphics-Pathtrace/blob/master/src/pathtrace.cpp

		// //REMOVE THAT SHIT !
	
		// glm::vec3 bbox = area->triangleIndices();

		//TO DO: MAKE IT UNIFORM !
		int s = area->triangleIndices().size();
		//int random_i = rand()%(s+1);
		int random_i = 0;
		glm::uvec3 random_triangle = area->triangleIndices()[random_i];


		return get_random_point_triangle(random_triangle);

	}


	//MOVE !!!!
	int random_int(int min, int max){
   		return rand()%(max-min+1) + min;
	}

	//MOVE !!!!
	inline glm::vec3& get_random_point_triangle(glm::uvec3 random_triangle) const
	{
	//https://math.stackexchange.com/questions/18686/uniform-random-point-in-triangle-in-3d
	glm::vec3 A = area->vertexPositions()[random_triangle[0]];
	glm::vec3 B = area->vertexPositions()[random_triangle[1]];
	glm::vec3 C = area->vertexPositions()[random_triangle[2]];
	//std::cout<<"A"<<A[0]<<", "<<A[1]<<std::endl;

	//MOVE THIS !!!!!
	std::random_device rd;
	std::default_random_engine eng(rd());
	std::uniform_real_distribution<> distr(0.f, 1.f);

	float r1 = (distr(eng));
	float r2 = (distr(eng));
	

	glm::vec3 p = (1-sqrt(r1))*A + (sqrt(r1)*(1-r2))*B + (r2*sqrt(r1))*C;
	
	return p;

	}

	bool is_area() {
		return m_is_area;
	}



	//MOEV THIS TO PRIVATE
	
	std::shared_ptr<Mesh> area;
	//change to another class

private:
	glm::vec3 m_direction;
	glm::vec3 m_color;
	float m_intensity;
	bool m_is_area;
	
	

};