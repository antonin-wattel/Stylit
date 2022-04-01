// ----------------------------------------------
// Polytechnique - INF584 "Image Synthesis"
//
// Base code for practical assignments.
//
// Copyright (C) 2022 Tamy Boubekeur
// All rights reserved.
// ----------------------------------------------
#pragma once

#include <cmath>
#include <algorithm>
#include <limits>
#include <memory>
#include <glm/glm.hpp>
#include <glm/ext.hpp>

#include "Renderer.h"
#include "Scene.h"
#include "BVH.h" 
#include "Image_multichannel.h"

using namespace std;

class RayTracer : public Renderer {
public:
	
	RayTracer();
	virtual ~RayTracer();

	inline void setResolution (int width, int height) { 
		m_image_multichannelPtr = make_shared<Image_multichannel>(width, height, 5);
	}

	inline std::shared_ptr<Image_multichannel> image_multichannel() { return m_image_multichannelPtr; }
	inline const std::shared_ptr<Image_multichannel> image_multichannel() const { return m_image_multichannelPtr; }

	void init (const std::shared_ptr<Scene> scenePtr);
	virtual void render (const std::shared_ptr<Scene> scenePtr) final;

private:
	template<typename T>
	inline T barycentricInterpolation (const T & p0, const T & p1, const T & p2, float w, float u, float v) const {	return w * p0 + u * p1 + v * p2; }

	struct Hit {
		size_t m_meshIndex;
		size_t m_triangleIndex;
		float m_uCoord;
		float m_vCoord;
		float m_distance;
	};

	//to do: use enum type instead
	struct ColorResponses{
		glm::vec3 full;//full global illumination
		glm::vec3 lde;//direct diffuse
		glm::vec3 lse;//direct specular
		glm::vec3 ld12e;//first two diffuse bounces
		glm::vec3 indirect;
	};


	bool rayTrace(const Ray& ray, const std::shared_ptr<Scene> scene, size_t originMeshIndex, size_t originTriangleIndex, Hit& hit, bool anyHit) const;
	bool pathTrace(const Ray& ray, const std::shared_ptr<Scene> scene, size_t originMeshIndex, size_t originTriangleIndex, Hit& hit, bool anyHit) const;
	inline bool rayTrace(const Ray& ray, const std::shared_ptr<Scene> scene, size_t originMeshIndex, size_t originTriangleIndex) { return rayTrace (ray, scene, originMeshIndex, originTriangleIndex, Hit(), true);}
	glm::vec3 lightRadiance (const std::shared_ptr<LightSource> lightPtr, const glm::vec3 & position) const;
	
	pair<glm::vec3, glm::vec3> materialReflectance (const std::shared_ptr<Scene> scenePtr, 
								   const std::shared_ptr<Material> material, 
								   const glm::vec3& wi, 
								   const glm::vec3& wo, 
								   const glm::vec3& n) const;

	glm::vec3 shade(ColorResponses& color_responses, const std::shared_ptr<Scene> scenePtr, const Ray & ray, const Hit& hit, size_t bounces);
	glm::vec3 sample (ColorResponses& color_responses, const std::shared_ptr<Scene> scenePtr, const Ray & ray, size_t originMeshIndex, size_t originTriangleIndex, size_t bounces);
	void RayTracer::buildBVH(const std::shared_ptr<Scene> scene);
	

	//move
	glm::vec3 uniformSampleHemisphere(const float& r1, const float& r2);

	std::shared_ptr<Image_multichannel> m_image_multichannelPtr;
	std::unique_ptr<BVH> m_bvh;
};



