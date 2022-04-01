// ----------------------------------------------
// Polytechnique - INF584 "Image Synthesis"
//
// Base code for practical assignments.
//
// Copyright (C) 2022 Tamy Boubekeur
// All rights reserved.
// ----------------------------------------------
#include "RayTracer.h"

#include<chrono>

#include "Console.h"
#include "Camera.h"
#include "BoundingBox.h"
#include "PBR.h"
#include "BVH.h"

//--------------------------------
#include <random>
std::random_device rd;
std::default_random_engine eng(rd());
std::uniform_real_distribution<> distr(0.f, 1.f);
size_t max_bounces = 2;
//----------------------------------

RayTracer::RayTracer():
	Renderer (), 
	m_image_multichannelPtr(std::make_shared<Image_multichannel>()) {}

RayTracer::~RayTracer() {}

void RayTracer::init (const std::shared_ptr<Scene> scenePtr) {
	buildBVH (scenePtr);

	// m_LPEs.full  = (std::make_shared<Image>());
	// m_LPEs.ld12e  = (std::make_shared<Image>());
	// m_LPEs.lde  = (std::make_shared<Image>());
	// m_LPEs.lse  = (std::make_shared<Image>());
	// m_LPEs.indirect = (std::make_shared<Image>());

}

void RayTracer::buildBVH (const std::shared_ptr<Scene> scene) {
	if (m_bvh != nullptr)
		m_bvh.reset();
	m_bvh = std::make_unique<BVH>(scene);
	std::cout<<"built BVH"<<std::endl;
}

void RayTracer::render (const std::shared_ptr<Scene> scenePtr) {

	size_t width = m_image_multichannelPtr->width();
	size_t height = m_image_multichannelPtr->height();
	m_image_multichannelPtr->clear (scenePtr->backgroundColor ());

	const auto cameraPtr = scenePtr->camera();
	int sample_count =  10;
	
#pragma omp parallel for collapse(3)
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {//loop over all pixels
			glm::vec3 colorResponse (0.f, 0.f, 0.f);

			//write this better !
			ColorResponses color_responses;
			color_responses.full = glm::vec3(0.f, 0.f, 0.f);
			color_responses.ld12e = glm::vec3(0.f, 0.f, 0.f);
			color_responses.lde = glm::vec3(0.f, 0.f, 0.f);
			color_responses.lse = glm::vec3(0.f, 0.f, 0.f);
			color_responses.indirect = glm::vec3(0.f, 0.f, 0.f);
//#pragma omp parallel for
			for (int i =0; i< sample_count; i++){//N samples for monte carlo path tracing

						float r1 = (distr(eng)-0.5)*0.6;
						float r2 = (distr(eng)-0.5)*0.6;
						
						//shoot a random ray from X
						Ray ray = cameraPtr->rayAt((float(x) + 0.5 + r1 )/ width , 1.f - (float(y) + 0.5+ r2) / height); 
						colorResponse += sample (color_responses, scenePtr, ray, 0, 0, max_bounces); 
			}
			colorResponse/=sample_count;

			//write this better !
			color_responses.full /= sample_count;
			color_responses.ld12e /= sample_count;
			color_responses.lde /= sample_count;
			color_responses.lse /= sample_count;
			color_responses.indirect /= sample_count;

			//TO DO: use enum type instead
			m_image_multichannelPtr->operator()(0, x, y) = color_responses.full;
			m_image_multichannelPtr->operator()(1, x, y) = color_responses.ld12e;
			m_image_multichannelPtr->operator()(2, x, y) = color_responses.lde;
			m_image_multichannelPtr->operator()(3, x, y) = color_responses.lse;
			m_image_multichannelPtr->operator()(4, x, y) = color_responses.indirect;
		}
	}

}

//pathtrace
bool RayTracer::rayTrace(const Ray& ray, const std::shared_ptr<Scene> scene, size_t originMeshIndex, size_t originTriangleIndex, Hit& hit, bool anyHit) const {
	float closest = std::numeric_limits<float>::max();
	bool intersectionFound = false;
	std::vector<std::pair<size_t,size_t>> candidateMeshTrianglePairs;

	//BVH ray intersection
	candidateMeshTrianglePairs.reserve (64);
	m_bvh->intersect(ray, candidateMeshTrianglePairs);

	//Naive ray-intersection routine
	// for (size_t mIndex = 0; mIndex < scene->numOfMeshes (); mIndex++) {
	// 	const auto& T = scene->mesh(mIndex)->triangleIndices();
	// 	for (size_t tIndex = 0; tIndex < T.size(); tIndex++) 
	// 		candidateMeshTrianglePairs.push_back (std::pair<size_t,size_t> (mIndex, tIndex));
	// }

	for (size_t i = 0; i < candidateMeshTrianglePairs.size(); i++) {
		size_t mIndex = candidateMeshTrianglePairs[i].first;
		size_t tIndex = candidateMeshTrianglePairs[i].second;
		if (anyHit && mIndex == originMeshIndex && tIndex == originTriangleIndex)
				continue;
		const auto& mesh = scene->mesh(mIndex);
		const auto& triangleIndices = mesh->triangleIndices();
		const glm::uvec3& triangle = triangleIndices[tIndex];
		const auto& P = mesh->vertexPositions();
		glm::mat4 modelMatrix = mesh->computeTransformMatrix ();
		float ut, vt, dt;
		if (ray.triangleIntersect(glm::vec3 (modelMatrix * glm::vec4 (P[triangle[0]], 1.0)), 
								  glm::vec3 (modelMatrix * glm::vec4 (P[triangle[1]], 1.0)), 
								  glm::vec3 (modelMatrix * glm::vec4 (P[triangle[2]], 1.0)), 
								  ut, vt, dt) == true) {
			if (dt > 0.f) {
				if (anyHit)
					return true;
				else if (dt < closest) {
					intersectionFound = true;
					closest = dt;
					hit.m_meshIndex = mIndex;
					hit.m_triangleIndex = tIndex;
					hit.m_uCoord = ut;
					hit.m_vCoord = vt;
					hit.m_distance = dt;
				}
			}
		}
	}
	return intersectionFound;
}

glm::vec3 RayTracer::lightRadiance (const std::shared_ptr<LightSource> lightPtr, const glm::vec3 & position) const {
	return lightPtr->color() * lightPtr->intensity() * glm::pi<float>(); 
}


pair<glm::vec3, glm::vec3> RayTracer::materialReflectance (const std::shared_ptr<Scene> scenePtr,
										  const std::shared_ptr<Material> materialPtr, 
										  const glm::vec3& wi, 
										  const glm::vec3& wo, 
										  const glm::vec3& n) const {

		return  BRDF (wi, wo, n, materialPtr->albedo (), materialPtr->roughness (), materialPtr->metallicness ()); 
}

glm::vec3 RayTracer::shade(ColorResponses& color_responses, const std::shared_ptr<Scene> scenePtr, const Ray & ray, const Hit& hit, size_t bounces) {

	const auto& mesh = scenePtr->mesh(hit.m_meshIndex);
	const std::shared_ptr<Material> materialPtr = scenePtr->material(scenePtr->mesh2material(hit.m_meshIndex));
	const auto& P = mesh->vertexPositions();
	const auto& N = mesh->vertexNormals();
	glm::mat4 modelMatrix = mesh->computeTransformMatrix ();
	const glm::uvec3 & triangle = mesh->triangleIndices()[hit.m_triangleIndex];
	float w = 1.f - hit.m_uCoord - hit.m_vCoord;
	glm::vec3 hitPosition = barycentricInterpolation(P[triangle[0]], P[triangle[1]], P[triangle[2]], w, hit.m_uCoord, hit.m_vCoord);
	hitPosition = glm::vec3 (modelMatrix * glm::vec4 (hitPosition, 1.0));
	glm::vec3 unormalizedHitNormal = barycentricInterpolation(N[triangle[0]], N[triangle[1]], N[triangle[2]], w, hit.m_uCoord, hit.m_vCoord);
	glm::mat4 normalMatrix = glm::transpose (glm::inverse (modelMatrix));
	glm::vec3 hitNormal = normalize (glm::vec3 (normalMatrix * glm::vec4 (normalize (unormalizedHitNormal), 1.0)));

	glm::vec3 wo = normalize(-ray.direction ());
	
	glm::vec3 colorResponse (0.f, 0.f, 0.f);

	std::pair<glm::vec3, glm::vec3> fs_fd;

	//DIRECT LIGHTING------------------------------------
	for (size_t i = 0; i < scenePtr->numOfLightSources(); ++i) {
		const std::shared_ptr<LightSource> light = scenePtr->lightSource(i);
			int Ni;
			glm::vec3 wi;

			if (light->is_area()) {
				Ni = 1;
			}
			else {
				//directional light (hard shadows only)
				Ni = 1;
				wi = normalize(-light->direction());
			}

			//monte carlo 
			//#pragma omp parallel for
			for (int j = 0; j < Ni; j++) {

				if (light->is_area()) {
					int s = light->area->triangleIndices().size();
					int random_i = rand() % s;
					glm::uvec3 random_triangle = light->area->triangleIndices()[random_i];
					glm::vec3 A = light->area->vertexPositions()[random_triangle[0]];
					glm::vec3 B = light->area->vertexPositions()[random_triangle[1]];
					glm::vec3 C = light->area->vertexPositions()[random_triangle[2]];
					float r1 = (distr(eng));
					float r2 = (distr(eng));
					glm::vec3 sample_light = (1 - sqrt(r1)) * A + (sqrt(r1) * (1 - r2)) * B + (r2 * sqrt(r1)) * C;
					//std::cout<<"sample = "<<sample_light[0]<<", "<<sample_light[1]<<", "<<sample_light[2]<<std::endl;
					wi = normalize(sample_light - hitPosition);
				}

				float wiDotN = max(0.f, dot(wi, hitNormal));

				if (wiDotN <= 0.f)
					continue;

				Ray to_light = Ray(hitPosition, wi);
				Hit hit2;
				bool shadow = rayTrace(to_light, scenePtr, 0, 0, hit2, true);

				if (!shadow) { 

						fs_fd = materialReflectance(scenePtr, materialPtr, wi, wo, hitNormal);
						glm::vec3 f = fs_fd.first + fs_fd.second;
						colorResponse += (lightRadiance(light, hitPosition) * f * wiDotN) / float(Ni);

						color_responses.full += (lightRadiance(light, hitPosition) * f * wiDotN) / float(Ni);

						if (bounces == max_bounces) {//first bounce
							color_responses.lde += (lightRadiance(light, hitPosition) * fs_fd.first * wiDotN) / float(Ni);
							color_responses.lse += (lightRadiance(light, hitPosition) * fs_fd.second * wiDotN) / float(Ni);
						}

						if (bounces >= max_bounces - 1) {//two diffuse bounces
							color_responses.ld12e += (lightRadiance(light, hitPosition) * fs_fd.first * wiDotN) / float(Ni);
						}
						
						if ((bounces < max_bounces) ) {//diffuse indirect illumination
							color_responses.indirect += (lightRadiance(light, hitPosition) * fs_fd.first * wiDotN) / (float(Ni)* 1 / (2 * PI));
						}
				}
			}

		//russian roulette
		if (bounces > 1) {
			float rr = distr(eng);
			if (rr > 0.5) {
				bounces = 0;
				return colorResponse;
			}
		}

		//INDIRECT LIGHTING----------------------------------
		//https://www.scratchapixel.com/lessons/3d-basic-rendering/global-illumination-path-tracing/global-illumination-path-tracing-practical-implementation

		if (bounces>1){

			float r1 = (distr(eng));
			float r2 = (distr(eng));
			glm::vec3 random_point = uniformSampleHemisphere(r1, r2);
			//transform to the local coordinate system

			Ray random_ray = Ray(hitPosition, random_point);
			float pdf = 1 / (2 * PI);
			fs_fd = materialReflectance(scenePtr, materialPtr, random_point, wo, hitNormal);
			glm::vec3 color = sample(color_responses, scenePtr, random_ray, 0, 0, bounces - 1) * (fs_fd.first + fs_fd.second) / pdf;
			color_responses.full += color;
		}
		return colorResponse;
	}
	return colorResponse;
}


glm::vec3 RayTracer::sample (ColorResponses& color_responses, const std::shared_ptr<Scene> scenePtr, const Ray & ray, size_t originMeshIndex, size_t originTriangleIndex, size_t bounces) {
	
	Hit hit;
	bool intersectionFound = rayTrace(ray, scenePtr, originMeshIndex, originTriangleIndex, hit, false);
	if (intersectionFound && hit.m_distance > 0.f && bounces >0) {
		return shade(color_responses, scenePtr, ray, hit, bounces);
	} else 
		return scenePtr->backgroundColor (); 
}

//move
glm::vec3 RayTracer::uniformSampleHemisphere(const float& r1, const float& r2)
{
	float sinTheta = sqrtf(1 - r1 * r1);
	float phi = 2 * PI * r2;
	float x = sinTheta * cosf(phi);
	float z = sinTheta * sinf(phi);
	return glm::vec3(x, r1, z); //wtf is u1 ???
}




