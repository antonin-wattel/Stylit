#pragma once

#include <glm/glm.hpp>
#include <glm/ext.hpp>

const float PI = 3.14159265358979323846;

inline float sqr (float x) { return x*x; }

inline float GGX (float NdotH, float roughness) {
	if (roughness >= 1.0f) 
		return glm::one_over_pi<float>();
	float alpha = sqr (roughness);
	float tmp = alpha / std::max(1e-8f,(NdotH*NdotH*(sqr (alpha)-1.0f)+1.0f));
	return sqr (tmp) * glm::one_over_pi<float>();
}

inline glm::vec3 SchlickSGFresnel (float VdotH, glm::vec3 F0) {
	float sphg = exp2 ((-5.55473f*VdotH - 6.98316f) * VdotH);
	return F0 + (glm::vec3(1.0f) - F0) * sphg;
}

inline float smithG_GGX (float NdotV, float alphaG) {
	return 2.0f/(1.0f + sqrt (1.0f + sqr (alphaG) * (1.0f - sqr (NdotV) / sqr(NdotV))));
}

inline float G1 (float D, float k) {
	return 1.0f / (D * (1.0f-k) + k);
}

inline float geometry (float NdotL, float NdotV, float roughness) {
	float k = roughness * roughness * 0.5f;
	return G1(NdotL,k) * G1(NdotV,k);
}


inline pair<glm::vec3, glm::vec3> BRDF (glm::vec3 L, glm::vec3 V, glm::vec3 N,  glm::vec3 albedo, float roughness, float metallic)  {
//inline glm::vec3 BRDF (glm::vec3 L, glm::vec3 V, glm::vec3 N,  glm::vec3 albedo, float roughness, float metallic)  {
	glm::vec3 diffuseColor = albedo * (1.0f - metallic);
	glm::vec3 specularColor = mix(glm::vec3(0.08f), albedo, metallic);

	float NdotL = std::max (0.0f, dot (N, L));
	float NdotV = std::max (0.0f, dot (N, V));

	if (NdotL <= 0.0f){
		glm::vec3 z(0.0f);
		return pair(z, z); 
	}
	

	glm::vec3 H = normalize (L + V);
	float NdotH = std::max (0.0f, dot (N, H));
	float VdotH = std::max (0.0f, dot (V, H));

	float D = GGX (NdotH, roughness);
	glm::vec3  F = SchlickSGFresnel (VdotH, specularColor);
	float G = geometry (NdotL, NdotV, roughness);

	glm::vec3 fd = diffuseColor * (glm::vec3(1.0f)-specularColor) / glm::pi<float>();
	glm::vec3 fs = F * D * G / (4.0f);
	//std::cout<<"fs: "<<fs[0]<<std::endl;

	//return fs;

	return 	pair(fd, fs);
	//return (fd+fs);
	//return (fd + fs); //speparate these into different buffers for LPE!
	//instead: return [(fd+fs, fd, fs)]
}


//importance sampling for the microfacet GGX BRDF
//https://hal.archives-ouvertes.fr/hal-01509746/document
inline glm::vec3 sampleGGXVNDF(glm::vec3 wo, float alpha_x, float alpha_y, float U1, float U2)
{
	// stretch view
	wo = normalize(glm::vec3(alpha_x * wo.x, alpha_y * wo.y, wo.z));
	// orthonormal basis
	glm::vec3 T1 = (wo.z < 0.9999) ? normalize(glm::cross(wo, glm::vec3(0,0,1))) : glm::vec3(1,0,0);
	glm::vec3 T2 = cross(T1, wo);
	// sample point with polar coordinates (r, phi)
	float a = 1.0 / (1.0 + wo.z);
	float r = sqrt(U1);
	float phi = (U2<a) ? U2/a * PI : PI + (U2-a)/(1.0-a) * PI; //TO DO: USE glm::pi<float> INSTEAD
	float P1 = r*cos(phi);
	float P2 = r*sin(phi)*((U2<a) ? 1.0 : wo.z);
	// compute normal
	float tmp = (sqrt(max(0.0, 1.0 - P1*P1 - P2*P2)));
	glm::vec3 N =  P1*T1 + P2*T2 + tmp*wo;
	// unstretch
	N = normalize(glm::vec3(alpha_x*N.x, alpha_y*N.y, max(0.f, N.z)));

	//glm::vec3 wi = wo.reflect(N); 
	return N;
}


void ImportanceSampleGgxVdn(glm::vec3  wo, glm::vec3 wm, glm::vec3 albedo, float roughness, float metallic,
                            glm::vec3 & wi, float r0, float r1){

	glm::vec3 diffuseColor = albedo * (1.0f - metallic);
	glm::vec3 specularColor = mix(glm::vec3(0.08f), albedo, metallic);

    float a = roughness; //change !
    float a2 = a * a;

    wm = sampleGGXVNDF(wo, roughness, roughness, r0, r1); // ?

    wi = 2.0f * dot(wo, wm) * wm - wo; // REWRITE !


	//change the reflectance !!!!!
    // if(BsdfNDot(wi) > 0.0f) {

    //     float3 F = SchlickFresnel(specularColor, Dot(wi, wm));
    //     float G1 = SmithGGXMasking(wi, wo, a2);
    //     float G2 = SmithGGXMaskingShadowing(wi, wo, a2);

    //     reflectance = F * (G2 / G1);
        
    // }
    // else {
    //     reflectance = float3::Zero_;
    // }
}


//https://computergraphics.stackexchange.com/questions/7656/importance-sampling-microfacet-ggx#7662
void BRDF_sample(glm::vec3 wo, float r0, float r1, float roughness, glm::vec3& wi, float& pdf) {
   
    float a = roughness*roughness;
	float a2 = a*a;
    float theta = acos(sqrt((1 - r0) / ((a2-1)*r0 + 1)));
    float phi = 2 * PI * r1;
    float x = sin(theta) * cos(phi);
    float y = cos(theta);
    float z = sin(theta) * sin(phi);
    glm::vec3 wm = {x, y, z};
	normalize(wm);
   	wi =  2.0f * dot(wo, wm) * wm - wo; // reflect 
	
 	float cosTheta = wm[2];
    float exp = (a2-1)*cosTheta*cosTheta + 1;
    float D = a2 / (PI * exp * exp);
    pdf = (D *cosTheta) / (4 * dot(wo, wm));
}


