#ifndef _Tonemapping_h
#define _Tonemapping_h

#include <glm/glm.hpp>

using namespace glm;

// Tonemapping
vec3 convertRGBtoXYZ(vec3 RGB)
{
	return vec3(
		dot(RGB, vec3(0.5141364, 0.3238786,  0.16036376)),
		dot(RGB, vec3(0.265068,  0.67023428, 0.06409157)),
		dot(RGB, vec3(0.0241188, 0.1228178,  0.8444266))
	);
}

vec3 convertXYZtoRGB(vec3 XYZ)
{
	return vec3(
		dot(XYZ, vec3(2.5651,-1.1665,-0.3986)),
		dot(XYZ, vec3(-1.0217, 1.9777, 0.0439)),
		dot(XYZ, vec3(0.0753, -0.2543, 1.1892))
	);
}

vec3 convertXYZtoYxy(vec3 XYZ)
{
	vec3 Yxy;
	Yxy.x = XYZ.y;
	Yxy.y = XYZ.x / (XYZ.x + XYZ.y + XYZ.z);
	Yxy.z = XYZ.y / (XYZ.x + XYZ.y + XYZ.z);
	return Yxy;
}

vec3 convertYxytoXYZ(vec3 Yxy)
{
	vec3 XYZ;
	XYZ.x = Yxy.x * Yxy.y / Yxy.z;
	XYZ.y = Yxy.x;
	XYZ.z = Yxy.x * (1.0f - Yxy.y - Yxy.z) / Yxy.z;
	return XYZ;
}

vec3 tonemap(vec3 input)
{
	//return input * 0.2f;

	/*float3 x = max(0.0f, input - 0.004f) * 0.05f;
	float3 output = (x*(6.2f*x + 0.5f)) / (x*(6.2f*x + 1.7f) + 0.06f);
	return output;*/

	vec3 Yxy = convertXYZtoYxy(convertRGBtoXYZ(input));

	float x = glm::max(0.0f, Yxy.x - 0.004f) * 0.1f;
	Yxy.x = (x*(6.2f*x + 0.5f)) / (x*(6.2f*x + 1.7f) + 0.06f);

	return convertXYZtoRGB(convertYxytoXYZ(Yxy));
}

#endif