#ifndef SCENE_H_
#define SCENE_H_

#include <glm/glm.hpp>

using namespace glm;

struct Scene
{
	vec3	ambientColour;
	vec3	lightColour;
	vec3	lightDirection;
	float	fogDensity;
};

#endif