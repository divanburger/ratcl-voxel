#ifndef _Noise_h_
#define _Noise_h_

#include <cmath>
#include <glm/glm.hpp>

#include "Math.h"

using namespace glm;

float	noise(int x, int seed);
float	noise(int x, int y, int seed);
float	noise(int x, int y, int z, int seed);

float	smoothNoise(int x, int seed);
float	smoothNoise(int x, int y, int seed);
float	smoothNoise(int x, int y, int z, int seed);

float	interpolatedNoise(float x, int seed);
float	interpolatedNoise(float x, float y, int seed);
float	interpolatedNoise(float x, float y, float z, int seed);

float	perlinNoise(float persistence, int octaves, float x, int seed);
float	perlinNoise(float persistence, int octaves, float x, float y, int seed);
float	perlinNoise(float persistence, int octaves, float x, float y, float z, int seed);

inline float perlinNoise(float persistence, int octaves, vec2 p, int seed) 
{
	return perlinNoise(persistence, octaves, p.x, p.y, seed);
}

inline float perlinNoise(float persistence, int octaves, vec3 p, int seed)
{
	return perlinNoise(persistence, octaves, p.x, p.y, p.z, seed);
}

#endif
