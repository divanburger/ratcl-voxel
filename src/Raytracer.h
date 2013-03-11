#ifndef RAYTRACER_H_
#define RAYTRACER_H_

#include <glm/glm.hpp>

#include "Scene.h"
#include "Ray.h"
#include "VoxelGrid.h"
#include "Noise.h"
#include "Camera.h"
#include "Tonemapping.h"

using namespace glm;

class Raytracer
{
	public:
		bool debug;

		Raytracer(int width, int height, Scene& scene, VoxelGrid& voxels, Camera& camera) : width(width), height(height), scene(scene), voxels(voxels), camera(camera)
		{

		}

		vec4 raytrace(Ray ray, int seed, int depth, float max)
		{
			const vec3 sunColour = scene.lightColour * 2.5f;
			const vec3 SKY_BASE_COLOUR = vec3(0.05f, 0.1f, 0.3f);
			vec3 skyColour = SKY_BASE_COLOUR * clamp(-ray.dir.y, 0.1f, 1.0f);
			vec3 sunFillColour = sunColour * clamp(-dot(ray.dir, scene.lightDirection) * 0.2f + 0.2f, 0.0f, 1.0f);

			vec3 fogColour = skyColour + sunFillColour + scene.ambientColour;
			vec3 colour = fogColour;

			VoxelGrid::Result result;

			//if (debug)
			//	result = voxels.intersectOctreeFast(ray, max);
			//else
				result = voxels.intersectOctree(ray, max);

			if (result.hit)
			{
				vec3 p = ray.pos + ray.dir * result.t;

				float fog = result.t * scene.fogDensity;
				fog = 1.0f / pow((float)M_E, fog * fog);

				if (fog > 0.01)
				{
					ivec3 blockPos = ivec3(result.pos);
					int cell = voxels.get(blockPos.x, blockPos.y, blockPos.z);
					vec3 blockColour;

					if (cell == 1)
						blockColour = vec3(0.65f, 0.9f, 0.2f);
					else if (cell == 2)
						blockColour = vec3(0.4f, 0.2f, 0.05f);
					else if (cell == 3)
						blockColour = vec3(0.4f, 0.4f, 0.4f);
					else if (cell == 4)
						blockColour = vec3(0.2f, 0.35f, 0.6f);
					else if (cell == 5)
						blockColour = vec3(0.25f, 0.13f, 0.03f);
					else if (cell == 6)
						blockColour = vec3(0.4f, 0.8f, 0.05f);
					else
						blockColour = vec3(0.1f, 0.15f, 0.2f);

					blockColour *= 0.4f;

					vec3 colour;

					// Diffuse Factor
					float diffuse = dot(scene.lightDirection, normalize(result.normal));

					// Shadowing
					if (diffuse > 0.0f)
					{
						Ray shadowRay = Ray(p - result.normal * 0.001f, -scene.lightDirection);
						//if (voxels.intersectShadow(shadowRay, max)) d = 0.0f;
						if (voxels.intersectShadowOctree(shadowRay, max)) diffuse = 0.0f;
					}

					colour += sunColour * clamp(diffuse, 0.0f, 1.0f) * 0.6f + scene.ambientColour;

					return vec4(tonemap(mix(fogColour, blockColour * colour, fog)), result.t);
				}
			}

			if (depth == 0)
			{
				vec3 sunColour = scene.lightColour * pow(clamp(-dot(ray.dir, scene.lightDirection), 0.0f, 1.0f) + 0.014f, 500.0f);
				colour += sunColour;
			}

			return vec4(tonemap(colour), INFINITY);
		}

		vec4 trace(float x, float y, int seed)
		{
			debug = (x == 128 && y == 128);

			const float fov = 75.0f;
			const float aspect = (float)width / height;
			const float fx = x * 2.0f / height - aspect;
			const float fy = y * 2.0f / height - 1.0f;
			const float r = tan(fov / 180.0f * M_PI * 0.5f);

			const vec3 dir = vec3(camera.getViewMatrix() * vec4(fx*r, fy*r, 1.0f, 1.0f));
			Ray ray(camera.position, normalize(dir));
			
			return raytrace(ray, seed, 0, 1000.0f);
		}

	private:
		int		width;
		int		height;

		Scene&		scene;
		VoxelGrid&	voxels;
		Camera&		camera;
};

#endif