#define OPENCL
#include "src/Entities.h"
#include "src/random.cl"
#include "src/intersection.cl"
#include "src/util.cl"

float3 tracePath(const Ray startRay, uint2* rngState, const float3 sunDirection, 
	global const uchar *voxels, const int3 size,
	global const Node *nodes, global const uint *ptrTable, 
	const float maxt, const int iterations, const int startMaterialType)
{
	const float3 skyBaseColour = (float3)(0.05f, 0.1f, 0.3f);
	const float3 sunColour = (float3)(0.63f, 0.55f, 0.5f);
	const float3 ambientColour = (float3)(0.10f, 0.16f, 0.20f);

	Ray ray = startRay;

	const float3 skyColour = skyBaseColour * clamp(-ray.dir.y, 0.2f, 1.0f);
	const float3 sunFillColour = sunColour * clamp(-dot(ray.dir, sunDirection) * 0.2f + 0.1f, 0.0f, 1.0f) * 2.5f;
	const float3 fogColour = skyColour + sunFillColour + ambientColour;

	float3 final = (float3)0.0f;
	float3 factor = (float3)1.0f; 
	int materialType = startMaterialType;

	float tmax = maxt;

	Result result;
	for (int iter = 0; iter < iterations; iter++)
	{
		if (materialType == AIR && maxt > 50.0f)
			result = intersectOctree(voxels, size, nodes, ptrTable, ray, tmax);
		else
			result = intersectVoxels(voxels, size, ray, tmax, materialType);

		const float fogf = result.t * FOG_DENSITY;
		const float fog = native_recip(native_exp(fogf * fogf));

		if (!result.hit || fog < 0.01)
		{
			if (materialType == WATER) factor *= (float3)(exp(-1.0f * result.t * (float3)(1.0f, 0.7f, 0.4f)));

			final += factor * fogColour;

			if (iter == 0)
			 	final += factor * sunColour * 10.0f * pow(clamp(-dot(ray.dir, sunDirection), 0.0f, 1.0f) + 0.014f, 500.0f);

			break;
		}

		const float3 p = ray.pos + ray.dir * result.t;

		const int3 blockPos = convert_int3(result.pos);
		const int cell = voxels[blockPos.x + (blockPos.y + blockPos.z * size.y) * size.x];
		const float3 blockColour = getCellColourSimple(cell);

		float3 colour = (float3)(0.0f, 0.0f, 0.0f);

		if (cell != 4)
		{
			float diffuse = dot(result.normal, sunDirection);

			if (diffuse > 0.0f)
			{
				float u1 = nextRandomFloat(rngState) * SUN_SIZE;
				float u2 = nextRandomFloat(rngState);

				Ray shadowRay;
				shadowRay.pos = p - result.normal * 0.001f;
				shadowRay.dir = sampleHemisphere(u1, u2, -sunDirection);
				shadowRay.invdir = native_recip(shadowRay.dir);

				if (intersectVoxelsShadow(voxels, size, shadowRay, tmax))
					diffuse = 0.0f;

				colour += sunColour * 10.0f * diffuse;
			}
		}

		if (materialType == WATER) // Apply Beer's law
			factor *= exp(-1.0f * result.t * (float3)(1.0f, 0.7f, 0.4f));
		else
			factor *= fog;

		final += factor * colour * blockColour;

		if (iter + 1 == iterations)	
			break;

		if (cell == 4)
		{
			const float airDensity = 1.0f;
			const float waterDensity = 1.33f;
			const float n = airDensity / waterDensity;
			const float k = 1.0f - dot(ray.dir, result.normal);

			const float a = waterDensity - airDensity;
			const float b = waterDensity + airDensity;
			const float refl0 = (a*a)/(b*b);

			const float reflectance = refl0 + (1.0f - refl0) * k * k * k * k;

			const float choice = nextRandomFloat(rngState);
			const float probability = clamp(reflectance * 0.8f + 0.1f, 0.0f, 1.0f);

			if (choice < probability)
			{
				float u1 = nextRandomFloat(rngState) * 0.0025f;
				float u2 = nextRandomFloat(rngState);

				ray.pos = p - result.normal * 0.001f;
				ray.dir = sampleHemisphere(u1, u2, reflect(result.normal, ray.dir));
				ray.invdir = native_recip(ray.dir);

				factor *= reflectance / probability;
			}
			else
			{
				const float cosi = dot(ray.dir, result.normal);
				const float cost2 = 1.0f - n * n * (1.0f - cosi * cosi);

				if (cost2 > 0.0f)
				{
					float3 refrDir = fast_normalize(ray.dir * n - result.normal * (cosi * n - native_sqrt(cost2)));
					float u1 = nextRandomFloat(rngState) * 0.0025f;
					float u2 = nextRandomFloat(rngState);

					ray.pos = p + result.normal * 0.001f;
					ray.dir = sampleHemisphere(u1, u2, refrDir);
					ray.invdir = native_recip(ray.dir);

					factor *= (1.0f - reflectance) / (1.0f - probability);
					materialType = WATER;
				}
				else
					break;
			}
		}
		else
		{
			float u1 = nextRandomFloat(rngState);
			float u2 = nextRandomFloat(rngState);

			ray.pos = p - result.normal * 0.001f;
			ray.dir = -sampleHemisphere(u1, u2, result.normal);
			ray.invdir = native_recip(ray.dir);

			factor *= blockColour;
		}

		tmax = tmax - result.t;
	}

	return final;
}

kernel void raytrace(const int initial_seed, write_only image2d_t output,
	global const float3 *cameraPos, global const float4 *cameraViewMatrix,
	global const float3 *sunDir, global const uchar *voxels, global const int3 *voxelsSize, 
	global const Node *nodes, global const uint *ptrTable)
{
	const int width = get_image_width(output);
	const int height = get_image_height(output);

	const int2 coords = (int2)(get_global_id(0), get_global_id(1));

	if (coords.x > width) return;
	if (coords.y > height) return;

	const uint seed = initial_seed + coords.x * 6971 + coords.y * 31337;

	uint2 rngState = {noise1i(0, seed), initial_seed};

	const int3 size = *voxelsSize;
	const float3 sunDirection = *sunDir;
	const float fov = 75.0f;

	const float r = native_tan(radians(fov * 0.5f));
	const float2 screenCoord = (float2)(coords.x, coords.y);
	const float inv_hheight = 2.0f / height;

	int p = 1;

	float3 result = (float3)(0.0f, 0.0f, 0.0f);

	float2 locations[4] = {(float2)(0.0f, 0.5f), (float2)(0.25f, 0.0f), (float2)(0.5f, 0.75f), (float2)(0.75f, 0.25f)};

	for (int s = 0; s < 4; s++)
	{
		const float x = (coords.x+locations[s].x) * inv_hheight - (float)width / height;
		const float y = (coords.y+locations[s].y) * inv_hheight - 1.0f;
		const float4 screen = (float4)(x*r, y*r, 1.0f, 1.0f);

		Ray ray;
		ray.pos = *cameraPos;
		ray.dir = fast_normalize((float3)(dot(cameraViewMatrix[0], screen), dot(cameraViewMatrix[1], screen), dot(cameraViewMatrix[2], screen)));
		ray.invdir = native_recip(ray.dir);

		result += tracePath(ray, &rngState, sunDirection, voxels, size, nodes, ptrTable, INITIAL_MAX, 3, AIR);
	}

	result *= 0.25f;

	write_imagef(output, coords, (float4)(result, 0.0f));
}