#define OPENCL
#include "src/Entities.h"
#include "src/random.cl"
#include "src/intersection.cl"
#include "src/util.cl"

float3 tracePath(const Ray startRay, uint2* rngState, const float3 sunDirection, 
	global const uchar *voxels, const int3 size,
	global const Node *nodes, global const uint *ptrTable, 
	const float maxt, const int iterations, const int startMaterialType, bool showSun)
{
	const float3 skyBaseColour = (float3)(0.05f, 0.1f, 0.3f);
	const float3 sunColour = (float3)(0.63f, 0.55f, 0.5f);
	const float3 ambientColour = (float3)(0.10f, 0.16f, 0.20f);

	Ray ray = startRay;

	const float3 skyColour = skyBaseColour * clamp(-ray.dir.y, 0.1f, 1.0f);
	const float3 sunFillColour = sunColour * clamp(-dot(ray.dir, sunDirection) * 0.2f + 0.2f, 0.0f, 1.0f) * 2.5f;
	const float3 fogColour = skyColour + sunFillColour + ambientColour;

	float3 final = (float3)0.0f;
	float3 factor = (float3)1.0f; 
	int materialType = startMaterialType;

	Result result;
	for (int iter = 0; iter < iterations; iter++)
	{
		if (materialType == AIR)
			result = intersectOctree(voxels, size, nodes, ptrTable, ray, maxt * 0.8f);
		else
			result = intersectVoxels(voxels, size, ray, maxt * 0.8f, materialType);

		if (result.hit)
		{
			const float fogf = result.t * FOG_DENSITY;
			const float fog = native_recip(native_exp(fogf * fogf));

			if (fog > 0.01)
			{
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

						if (intersectVoxelsShadow(voxels, size, shadowRay, maxt))
							diffuse = 0.0f;

						colour += sunColour * 8.0f * diffuse;
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
					showSun = false;
				}

				continue;
			}
		}

		float3 sky = fogColour;
		if (showSun) sky += sunColour * 10.0f * pow(clamp(-dot(ray.dir, sunDirection), 0.0f, 1.0f) + 0.014f, 500.0f);

		if (materialType == WATER)
			final += factor * sky * (float4)(exp(-1.0f * result.t * (float3)(1.0f, 0.7f, 0.4f)), 1.0f);
		else
			final += factor * sky;

		break;
	}

	return final;
}

float3 trace(const Ray ray, uint2* rngState, const float3 sunDirection, 
	global const uchar *voxels, const int3 size,
	global const Node *nodes, global const uint *ptrTable, 
	const float maxt)
{
	const float3 skyBaseColour = (float3)(0.05f, 0.1f, 0.3f);
	const float3 sunColour = (float3)(0.63f, 0.55f, 0.5f);
	const float3 ambientColour = (float3)(0.10f, 0.16f, 0.20f);

	const Result result = intersectOctree(voxels, size, nodes, ptrTable, ray, maxt);

	const float3 skyColour = skyBaseColour * clamp(-ray.dir.y, 0.1f, 1.0f);
	const float3 sunFillColour = sunColour * clamp(-dot(ray.dir, sunDirection) * 0.2f + 0.2f, 0.0f, 1.0f) * 2.5f;
	const float3 fogColour = skyColour + sunFillColour + ambientColour;

	float3 colour = (float3)0.0f;

	if (result.hit)
	{
		const float fogf = result.t * FOG_DENSITY;
		const float fog = native_recip(native_exp(fogf * fogf));

		if (fog > 0.01)
		{
			const float3 p = ray.pos + ray.dir * result.t;

			const int3 blockPos = convert_int3(result.pos);
			const int cell = voxels[blockPos.x + (blockPos.y + blockPos.z * size.y) * size.x];
			const float3 blockColour = getCellColour(cell, p);

			colour = (float3)(0.0f, 0.0f, 0.0f);

			if (cell != 4)
			{
				float diffuse = dot(result.normal, sunDirection);

				if (diffuse > 0.0f)
				{
					float light = 1.0f;

					const int shdwSamplesU = 1;
					const int shdwSamplesV = 1;
					const float invSHDWu = 1.0f / shdwSamplesU;
					const float invSHDWv = 1.0f / shdwSamplesV;

					#pragma unroll
					for (int i = 0; i < shdwSamplesU; i++)
					{
						for (int j = 0; j < shdwSamplesV; j++)
						{
							float u1 = (i + nextRandomFloat(rngState)) * invSHDWu * SUN_SIZE;
							float u2 = (j + nextRandomFloat(rngState)) * invSHDWv;

							Ray shadowRay;
							shadowRay.pos = p - result.normal * 0.001f;
							shadowRay.dir = sampleHemisphere(u1, u2, -sunDirection);
							shadowRay.invdir = native_recip(shadowRay.dir);
							if (intersectOctreeShadow(voxels, size, nodes, ptrTable, shadowRay, maxt)) light -= invSHDWu * invSHDWv;
						}
					}

					diffuse *= light;
				}

				colour += sunColour * 8.0f * clamp(diffuse, 0.0f, 1.0f);
			}

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

				const int reflSamplesU = 1;
				const float invRFLu = 1.0f / reflSamplesU;

				const int reflSamplesV = 1;
				const float invRFLv = 1.0f / reflSamplesV;

				const float jitterU = nextRandomFloat(rngState);
				const float jitterV = nextRandomFloat(rngState);

				for (int i = 0; i < reflSamplesU; i++)
					for (int j = 0; j < reflSamplesV; j++)
					{
						float u1 = (i + jitterU) * invRFLu * 0.001f;
						float u2 = (j + jitterV) * invRFLv;

						Ray giRay;
						giRay.pos = p - result.normal * 0.001f;
						giRay.dir = sampleHemisphere(u1, u2, reflect(result.normal, ray.dir));
						giRay.invdir = native_recip(giRay.dir);

						colour += tracePath(giRay, rngState, sunDirection, voxels, size, nodes, ptrTable, maxt, 3, AIR, true) * invRFLu * invRFLv * reflectance;
					}

				const int refrSamples = 1;
				const float invRFR = 1.0f / refrSamples;

				for (int i = 0; i < refrSamples; i++)
				{
					const float cosi = dot(ray.dir, result.normal);
					const float cost2 = 1.0f - n * n * (1.0f - cosi * cosi);
					if (cost2 > 0.0f)
					{
						float3 refrDir = fast_normalize(ray.dir * n - result.normal * (cosi * n - native_sqrt(cost2)));
						float u1 = nextRandomFloat(rngState) * 0.001f;
						float u2 = (i + nextRandomFloat(rngState)) * invRFR;

						Ray tRay;
						tRay.pos = p + result.normal * 0.001f;
						tRay.dir = sampleHemisphere(u1, u2, refrDir);
						tRay.invdir = native_recip(tRay.dir);

						colour += tracePath(tRay, rngState, sunDirection, voxels, size, nodes, ptrTable, maxt, 3, WATER, true) * invRFR * (1.0f - reflectance);
					}
				}

				colour = mix(fogColour, colour, fog);
			}
			else
			{
				const int giSamplesU = 1;
				const float invGIu = 1.0f / giSamplesU;

				const int giSamplesV = 1;
				const float invGIv = 1.0f / giSamplesV;

				for (int i = 0; i < giSamplesU; i++)
					for (int j = 0; j < giSamplesV; j++)
					{
						float u1 = nextRandomFloat(rngState);
						float u2 = nextRandomFloat(rngState);

						Ray giRay;
						giRay.pos = p - result.normal * 0.001f;
						giRay.dir = -sampleHemisphere(u1, u2, result.normal);
						giRay.invdir = native_recip(giRay.dir);

						colour += tracePath(giRay, rngState, sunDirection, voxels, size, nodes, ptrTable, 64.0f, 2, AIR, false) * invGIu * invGIv;
					}

				//colour *= blockColour;
				colour = mix(fogColour, colour * blockColour, fog);
			}

			return fmax(colour, (float3)0.0f);
		}
	}

	return fogColour + sunColour * 10.0f * pow(clamp(-dot(ray.dir, sunDirection), 0.0f, 1.0f) + 0.014f, 500.0f);
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

    /*

	const float x = coords.x * inv_hheight - (float)width / height;
	const float y = coords.y * inv_hheight - 1.0f;
	const float4 screen = (float4)(x*r, y*r, 1.0f, 1.0f);

	Ray ray;
	ray.pos = *cameraPos;
	ray.dir = fast_normalize((float3)(dot(cameraViewMatrix[0], screen), dot(cameraViewMatrix[1], screen), dot(cameraViewMatrix[2], screen)));
	ray.invdir = native_recip(ray.dir);
	result = trace(ray, &rngState, sunDirection, voxels, size, nodes, ptrTable, INITIAL_MAX);
	//result = tracePath(ray, &rngState, sunDirection, voxels, size, nodes, ptrTable, INITIAL_MAX, 3, AIR);

    /*/

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

		result += trace(ray, &rngState, sunDirection, voxels, size, nodes, ptrTable, INITIAL_MAX);
	    //result += tracePath(ray, &rngState, sunDirection, voxels, size, nodes, ptrTable, INITIAL_MAX, 3, AIR, true);
	}

	result *= 0.25f;

	//*/

	write_imagef(output, coords, (float4)(result, 0.0f));
}