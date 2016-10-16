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
	const float3 ambientColour = (float3)(0.10f, 0.16f, 0.20f);

	Ray ray = startRay;

	const float3 skyColour = SKY_BASE_COLOUR * clamp(-ray.dir.y, 0.1f, 1.0f);
	const float3 sunFillColour = SUN_COLOUR * clamp(-dot(ray.dir, sunDirection) * 0.2f + 0.2f, 0.0f, 1.0f) * 2.5f;
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
						Ray shadowRay;
						shadowRay.pos = p - sunDirection * 0.001f;
						shadowRay.dir = -sunDirection;
						shadowRay.invdir = native_recip(shadowRay.dir);

						if (intersectVoxelsShadow(voxels, size, shadowRay, INITIAL_MAX))
							diffuse = 0.0f;
					}

					colour += SUN_COLOUR * 2.5f * clamp(diffuse, 0.0f, 1.0f);

					float ambientFactor = 1.0f;

					const int ambientOcclusionSamples = 1;
					const float invAO = 1.0f / ambientOcclusionSamples;

					#pragma unroll
					for (int i = 0; i < ambientOcclusionSamples; i++)
					{
						float u1 = nextRandomFloat(rngState);
						float u2 = (i + nextRandomFloat(rngState)) * invAO;

						Ray aoRay;
						aoRay.pos = p - result.normal * 0.001f;
						aoRay.dir = -sampleHemisphere(u1, u2, result.normal);
						aoRay.invdir = native_recip(aoRay.dir);

						float u3 = nextRandomFloat(rngState) * 6.5f + 0.5f;

						if (intersectVoxelsShadow(voxels, size, aoRay, u3)) ambientFactor -= invAO;
					}

					colour += ambientColour * ambientFactor;
				}

				if (materialType == WATER) // Apply Beer's law
					factor *= exp(-1.0f * result.t * (float3)(1.0f, 0.7f, 0.4f));

				final += factor * mix(fogColour, colour * blockColour, fog);

				if (materialType != WATER)
					factor *= fog;

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
						ray.pos = p - result.normal * 0.001f;
						ray.dir = reflect(result.normal, ray.dir);
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

							ray.pos = p + result.normal * 0.001f;
							ray.dir = refrDir;
							ray.invdir = native_recip(ray.dir);

							factor *= (1.0f - reflectance) / (1.0f - probability);
							materialType = WATER;
						}
						else
							break;
					}
				}
				else if (cell == 0)
				{
					ray.pos = p + result.normal * 0.001f;
					materialType = AIR;
				}

				continue;
			}
		}

		float3 colour = factor * (fogColour + SUN_COLOUR * 2.5f * pow(clamp(-dot(ray.dir, sunDirection), 0.0f, 1.0f) + 0.014f, 500.0f));

		if (materialType == WATER) 
			final += colour * (float3)(exp(-1.0f * result.t * (float3)(1.0f, 0.7f, 0.4f) * 0.8f));
		else		
			final += colour;

		break;
	}

	return final;
}

float3 traceOne(const Ray startRay, uint2* rngState, const float3 sunDirection, 
	global const uchar *voxels, const int3 size,
	global const Node *nodes, global const uint *ptrTable, 
	const float maxt, const int iterations, const int startMaterialType)
{
	const float3 ambientColour = (float3)(0.10f, 0.16f, 0.20f);

	Ray ray = startRay;

	const float3 skyColour = SKY_BASE_COLOUR * clamp(-ray.dir.y, 0.1f, 1.0f);
	const float3 sunFillColour = SUN_COLOUR * clamp(-dot(ray.dir, sunDirection) * 0.2f + 0.2f, 0.0f, 1.0f) * 2.5f;
	const float3 fogColour = skyColour + sunFillColour + ambientColour;

	float3 final = (float3)0.0f;
	float3 factor = (float3)1.0f; 
	int materialType = startMaterialType;

	Result result;

	if (materialType == AIR)
		result = intersectOctree(voxels, size, nodes, ptrTable, ray, maxt);
	else
		result = intersectVoxels(voxels, size, ray, maxt, materialType);

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

			float3 colour = (float3)(0.0f, 0.0f, 0.0f);

			if (cell != 4)
			{
				float diffuse = dot(result.normal, sunDirection);

				if (diffuse > 0.0f)
				{
					Ray shadowRay;
					shadowRay.pos = p - sunDirection * 0.001f;
					shadowRay.dir = -sunDirection;
					shadowRay.invdir = native_recip(shadowRay.dir);

					if (intersectVoxelsShadow(voxels, size, shadowRay, INITIAL_MAX))
						diffuse = 0.0f;
				}

				colour += SUN_COLOUR * 2.5f * clamp(diffuse, 0.0f, 1.0f);

				float u1 = nextRandomFloat(rngState);
				float u2 = nextRandomFloat(rngState);

				Ray aoRay;
				aoRay.pos = p - result.normal * 0.001f;
				aoRay.dir = -sampleHemisphere(u1, u2, result.normal);
				aoRay.invdir = native_recip(aoRay.dir);

				if (!intersectVoxelsShadow(voxels, size, aoRay, 8.0f)) 
					colour += ambientColour;
			}

			if (materialType == WATER) // Apply Beer's law
				factor *= exp(-1.0f * result.t * (float3)(1.0f, 0.7f, 0.4f) * 0.8f);

			final += factor * mix(fogColour, colour * blockColour, fog);
		}

		return final;
	}
	else
	{
		float3 colour = fogColour + SUN_COLOUR * 2.5f * pow(clamp(-dot(ray.dir, sunDirection), 0.0f, 1.0f) + 0.014f, 500.0f);
		
		if (materialType == WATER) 
			return colour * (float3)(exp(-1.0f * result.t * (float3)(1.0f, 0.7f, 0.4f) * 0.8f));
		else		
			return colour;
	}
}

float3 traceNoGI(const Ray ray, uint2* rngState, const float3 sunDirection, 
	global const uchar *voxels, const int3 size,
	global const Node *nodes, global const uint *ptrTable, 
	const float maxt)
{
	const float3 skyBaseColour = (float3)(0.05f, 0.1f, 0.3f);
	const float3 sunColour = (float3)(0.63f, 0.55f, 0.5f);
	const float3 ambientColour = (float3)(0.10f, 0.16f, 0.20f);

	const float3 skyColour = skyBaseColour * clamp(-ray.dir.y, 0.1f, 1.0f);
	const float3 sunFillColour = sunColour * clamp(-dot(ray.dir, sunDirection) * 0.2f + 0.2f, 0.0f, 1.0f) * 2.5f;
	const float3 fogColour = skyColour + sunFillColour + ambientColour;

	const Result result = intersectOctree(voxels, size, nodes, ptrTable, ray, maxt);
	//const Result result = intersectVoxels(voxels, size, ray, maxt, AIR);

    float3 colour = fogColour;

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

			/*float2 coord = (float2)(p.x, p.z);
			if (result.side == 0) coord = (float2)(p.y, p.z);
			if (result.side == 2) coord = (float2)(p.x, p.y);

			float3 blockColour = getCellColourSimple(cell);
			if (cell == 1) blockColour = pow(read_imagef(grassTexture, textureSampler, coord * 0.5f).xyz, 2.2f);
			*/

			colour = (float3)(0.0f, 0.0f, 0.0f);

			if (cell != 4)
			{
		    	float diffuse = dot(result.normal, sunDirection);
		    	if (diffuse > 0.0f)
				{
					Ray shadowRay;
					shadowRay.pos = p - sunDirection * 0.001f;
					shadowRay.dir = -sunDirection;
					shadowRay.invdir = native_recip(shadowRay.dir);
					if (intersectOctreeShadow(voxels, size, nodes, ptrTable, shadowRay, maxt)) diffuse = 0.0f;
				}

				colour += sunColour * 2.5f * clamp(diffuse, 0.0f, 1.0f);
			}

			if (cell == 4)
			{
				const float airDensity = 1.0f;
				const float waterDensity = 1.33f;

				const float k = 1.0f - dot(ray.dir, result.normal);

				const float a = waterDensity - airDensity;
				const float b = waterDensity + airDensity;
				const float refl0 = (a*a)/(b*b);

				const float reflectance = refl0 + (1.0f - refl0) * k * k * k * k;

				const int reflectionSamples = 2;
				const float invRFL = 1.0f / reflectionSamples;

				#pragma unroll
				for (int i = 0; i < reflectionSamples; i++)
				{
					float u1 = nextRandomFloat(rngState) * 0.001f;
					float u2 = (i + nextRandomFloat(rngState)) * invRFL;

					Ray giRay;
					giRay.pos = p - result.normal * 0.001f;
					giRay.dir = sampleHemisphere(u1, u2, reflect(result.normal, ray.dir));
					giRay.invdir = native_recip(giRay.dir);

					colour += traceOne(giRay, rngState, sunDirection, voxels, size, nodes, ptrTable, maxt - result.t, 2, AIR) * invRFL * reflectance;
				}

				const int refractionSamples = 1;
				const float invRFR = 1.0f / refractionSamples;

				#pragma unroll
				for (int i = 0; i < refractionSamples; i++)
				{
					float u1 = nextRandomFloat(rngState) * 0.001f;
					float u2 = (i + nextRandomFloat(rngState)) * invRFR;

					Ray tRay;
					tRay.pos = p + result.normal * 0.001f;
					tRay.dir = sampleHemisphere(u1, u2, ray.dir);
					tRay.invdir = native_recip(tRay.dir);

					colour += traceOne(tRay, rngState, sunDirection, voxels, size, nodes, ptrTable, 20.0f, 2, WATER) * invRFR * (1.0f - reflectance);
				}

				colour = mix(fogColour, colour, fog);
			}
	   		else
	   		{
				float ambientFactor = 1.0f;

				const int aoSamplesU = 3;
				const float invAOu = 1.0f / aoSamplesU;

				const int aoSamplesV = 3;
				const float invAOv = 1.0f / aoSamplesV;

				for (int i = 0; i < aoSamplesU; i++)
					for (int j = 0; j < aoSamplesV; j++)
					{
						float u1 = (i + nextRandomFloat(rngState)) * invAOu;
						float u2 = (j + nextRandomFloat(rngState)) * invAOv;

						Ray aoRay;
						aoRay.pos = p - result.normal * 0.01f;
						aoRay.dir = -sampleHemisphere(u1, u2, result.normal);
						aoRay.invdir = native_recip(aoRay.dir);

						//if (intersectOctreeShadow(voxels, size, nodes, ptrTable, aoRay, 3.0f)) ambientFactor -= invAOu * invAOv;
						if (intersectVoxelsShadow(voxels, size, aoRay, 12.0f)) ambientFactor -= invAOu * invAOv;
					}

	    		colour += ambientColour * ambientFactor;
	    		colour = mix(fogColour, colour * blockColour, fog);
			}
		}
	}
	else
		colour += sunColour * 2.5f * pow(clamp(-dot(ray.dir, sunDirection), 0.0f, 1.0f) + 0.014f, 500.0f);

	return colour;
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

    const int seed = initial_seed + coords.x * 6971 + coords.y * 31337;

	uint2 rngState = {noise1i(0, seed), initial_seed};

    const int3 size = *voxelsSize;
    const float3 sunDirection = *sunDir;
    const float fov = 75.0f;

    const float inv_hheight = 2.0f / height;
    const float x = coords.x * inv_hheight - (float)width / height;
    const float y = coords.y * inv_hheight - 1.0f;
    const float r = native_tan(radians(fov * 0.5f));

    float4 screen = (float4)(x*r, y*r, 1.0f, 1.0f);

    int p = 4;

	Ray ray;
    ray.pos = *cameraPos;
    ray.dir = fast_normalize((float3)(dot(cameraViewMatrix[0], screen), dot(cameraViewMatrix[1], screen), dot(cameraViewMatrix[2], screen)));
    ray.invdir = native_recip(ray.dir);

	float3 result = traceNoGI(ray, &rngState, sunDirection, voxels, size, nodes, ptrTable, INITIAL_MAX);

    write_imagef(output, coords, (float4)(result, 0.0f));
}