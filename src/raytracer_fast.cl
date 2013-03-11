#define OPENCL
#include "src/Entities.h"
#include "src/random.cl"
#include "src/intersection.cl"
#include "src/util.cl"

float3 trace(const Ray ray, uint2* rngState, const float3 sunDirection, 
	global const uchar *voxels, const int3 size,
	global const Node *nodes, global const uint *ptrTable,
	const float maxt)
{
	const float3 sunColour = SUN_COLOUR * 2.5f;
	const float3 ambientColour = (float3)(0.10f, 0.16f, 0.20f);

	const float3 skyColour = SKY_BASE_COLOUR * clamp(-ray.dir.y, 0.1f, 1.0f);
	const float3 sunFillColour = sunColour * clamp(-dot(ray.dir, sunDirection) * 0.2f + 0.2f, 0.0f, 1.0f);
	const float3 fogColour = skyColour + sunFillColour + ambientColour;

	//const Result result = intersectVoxels(voxels, size, ray, maxt, AIR);
	const Result result = intersectOctree(voxels, size, nodes, ptrTable, ray, maxt);

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

			float diffuse = dot(result.normal, sunDirection);
			if (diffuse > 0.0f)
			{
				Ray shadowRay;
				shadowRay.pos = p - result.normal * 0.001f;
				shadowRay.dir = -sunDirection;
				shadowRay.invdir = native_recip(shadowRay.dir);
				if (intersectOctreeShadow(voxels, size, nodes, ptrTable, shadowRay, maxt)) diffuse = 0.0f;
				//if (intersectVoxelsShadow(voxels, size, shadowRay, maxt)) diffuse = 0.0f;
			}

			colour = sunColour * clamp(diffuse, 0.0f, 1.0f) * 0.6f + ambientColour;
			colour = mix(fogColour, colour * blockColour, fog);
		}
	}
	else
		colour += sunColour * pow(clamp(-dot(ray.dir, sunDirection), 0.0f, 1.0f) + 0.014f, 500.0f);

	return colour;
}

kernel void raytrace(const uint initial_seed, write_only image2d_t output,
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

	const float inv_hheight = 2.0f / height;
	const float x = coords.x * inv_hheight - (float)width / height;
	const float y = coords.y * inv_hheight - 1.0f;
	const float r = native_tan(radians(FOV * 0.5f));

	const float4 screen = (float4)(x*r, y*r, 1.0f, 1.0f);

	int p = 1;

	Ray ray;
	ray.pos = *cameraPos;
	ray.dir = fast_normalize((float3)(dot(cameraViewMatrix[0], screen), dot(cameraViewMatrix[1], screen), dot(cameraViewMatrix[2], screen)));
	ray.invdir = native_recip(ray.dir);

	float3 result = trace(ray, &rngState, sunDirection, voxels, size, nodes, ptrTable, INITIAL_MAX);

	write_imagef(output, coords, (float4)(result, 0.0f));
}