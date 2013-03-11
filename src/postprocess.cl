#define OPENCL
#include "src/random.cl"

__kernel void main(const int initial_seed, const float mixFactor, read_only image2d_t input1, read_only image2d_t input2, write_only image2d_t output)
{
	const int width = get_image_width(output);
	const int height = get_image_height(output);

	const int2 coords = (int2)(get_global_id(0), get_global_id(1));

	if (coords.x > width) return;
	if (coords.y > height) return;

	const uint seed = initial_seed + coords.x * 6971 + coords.y * 31337;
	uint2 rngState = {noise1i(0, seed), initial_seed};

	const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;

	float ditherNoise = noise2(coords.x, coords.y, seed) * 0.004f - 0.002f;

	float3 colour = read_imagef(input1, sampler, (float2)(coords.x+0.5f, coords.y+0.5f)).xyz + (float3)(ditherNoise);

	/*
	float off = noise2(coords.x, coords.y, seed) * M_PI * 2.0f;
	float gauss[11] = {0.0085, 0.0278, 0.0667, 0.1222, 0.1746, 0.1964, 0.1746, 0.1222, 0.0667, 0.0278, 0.0085};

	for (int x = -5; x <= 5; x++)
		for (int y = -5; y <= 5; y++)
			colour += max(read_imagef(input1, sampler, coords + (int2)(x * sin(off) * 4.0f, y * cos(off) * 4.0f)).xyz, (float3)(0.0f)) * gauss[x+5] * gauss[y+5] * 0.2f;
	*/

	if (mixFactor > 0.0f)
		colour = colour * (1.0f - mixFactor) + read_imagef(input2, sampler, (float2)(coords.x+0.5f, coords.y+0.5f)).xyz * mixFactor;

	write_imagef(output, coords, (float4)(colour, 0.0f));
}