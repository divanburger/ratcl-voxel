#define OPENCL
#include "src/random.cl"

// Tonemapping
float3 convertRGBtoXYZ(float3 RGB)
{
	return (float3)(
		dot(RGB, (float3)(0.5141364f, 0.3238786f,  0.16036376f)),
		dot(RGB, (float3)(0.265068f,  0.67023428f, 0.06409157f)),
		dot(RGB, (float3)(0.0241188f, 0.1228178f,  0.8444266f))
	);
}

float3 convertXYZtoRGB(float3 XYZ)
{
	return (float3)(
		dot(XYZ, (float3)(2.5651f,-1.1665f,-0.3986f)),
		dot(XYZ, (float3)(-1.0217f, 1.9777f, 0.0439f)),
		dot(XYZ, (float3)(0.0753f, -0.2543f, 1.1892f))
	);
}

float3 convertXYZtoYxy(float3 XYZ)
{
	float f = native_recip(XYZ.x + XYZ.y + XYZ.z);
	return (float3)(XYZ.y, XYZ.xy * f);
}

float3 convertYxytoXYZ(float3 Yxy)
{
	float f = native_recip(Yxy.z);
	return (float3)(Yxy.y * f, 1.0f, (1.0f - Yxy.y - Yxy.z) * f) * Yxy.x;
}

float3 tonemap(const float3 input, const float exposure)
{
	float3 Yxy = convertXYZtoYxy(convertRGBtoXYZ(input));

	float x = max(0.0f, Yxy.s0 - 0.004f) * exposure;
	Yxy.s0 = native_divide(x*(6.2f*x + 0.5f), x*(6.2f*x + 1.7f) + 0.06f);

	return convertXYZtoRGB(convertYxytoXYZ(Yxy));
}

__kernel void main(read_only image2d_t input, write_only image2d_t output, const float exposure)
{
	const int width = get_image_width(output);
	const int height = get_image_height(output);

	const int2 coords = (int2)(get_global_id(0), get_global_id(1));

	if (coords.x > width) return;
	if (coords.y > height) return;

	const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

	float3 colour = read_imagef(input, sampler, coords).xyz;

	write_imagef(output, coords, (float4)(tonemap(colour, exposure), 0.0f));
}