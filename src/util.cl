constant const float FOG_DENSITY = 0.003f;
constant const float INITIAL_MAX = 512.0f;
constant const float SUN_SIZE = 0.00925f;
constant const float FOV = 75.0f;

constant const float3 SKY_BASE_COLOUR = (float3)(0.05f, 0.1f, 0.3f);
constant const float3 SUN_COLOUR = (float3)(0.63f, 0.55f, 0.5f);

constant const sampler_t textureSampler = CLK_NORMALIZED_COORDS_TRUE | CLK_ADDRESS_REPEAT | CLK_FILTER_NEAREST;

float3 reflect(float3 normal, float3 in)
{
	return in - normal * dot(normal, in) * 2.0f;
}

float3 getCellColourSimple(int cell)
{
	float3 blockColour;

	if (cell == 1)
		blockColour = (float3)(0.65f, 0.9f, 0.2f);
	else if (cell == 2)
		blockColour = (float3)(0.4f, 0.2f, 0.05f);
	else if (cell == 3)
		blockColour = (float3)(0.4f, 0.4f, 0.4f);
	else if (cell == 4)
		blockColour = (float3)(0.2f, 0.35f, 0.6f);
	else if (cell == 5)
		blockColour = (float3)(0.25f, 0.13f, 0.03f);
	else if (cell == 6)
		blockColour = (float3)(0.4f, 0.8f, 0.05f);
	else
		blockColour = (float3)(0.0f, 0.0f, 0.0f);

	return blockColour * 0.4f;
}

float3 getCellColour(int cell, float3 p)
{
	const float tl = perlinNoise(0.7, 7, p * 0.33f, 1) * 0.7f + 0.15f;

	float3 blockColour;

	if (cell == 1)
		blockColour = mix((float3)(0.6f, 0.95f, 0.1f), (float3)(0.7f, 0.85f, 0.3f), perlinNoise(0.5, 5, p, 2)) * tl;
	else if (cell == 2)
		blockColour = (float3)(0.4f, 0.2f, 0.05f) * tl;
	else if (cell == 3)
		blockColour = (float3)(0.4f, 0.4f, 0.4f) * tl;
	else if (cell == 4)
		blockColour = (float3)(0.2f, 0.35f, 0.6f) * tl;
	else if (cell == 5)
		blockColour = (float3)(0.25f, 0.13f, 0.03f) * tl;
	else if (cell == 6)
		blockColour = (float3)(0.4f, 0.8f, 0.05f) * tl;
	else
		blockColour = (float3)(0.0f, 0.0f, 0.0f);

	return blockColour;
}

float3 sampleHemisphere(float u1, float u2, float3 normal)
{
	const float3 u = normal;
	float3 v = (fabs(u.x) > 0.1f) ? (float3)(0.0f, 1.0f, 0.0f) : (float3)(1.0f, 0.0f, 0.0f);
	const float3 w = fast_normalize(cross(u, v));
	v = cross(w, u);

	const float r = native_sqrt(u1);
	const float theta = 2.0f * M_PI * u2;
	const float x = r * native_cos(theta);
	const float y = r * native_sin(theta);
	return fast_normalize(u * native_sqrt(1.0f - u1) + v * x + w * y);
}

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

float3 tonemap(float3 input)
{
	//return input * 0.2f;

	/*float3 x = max(0.0f, input - 0.004f) * 0.05f;
	float3 output = (x*(6.2f*x + 0.5f)) / (x*(6.2f*x + 1.7f) + 0.06f);
	return output;*/

	float3 Yxy = convertXYZtoYxy(convertRGBtoXYZ(input));

	float x = max(0.0f, Yxy.s0 - 0.004f) * 0.1f;
	Yxy.s0 = native_divide(x*(6.2f*x + 0.5f), x*(6.2f*x + 1.7f) + 0.06f);

	return convertXYZtoRGB(convertYxytoXYZ(Yxy));
}



