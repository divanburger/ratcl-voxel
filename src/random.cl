uint MWC64X(uint2* state)
{
	const uint A = 4294883355U;
	uint x = (*state).x, c = (*state).y;
	uint res = x^c;
	uint hi = mul_hi(x,A);
	x = x*A + c;
	c = hi + (x<c);
	*state = (uint2)(x, c);
	return res;
}

float nextRandomFloat(uint2* state)
{
	return (float)(MWC64X(state)) / (float)(0xFFFFFFFF);
}

int noise1i(const int x, const uint seed)
{
	int n = x * 1619 + seed * 6971;
	n = (n>>8)^n;
	return n * (n * n * 60493 + 19990303) + 1376312589;
}

float noise1(const int x, const uint seed)
{
	int n = x * 1619 + seed * 6971;
	n = (n>>8)^n;
	return convert_float((n * (n * n * 60493 + 19990303) + 1376312589) & 0x7fffffff) / 2147483647.0f;
}

float noise2(const int x, const int y, const uint seed)
{
	int n = x * 1619 + y * 31337 + seed * 6971;
	n = (n>>8)^n;
	return convert_float((n * (n * n * 60493 + 19990303) + 1376312589) & 0x7fffffff) / 2147483647.0f;
}

float noise3(const int x, const int y, const int z, const uint seed)
{
	int n = x * 1619 + y * 31337 + z + seed * 6971;
	n = (n>>8)^n;
	return convert_float((n * (n * n * 60493 + 19990303) + 1376312589) & 0x7fffffff) / 2147483647.0f;
}

float interpolatedNoise(const float3 pos, const uint seed)
{
	const int3 ipos = convert_int3(pos);
	const float3 fpos = pos - convert_float3(ipos);

	float NY1 = mix(noise3(ipos.x,ipos.y  ,ipos.z,seed), noise3(ipos.x+1,ipos.y  ,ipos.z,seed), fpos.x);
	float NY2 = mix(noise3(ipos.x,ipos.y+1,ipos.z,seed), noise3(ipos.x+1,ipos.y+1,ipos.z,seed), fpos.x);

	const float NZ1 = mix(NY1, NY2, fpos.y);

	NY1 = mix(noise3(ipos.x,ipos.y  ,ipos.z+1,seed), noise3(ipos.x+1,ipos.y  ,ipos.z+1,seed), fpos.x);
	NY2 = mix(noise3(ipos.x,ipos.y+1,ipos.z+1,seed), noise3(ipos.x+1,ipos.y+1,ipos.z+1,seed), fpos.x);

	const float NZ2 = mix(NY1, NY2, fpos.y);

	return mix(NZ1, NZ2, fpos.z);
}

float perlinNoise(const float persistence, const int octaves, const float3 pos, const uint seed)
{
	float frequency = 1;
	float amplitude = 1;
	float total = 0;
	float max = 0;

	for (int i = 0; i < (octaves - 1); i++)
	{
		if (i > 0)
		{
			frequency *= 2;
			amplitude *= persistence;
		}

		max += amplitude;
		total += interpolatedNoise(pos * frequency, seed + i) * amplitude;
	}

	return (total / max);
}