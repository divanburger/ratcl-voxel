#include "Noise.h"

float noise(int x, int seed)
{
	int n = x * 1619 + seed * 6971;
	n = (n>>8)^n;
	return ((n * (n * n * 60493 + 19990303) + 1376312589) & 0x7fffffff) / 2147483647.0f; 
}

float noise(int x, int y, int seed)
{
	int n = x * 1619 + y * 31337 + seed * 6971;
	n = (n>>8)^n;
	return ((n * (n * n * 60493 + 19990303) + 1376312589) & 0x7fffffff) / 2147483647.0f;
}

float noise(int x, int y, int z, int seed)
{
	int n = x * 1619 + y * 31337 + z + seed * 6971;
	n = (n>>8)^n;
	return ((n * (n * n * 60493 + 19990303) + 1376312589) & 0x7fffffff) / 2147483647.0f;
}

inline float smoothNoise(int x, int seed)
{
	return (noise(x-1,seed) + noise(x,seed) * 2 + noise(x+1,seed)) / 4;
}

inline float smoothNoise(int x, int y, int seed)
{
	float Total = 0;
	
	for (int ix = -1; ix < 1; ix++)
		for (int iy = -1; iy < 1; iy++)
		{
			int zerocount = 0;
			if (ix == 0) zerocount++;
			if (iy == 0) zerocount++;
			
			switch (zerocount)
			{
				case 0: Total += noise(x+ix,y+iy,seed); break;      
				case 1: Total += noise(x+ix,y+iy,seed) * 2; break;
				case 2: Total += noise(x+ix,y+iy,seed) * 4; break;   
			}
		}
	
	return Total / 16;
}

inline float smoothNoise(int x, int y, int z, int seed)
{
	float Total = 0;
	
	for (int ix = -1; ix < 1; ix++)
		for (int iy = -1; iy < 1; iy++)
			for (int iz = -1; iz < 1; iz++)
			{
				int zerocount = 0;
				if (ix == 0) zerocount++;
				if (iy == 0) zerocount++;
				if (iz == 0) zerocount++;
				
				switch (zerocount)
				{
					case 0: Total += noise(x+ix,y+iy,z+iz,seed); break;      
					case 1: Total += noise(x+ix,y+iy,z+iz,seed) * 2; break;
					case 2: Total += noise(x+ix,y+iy,z+iz,seed) * 4; break;   
					case 3: Total += noise(x+ix,y+iy,z+iz,seed) * 8; break;
				}
			}
	
	return Total / 48;
}

float interpolatedNoise(float x, int seed)
{
	int ix = int(x);

	float gx = x - ix;
	
	gx = 3*gx*gx - 2*gx*gx*gx;
	
	float N1 = noise(ix,seed);
	float N2 = noise(ix+1,seed);
	
	return smoothmix(N1, N2, gx);
}

float interpolatedNoise(float x, float y, int seed)
{
	int ix = int(x);
	int iy = int(y);

	float gx = x - ix;
	float gy = y - iy;
	
	gx = 3*gx*gx - 2*gx*gx*gx;
	gy = 3*gy*gy - 2*gy*gy*gy;
	
	float NY1 = smoothmix(noise(ix,iy,seed), noise(ix+1,iy,seed), gx);
	float NY2 = smoothmix(noise(ix,iy+1,seed), noise(ix+1,iy+1,seed), gx);
	
	return smoothmix(NY1, NY2, gy);
}

float interpolatedNoise(float x, float y, float z, int seed)
{
	int ix = int(x);
	int iy = int(y);
	int iz = int(z);
	
	float NY1 = smoothmix(noise(ix,iy,iz,seed), noise(ix+1,iy,iz,seed), x - ix);
	float NY2 = smoothmix(noise(ix,iy+1,iz,seed), noise(ix+1,iy+1,iz,seed), x - ix);
	
	float NZ1 = smoothmix(NY1, NY2, y - iy);

	NY1 = smoothmix(noise(ix,iy,iz+1,seed), noise(ix+1,iy,iz+1,seed), x - ix);
	NY2 = smoothmix(noise(ix,iy+1,iz+1,seed), noise(ix+1,iy+1,iz+1,seed), x - ix);
	
	float NZ2 = smoothmix(NY1, NY2, y - iy);
	
	return smoothmix(NZ1, NZ2, z - iz);
}

float perlinNoise(float persistence, int octaves, float x, int seed)
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
		total += interpolatedNoise(x * frequency, seed + i) * amplitude;
	}
	
	return (total / max);
}

float perlinNoise(float persistence, int octaves, float x, float y, int seed)
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
		total += interpolatedNoise(x * frequency, y * frequency, seed + i) * amplitude;
	}
	
	return (total / max);
}

float perlinNoise(float persistence, int octaves, float x, float y, float z, int seed)
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
		total += interpolatedNoise(x * frequency, y * frequency, z * frequency, seed + i) * amplitude;
	}
	
	return (total / max);
}
