typedef struct
{
	float3 pos;
	float d1;
	float3 dir;
	float d2;
	float3 invdir;
} Ray;

typedef struct
{
	float3 pos;
	float d;
	float yaw;
	float pitch;
} Camera;

typedef struct Scene
{
	float3	ambientColour;
	float3	lightColour;
	float3	lightDirection;
	float	fogDensity;
} Scene;