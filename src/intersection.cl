constant const int AIR = 0;
constant const int WATER = 4;
constant const int WOOD = 5;
constant const int LEAVES = 6;

typedef struct {
	bool hit;
	int side;
	float t;
	float3 normal;
	float3 pos;
} Result;

typedef struct {
	ushort	children;
	uchar	child;
	uchar	leaf;
} Node;

Result intersectOctree(global const uchar *voxels, const int3 size, global const Node *nodes, global const uint *ptrTable, Ray ray, float maxt)
{
	Result result;
	result.normal = (float3)(1.0f, 1.0f, 1.0f);
	result.hit = false;

	float tmin = 0.0f;
	float tmax = maxt;

	int temp = size.x;
	int scale = 0;

	while (temp >>= 1) ++scale; 

	int hscale = scale - 1;

	float3 p = ray.pos;

	int x = p.x, y = p.y, z = p.z;
	x = clamp(x, 0, size.x);
	y = clamp(y, 0, size.y);
	z = clamp(z, 0, size.z);

	uchar	rayMask = 0;
	if (ray.dir.x > 0) rayMask |= 1;
	if (ray.dir.y > 0) rayMask |= 2;
	if (ray.dir.z > 0) rayMask |= 4;

	typedef struct Entry
	{
		uint	id;
		uchar	octantMask;
	} Entry;

	int sp = 0;
	Entry stck[12];
	{
		uchar	octantMask = 0;
		if (x & (1<<hscale)) octantMask |= 1;
		if (y & (1<<hscale)) octantMask |= 2;
		if (z & (1<<hscale)) octantMask |= 4;

		stck[0].id = 0U;
		stck[0].octantMask = octantMask;
	}

	int side = 0;

	Entry* entry = &stck[sp];

	for (int iter = 0; iter < 200 && sp < 12; iter++)
	{
		const Node node = nodes[entry->id];

		// Check if the current octant has child
		uchar childID = 1U << entry->octantMask;

		if (node.child & childID)
		{
			// Check if it is a leaf node
			if (node.leaf & childID) 
			{
				float sgn = ((rayMask >> side) & 1) ? 1.0f : -1.0f;

				result.hit = true;
				result.t = tmin;
				result.side = side;
				result.normal = (float3)(side == 0 ? 1.0f : 0.0f, side == 1 ? 1.0f : 0.0f, side == 2 ? 1.0f : 0.0f) * sgn;
				result.pos = ray.pos + ray.dir * tmin + result.normal * 0.0001f;
				result.pos.x = min(max((int)result.pos.x, 0), size.x);
				result.pos.y = min(max((int)result.pos.y, 0), size.y);
				result.pos.z = min(max((int)result.pos.z, 0), size.z);
				return result;
			}
			
			//// Push
			int hhscale = hscale - 1;

			uchar	octantMask = 0;
			if (x & (1<<hhscale)) octantMask |= 1;
			if (y & (1<<hhscale)) octantMask |= 2;
			if (z & (1<<hhscale)) octantMask |= 4;

			// Decode pointer
			uint children = node.children;
			if (node.children & 0x8000)
				children = ptrTable[node.children & 0x7FFF];

			// Add to stack
			++sp;
			stck[sp].id = entry->id + children + entry->octantMask;
			stck[sp].octantMask = octantMask;
			entry = &stck[sp];

			scale = hscale;
			hscale = scale - 1;

			continue;
		}

		//// Advance
		const int nx = ((x >> hscale) << hscale) + (rayMask & 1 ? (1<<hscale) : 0);
		const int ny = ((y >> hscale) << hscale) + (rayMask & 2 ? (1<<hscale) : 0);
		const int nz = ((z >> hscale) << hscale) + (rayMask & 4 ? (1<<hscale) : 0);

		const float3 nextt = ((float3)(nx, ny, nz) - ray.pos) * ray.invdir;

		tmin = max(min(min(nextt.x, nextt.y), nextt.z), tmin);
		if (tmin >= tmax) break;

		const float3 cp = ray.pos + ray.dir * tmin;
		x = clamp((int)(cp.x), 0, size.x);
		y = clamp((int)(cp.y), 0, size.y);
		z = clamp((int)(cp.z), 0, size.z);

		uchar stepMask = 0;
		if (nextt.x <= tmin) {stepMask |= 1; side = 0; x = nx - (rayMask & 1 ? 0 : 1);}
		if (nextt.y <= tmin) {stepMask |= 2; side = 1; y = ny - (rayMask & 2 ? 0 : 1);}
		if (nextt.z <= tmin) {stepMask |= 4; side = 2; z = nz - (rayMask & 4 ? 0 : 1);}
		
		//// Pop
		while (true)
		{
			const uchar endMask = ~(entry->octantMask ^ rayMask);
			const uchar popMask = endMask & stepMask;
			entry->octantMask ^= stepMask;

			if (!popMask) break;

			// If we exit the top level we done
			if (sp == 0) return result;

			// Pop of the stack
			entry = &stck[--sp];

			hscale = scale;
			scale = hscale + 1;
		}
	}

	return result;
}

bool intersectOctreeShadow(global const uchar *voxels, const int3 size, global const Node *nodes, global const uint *ptrTable, const Ray ray, const float maxt)
{
	float tmin = 0.0f;
	float tmax = maxt;

	int temp = size.x;
	int scale = 0;

	while (temp >>= 1) ++scale; 

	int hscale = scale - 1;

	float3 p = ray.pos;
	int x = p.x, y = p.y, z = p.z;
	x = clamp(x, 0, size.x);
	y = clamp(y, 0, size.y);
	z = clamp(z, 0, size.z);

	uchar	rayMask = 0;
	if (ray.dir.x > 0) rayMask |= 1;
	if (ray.dir.y > 0) rayMask |= 2;
	if (ray.dir.z > 0) rayMask |= 4;

	typedef struct Entry
	{
		uint	id;
		uchar	octantMask;
	} Entry;

	int sp = 0;
	Entry stck[12];
	{
		uchar	octantMask = 0;
		if (x & (1<<hscale)) octantMask |= 1;
		if (y & (1<<hscale)) octantMask |= 2;
		if (z & (1<<hscale)) octantMask |= 4;

		stck[0].id = 0U;
		stck[0].octantMask = octantMask;
	}

	Entry* entry = &stck[sp];

	for (int iter = 0; iter < 200 && sp < 12; iter++)
	{
		const Node node = nodes[entry->id];

		// Check if the current octant has child
		uchar childID = 1U << entry->octantMask;

		if (node.child & childID)
		{
			// Check if it is a leaf node
			if (node.leaf & childID) return true;
			
			int hhscale = hscale - 1;

			//// Push
			uchar	octantMask = 0;
			if (x & (1<<hhscale)) octantMask |= 1;
			if (y & (1<<hhscale)) octantMask |= 2;
			if (z & (1<<hhscale)) octantMask |= 4;

			// Decode pointer
			uint children = node.children;
			if (node.children & 0x8000)
				children = ptrTable[node.children & 0x7FFF];

			// Add to stack
			++sp;
			stck[sp].id = entry->id + children + entry->octantMask;
			stck[sp].octantMask = octantMask;
			entry = &stck[sp];

			scale = hscale;
			hscale = scale - 1;

			continue;
		}

		//// Advance
		int nx = ((x >> hscale) << hscale) + (rayMask & 1 ? (1<<hscale) : 0);
		int ny = ((y >> hscale) << hscale) + (rayMask & 2 ? (1<<hscale) : 0);
		int nz = ((z >> hscale) << hscale) + (rayMask & 4 ? (1<<hscale) : 0);

		float3 nextt = ((float3)(nx, ny, nz) - ray.pos) * ray.invdir;

		tmin = max(min(min(nextt.x, nextt.y), nextt.z), tmin);
		if (tmin >= tmax) break;

		float3 cp = ray.pos + ray.dir * tmin;
		x = clamp((int)(cp.x), 0, size.x);
		y = clamp((int)(cp.y), 0, size.y);
		z = clamp((int)(cp.z), 0, size.z);

		uchar stepMask = 0;
		if (nextt.x <= tmin) {stepMask |= 1; x = nx - (rayMask & 1 ? 0 : 1);}
		if (nextt.y <= tmin) {stepMask |= 2; y = ny - (rayMask & 2 ? 0 : 1);}
		if (nextt.z <= tmin) {stepMask |= 4; z = nz - (rayMask & 4 ? 0 : 1);}

		//// Pop
		while (true)
		{
			const uchar endMask = ~(entry->octantMask ^ rayMask);
			const uchar popMask = endMask & stepMask;
			entry->octantMask ^= stepMask;

			if (!popMask) break;

			// If we exit the top level we done
			if (sp == 0 && popMask) return false;

			// Pop of the stack
			entry = &stck[--sp];

			hscale = scale;
			scale = hscale + 1;
		}
	}

	return false;
}

bool intersectVoxelsShadow(global const uchar *voxels, const int3 size, const Ray ray, const float maxt)
{
	int3 p = clamp(convert_int(ray.pos), 0, size);
	float t = 0.0f;

	float3 dsign = (float3)(sign(ray.dir.x), sign(ray.dir.y), sign(ray.dir.z));
	const float3 next = (float3)(ray.dir.x >= 0 ? p.x + 1 : p.x, ray.dir.y >= 0 ? p.y + 1 : p.y, ray.dir.z >= 0 ? p.z + 1 : p.z);

	float3 tmax = (next - ray.pos) * ray.invdir, tdelta = dsign * ray.invdir;

	for (int l = 0; l < maxt * 2; l++)
	{
		int cell = voxels[p.x + (p.y + p.z * size.y) * size.x];

		if (cell != AIR && cell != WATER) return true;

		float ttmax = min(min(tmax.x, tmax.y), tmax.z);
		if (ttmax > maxt) break;

		if (tmax.x == ttmax)
		{
			p.x += dsign.x;
			tmax.x += tdelta.x;
		}
		else if (tmax.y == ttmax)
		{
			p.y += dsign.y;
			tmax.y += tdelta.y;
		}
		else if (tmax.z == ttmax)
		{
			p.z += dsign.z;
			tmax.z += tdelta.z;
		}

		if (any(p < 0) || any(p >= size)) break;
	}

	return false;
}

Result intersectVoxels(global const uchar *voxels, const int3 size, const Ray ray, const float maxt, const int skipType)
{
	Result result;
	result.t = maxt;
	result.normal = (float3)(1.0f, 1.0f, 1.0f);
	result.hit = false;

	int3 p = clamp(convert_int(ray.pos), 0, size);
	float t = 0.0f;

	const float3 dsign = (float3)(sign(ray.dir.x), sign(ray.dir.y), sign(ray.dir.z));
	const float3 next = (float3)(ray.dir.x >= 0 ? p.x + 1 : p.x, ray.dir.y >= 0 ? p.y + 1 : p.y, ray.dir.z >= 0 ? p.z + 1 : p.z);

	float3 tmax = (next - ray.pos) * ray.invdir, tdelta = dsign * ray.invdir;

	int side = 0;

	for (int l = 0; l < convert_int(maxt * 1.5f); l++)
	{
		const int cell = voxels[p.x + (p.y + p.z * size.y) * size.x];

		if (cell != skipType)
		{
			float sgn = (side == 0) ? dsign.s0 : (side == 1 ? dsign.s1 : dsign.s2);

			result.pos = convert_float3(p);
			result.hit = true;
			result.side = side;
			result.normal = (float3)(side == 0 ? 1.0f : 0.0f, side == 1 ? 1.0f : 0.0f, side == 2 ? 1.0f : 0.0f) * dsign;
			return result;
		}

		result.t = min(min(tmax.x, tmax.y), tmax.z);
		if (result.t > maxt) break;

		if (tmax.x == result.t)
		{
			p.x += dsign.x;
			tmax.x += tdelta.x;
			side = 0;
		}
		else if (tmax.y == result.t)
		{
			p.y += dsign.y;
			tmax.y += tdelta.y;
			side = 1;
		}
		else if (tmax.z == result.t)
		{
			p.z += dsign.z;
			tmax.z += tdelta.z;
			side = 2;
		}

		if (any(p < 0) || any(p >= size)) break;
	}

	return result;
}