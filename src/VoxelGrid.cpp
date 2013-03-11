/*
 * VoxelGrid.cpp
 *
 *  Created on: 6 Nov 2011
 *      Author: dburger
 */

#include "VoxelGrid.h"

const vec3 VoxelGrid::normals[3] = {vec3(1.0f, 0.0f, 0.0f), vec3(0.0f, 1.0f, 0.0f), vec3(0.0f, 0.0f, 1.0f)};

VoxelGrid::VoxelGrid(ivec3 size) : size(size), aabb(vec3(), vec3(size))
{
	voxels = new uint8_t[size.x * size.y * size.z];

	nodes = new Node[1024*1024*1024];
	ptrTable = new uint32_t[0x7FFF];
}

VoxelGrid::~VoxelGrid()
{
	delete [] voxels;
}

bool VoxelGrid::intersectShadowOctree(Ray ray, float maxt) const
{
	return intersectOctree(ray, maxt).hit;
}

VoxelGrid::Result VoxelGrid::intersectOctree(Ray ray, float maxt, uint8_t materialType) const
{
	Result result;
	result.normal = vec3(1.0f, 1.0f, 1.0f);
	result.hit = false;

	//cout << endl << "========================" << endl << endl;

	float tmin = 0.0f;
	float tmax = numeric_limits<float>::max();//maxt;

	vec3 p = ray.pos;
	int x = p.x, y = p.y, z = p.z;
	x = clamp(x, 0, size.x);
	y = clamp(y, 0, size.y);
	z = clamp(z, 0, size.z);

	uint8_t	rayMask = 0;
	if (ray.dir.x > 0) rayMask |= 1;
	if (ray.dir.y > 0) rayMask |= 2;
	if (ray.dir.z > 0) rayMask |= 4;

	struct Entry
	{
		uint32_t	id;
		uint8_t		octantMask;
		uint8_t		scale;
	};

	int sp = 0;
	Entry stck[32];
	{
		int temp = size.x;
		int scale = 0;

		while (temp >>= 1) ++scale; 

		int hscale = scale - 1;

		uint8_t	octantMask = 0;
		if (x & (1<<hscale)) octantMask |= 1;
		if (y & (1<<hscale)) octantMask |= 2;
		if (z & (1<<hscale)) octantMask |= 4;

		//cout << "At: <" << x << ", " << y << ", " << z<< ">" << endl;
		//cout << "Octant Mask: " << (octantMask & 1 ? "1 " : "0 ") << (octantMask & 2 ? "2 " : "0 ") << (octantMask & 4 ? "4 " : "0 ") << endl;

		stck[0].id = materialType == 0 ? 0U : waterNodeIndex;
		stck[0].octantMask = octantMask;
		stck[0].scale = scale;
	}

	//cout << "Ray Mask: " << (rayMask & 1 ? "1 " : "0 ") << (rayMask & 2 ? "2 " : "0 ") << (rayMask & 4 ? "4 " : "0 ") << endl;

	int side = 0;
	iters = 0;

	Entry* entry = &stck[sp];
	while (sp < 32)
	{
		iters++;

		Node& node = nodes[entry->id];
		int hscale = entry->scale - 1;

		/*cout << endl;
		for (int k = 0; k < sp; k++)
		{
			Entry* e = &stck[k];
			cout << k << ": " << e->id << " (" << (int)e->scale << ")" << endl;
			cout << " - Octant Mask: " << (e->octantMask & 1 ? "1 " : "0 ") << (e->octantMask & 2 ? "2 " : "0 ") << (e->octantMask & 4 ? "4 " : "0 ") << endl;
		}*/

		/*cout << sp << ": " << entry->id << " (" << (int)entry->scale << ", <" << x << ", " << y << ", " << z<< ">)" << endl;
		cout << " - Octant Mask: " << (entry->octantMask & 1 ? "1 " : "0 ") << (entry->octantMask & 2 ? "2 " : "0 ") << (entry->octantMask & 4 ? "4 " : "0 ") << endl;
		cout << " - Value: " << get(
			x + (entry->octantMask & 1 ? 1: 0), 
			y + (entry->octantMask & 2 ? 1: 0), 
			z + (entry->octantMask & 4 ? 1: 0)
		) << endl;*/

		// Check if the current octant has child
		int8_t childID = 1U << entry->octantMask;
		if (node.child & childID)
		{
			// Check if it is a leaf node
			if (node.leaf & childID) 
			{
				result.hit = true;
				result.t = tmin;
				result.pos = ray.pos + ray.dir * (tmin + 1e-3f);
				result.pos.x = std::min(std::max((int)result.pos.x, 0), size.x);
				result.pos.y = std::min(std::max((int)result.pos.y, 0), size.y);
				result.pos.z = std::min(std::max((int)result.pos.z, 0), size.z);

				result.side = side;
				result.sign = ((rayMask >> side) & 1) ? 1.0f : -1.0f;
				result.normal = normals[side] * result.sign;

				return result;
			}
			else
			{
				int hhscale = hscale - 1;

				//cout << " - Push - " << (1<<hhscale) << endl;
				//cout << x << ", " << y << ", " << z << endl;

				//// Push
				uint8_t	octantMask = 0;
				if (x & (1<<hhscale)) octantMask |= 1;
				if (y & (1<<hhscale)) octantMask |= 2;
				if (z & (1<<hhscale)) octantMask |= 4;
				//cout << " - Octant Mask: " << (octantMask & 1 ? "1 " : "0 ") << (octantMask & 2 ? "2 " : "0 ") << (octantMask & 4 ? "4 " : "0 ") << endl;

				// Decode pointer
				uint32_t children = node.children;
				if (node.children & 0x8000)
					children = ptrTable[node.children & 0x7FFF];

				// Add to stack
				++sp;
				stck[sp].id = entry->id + children + entry->octantMask;
				stck[sp].octantMask = octantMask;
				stck[sp].scale = hscale;
				entry = &stck[sp];

				//cout << " -- At " << sp << endl;
				continue;
			}
		}

		//cout << " - Advance" << endl;

		//// Advance
		int nx = ((x >> hscale) << hscale) + (rayMask & 1 ? (1<<hscale) : 0);
		int ny = ((y >> hscale) << hscale) + (rayMask & 2 ? (1<<hscale) : 0);
		int nz = ((z >> hscale) << hscale) + (rayMask & 4 ? (1<<hscale) : 0);

		//cout << " -- Scale: " << hscale << endl;
		//cout << " -- Pos: <" << x << ", " << y << ", " << z << ">" << endl;
		//cout << " -- Next: <" << nx << ", " << ny << ", " << nz << ">" << endl;

		vec3 nextt = (vec3(nx, ny, nz) - ray.pos) * ray.invDir;
		//cout << " -- NextT: <" << nextt.x << ", " << nextt.y << ", " << nextt.z << ">" << endl;

		tmin = std::max(std::min(std::min(nextt.x, nextt.y), nextt.z), tmin);

		uint8_t stepMask = 0;
		if (nextt.x <= tmin) {stepMask |= 1; side = 0;}
		if (nextt.y <= tmin) {stepMask |= 2; side = 1;}
		if (nextt.z <= tmin) {stepMask |= 4; side = 2;}

		if (tmin >= tmax) break;

		//tmin += 1e-4f;

		//cout << " -- Step Mask: " << (stepMask & 1 ? "1 " : "0 ") << (stepMask & 2 ? "2 " : "0 ") << (stepMask & 4 ? "4 " : "0 ") << endl;

		vec3 cp = ray.pos + ray.dir * tmin;
		x = clamp((int)(cp.x), 0, size.x);
		y = clamp((int)(cp.y), 0, size.y);
		z = clamp((int)(cp.z), 0, size.z);

		if (nextt.x == tmin) x = nx - (rayMask & 1 ? 0 : 1);
		if (nextt.y == tmin) y = ny - (rayMask & 2 ? 0 : 1);
		if (nextt.z == tmin) z = nz - (rayMask & 4 ? 0 : 1);

		//// Pop
		while (sp >= 0)
		{
			uint8_t endMask = ~(entry->octantMask ^ rayMask);
			uint8_t popMask = endMask & stepMask;
			//cout << " -- Octant Mask: " << (entry->octantMask & 1 ? "1 " : "0 ") << (entry->octantMask & 2 ? "2 " : "0 ") << (entry->octantMask & 4 ? "4 " : "0 ") << endl;
			//cout << " -- End Mask: " << (endMask & 1 ? "1 " : "0 ") << (endMask & 2 ? "2 " : "0 ") << (endMask & 4 ? "4 " : "0 ") << endl;
			//cout << " -- Pop Mask: " << (popMask & 1 ? "1 " : "0 ") << (popMask & 2 ? "2 " : "0 ") << (popMask & 4 ? "4 " : "0 ") << endl;
			entry->octantMask ^= stepMask;
			//cout << " -- New Octant Mask: " << (entry->octantMask & 1 ? "1 " : "0 ") << (entry->octantMask & 2 ? "2 " : "0 ") << (entry->octantMask & 4 ? "4 " : "0 ") << endl;

			if (!popMask) break;

			//cout << " - Pop" << endl;

			// If we exit the top level we done
			if (sp == 0 && popMask) return result;

			// Pop of the stack
			entry = &stck[--sp];
			//cout << " -- " << sp << ": " << entry->id << " (" << (int)entry->scale << ", <" << x << ", " << y << ", " << z<< ">)" << endl;
		}
	}

	return result;
}

VoxelGrid::Result VoxelGrid::intersectOctreeFast(Ray ray, float maxt, uint8_t materialType) const
{
	Result result;
	result.normal = vec3(1.0f, 1.0f, 1.0f);
	result.hit = false;

	cout << endl << "========================" << endl << endl;

	int temp = size.x;
	int scale = 0;

	while (temp >>= 1) ++scale; 

	int hscale = scale - 1;

	float tmin = 0.0f;
	float tmax = numeric_limits<float>::max();//maxt;

	vec3 p = ray.pos;
	vec3 lp = p / float(1 << scale);

	int x = p.x, y = p.y, z = p.z;
	x = clamp(x, 0, size.x);
	y = clamp(y, 0, size.y);
	z = clamp(z, 0, size.z);

	uint8_t	rayMask = 0;
	if (ray.dir.x > 0) rayMask |= 1;
	if (ray.dir.y > 0) rayMask |= 2;
	if (ray.dir.z > 0) rayMask |= 4;

	struct Entry
	{
		uint32_t	id;
		uint8_t		octantMask;
		ivec3		off;
	};

	int sp = 0;
	Entry stck[32];
	{
		ivec3 off = ivec3(lp * 2.0f);

		uint8_t	octantMask = off.x + off.y * 2 + off.z * 4;

		stck[0].id = materialType == 0 ? 0U : waterNodeIndex;
		stck[0].octantMask = octantMask;
		stck[0].off = off;

		cout << "Off: " << off.x << ", " << off.y << ", " << off.z << endl;
		cout << "LP: " << lp.x << ", " << lp.y << ", " << lp.z << endl;
	}

	//cout << "Ray Mask: " << (rayMask & 1 ? "1 " : "0 ") << (rayMask & 2 ? "2 " : "0 ") << (rayMask & 4 ? "4 " : "0 ") << endl;

	int side = 0;
	iters = 0;

	Entry* entry = &stck[sp];
	while (sp < 32)
	{
		iters++;

		Node& node = nodes[entry->id];

		cout << endl;
		for (int k = 0; k < sp; k++)
		{
			Entry* e = &stck[k];
			cout << k << ": " << e->id << endl;
			cout << " - Octant Mask: " << (e->octantMask & 1 ? "1 " : "0 ") << (e->octantMask & 2 ? "2 " : "0 ") << (e->octantMask & 4 ? "4 " : "0 ") << endl;
		}

		cout << sp << ": " << entry->id << " (" << (1 << (int)scale) << ", <" << x << ", " << y << ", " << z<< ">)" << endl;
		cout << " - Octant Mask: " << (entry->octantMask & 1 ? "1 " : "0 ") << (entry->octantMask & 2 ? "2 " : "0 ") << (entry->octantMask & 4 ? "4 " : "0 ") << endl;
		cout << " - Value: " << get(
			x + (entry->octantMask & 1 ? 1: 0), 
			y + (entry->octantMask & 2 ? 1: 0), 
			z + (entry->octantMask & 4 ? 1: 0)
		) << endl;

		// Check if the current octant has child
		int8_t childID = 1U << entry->octantMask;
		if (node.child & childID)
		{
			// Check if it is a leaf node
			if (node.leaf & childID) 
			{
				result.hit = true;
				result.t = tmin;
				result.pos = ray.pos + ray.dir * (tmin + 1e-3f);
				result.pos.x = std::min(std::max((int)result.pos.x, 0), size.x);
				result.pos.y = std::min(std::max((int)result.pos.y, 0), size.y);
				result.pos.z = std::min(std::max((int)result.pos.z, 0), size.z);

				result.side = side;
				result.sign = ((rayMask >> side) & 1) ? 1.0f : -1.0f;
				result.normal = normals[side] * result.sign;

				return result;
			}

			int hhscale = hscale - 1;

			cout << " - Push - " << (1<<hhscale) << endl;
			cout << " - " << x << ", " << y << ", " << z << endl;

			lp = lp * 2.0f - vec3(entry->off);
			ivec3 off = ivec3(lp * 2.0f);

			cout << " - Off: " << off.x << ", " << off.y << ", " << off.z << endl;
			cout << " - LP: " << lp.x << ", " << lp.y << ", " << lp.z << endl;

			//// Push
			uint8_t	octantMask = 0;
			if (x & (1<<hhscale)) octantMask |= 1;
			if (y & (1<<hhscale)) octantMask |= 2;
			if (z & (1<<hhscale)) octantMask |= 4;
			cout << " - Octant Mask: " << (octantMask & 1 ? "1 " : "0 ") << (octantMask & 2 ? "2 " : "0 ") << (octantMask & 4 ? "4 " : "0 ") << endl;

			// Decode pointer
			uint32_t children = node.children;
			if (node.children & 0x8000)
				children = ptrTable[node.children & 0x7FFF];

			// Add to stack
			++sp;
			stck[sp].id = entry->id + children + entry->octantMask;
			stck[sp].octantMask = octantMask;
			stck[sp].off = off;
			entry = &stck[sp];

			scale = hscale;
			hscale = hhscale;

			continue;
		}

		//cout << " - Advance" << endl;

		//// Advance
		int nx = ((x >> hscale) << hscale) + (rayMask & 1 ? (1<<hscale) : 0);
		int ny = ((y >> hscale) << hscale) + (rayMask & 2 ? (1<<hscale) : 0);
		int nz = ((z >> hscale) << hscale) + (rayMask & 4 ? (1<<hscale) : 0);

		//cout << " -- Scale: " << hscale << endl;
		//cout << " -- Pos: <" << x << ", " << y << ", " << z << ">" << endl;
		//cout << " -- Next: <" << nx << ", " << ny << ", " << nz << ">" << endl;

		vec3 nextt = (vec3(nx, ny, nz) - ray.pos) * ray.invDir;
		//cout << " -- NextT: <" << nextt.x << ", " << nextt.y << ", " << nextt.z << ">" << endl;

		tmin = std::max(std::min(std::min(nextt.x, nextt.y), nextt.z), tmin);

		uint8_t stepMask = 0;
		if (nextt.x <= tmin) {stepMask |= 1; side = 0;}
		if (nextt.y <= tmin) {stepMask |= 2; side = 1;}
		if (nextt.z <= tmin) {stepMask |= 4; side = 2;}

		if (tmin >= tmax) break;

		//tmin += 1e-4f;

		//cout << " -- Step Mask: " << (stepMask & 1 ? "1 " : "0 ") << (stepMask & 2 ? "2 " : "0 ") << (stepMask & 4 ? "4 " : "0 ") << endl;

		vec3 cp = ray.pos + ray.dir * tmin;
		x = clamp((int)(cp.x), 0, size.x);
		y = clamp((int)(cp.y), 0, size.y);
		z = clamp((int)(cp.z), 0, size.z);

		if (nextt.x <= tmin) x = nx - (rayMask & 1 ? 0 : 1);
		if (nextt.y <= tmin) y = ny - (rayMask & 2 ? 0 : 1);
		if (nextt.z <= tmin) z = nz - (rayMask & 4 ? 0 : 1);

		lp = cp / float(1 << scale);

		//// Pop
		while (sp >= 0)
		{
			uint8_t endMask = ~(entry->octantMask ^ rayMask);
			uint8_t popMask = endMask & stepMask;
			//cout << " -- Octant Mask: " << (entry->octantMask & 1 ? "1 " : "0 ") << (entry->octantMask & 2 ? "2 " : "0 ") << (entry->octantMask & 4 ? "4 " : "0 ") << endl;
			//cout << " -- End Mask: " << (endMask & 1 ? "1 " : "0 ") << (endMask & 2 ? "2 " : "0 ") << (endMask & 4 ? "4 " : "0 ") << endl;
			//cout << " -- Pop Mask: " << (popMask & 1 ? "1 " : "0 ") << (popMask & 2 ? "2 " : "0 ") << (popMask & 4 ? "4 " : "0 ") << endl;
			entry->octantMask ^= stepMask;
			//cout << " -- New Octant Mask: " << (entry->octantMask & 1 ? "1 " : "0 ") << (entry->octantMask & 2 ? "2 " : "0 ") << (entry->octantMask & 4 ? "4 " : "0 ") << endl;

			if (!popMask) break;

			//cout << " - Pop" << endl;

			// If we exit the top level we done
			if (sp == 0 && popMask) return result;

			// Pop of the stack
			entry = &stck[--sp];

			lp = (lp + vec3(entry->off)) * 0.5f;

			hscale = scale;
			scale = hscale + 1;

			//cout << " -- " << sp << ": " << entry->id << " (" << (int)entry->scale << ", <" << x << ", " << y << ", " << z<< ">)" << endl;
		}
	}

	return result;
}

bool VoxelGrid::intersectShadow(Ray ray, float maxt) const
{
	float	t = 0.0f;
	vec3 p = ray.pos + ray.dir * t;
	int x = p.x, y = p.y, z = p.z;
	x = clamp(x, 0, size.x);
	y = clamp(y, 0, size.y);
	z = clamp(z, 0, size.z);

	vec3 dsign = vec3(
		sign(ray.dir.x),
		sign(ray.dir.y),
		sign(ray.dir.z)
	);

	vec3 next;
	next.x = 0.0f + (ray.dir.x >= 0 ? x + 1 : x);
	next.y = 0.0f + (ray.dir.y >= 0 ? y + 1 : y);
	next.z = 0.0f + (ray.dir.z >= 0 ? z + 1 : z);

	vec3 tmax = (next - p) * ray.invDir, tdelta = dsign * ray.invDir;

	for (int l = 0; l < 1000; l++)
	{
		int cell = get(x, y, z);

		if (cell != 0) return true;

		if (tmax.x < tmax.y)
		{
			if (tmax.x < tmax.z)
			{
				x = x + dsign.x;
				if (tmax.x > maxt || x < 0 || x >= size.x) break;
				tmax.x += tdelta.x;
			}
			else
			{
				z = z + dsign.z;
				if (tmax.z > maxt || z < 0 || z >= size.z) break;
				tmax.z += tdelta.z;
			}
		}
		else
		{
			if (tmax.y < tmax.z)
			{
				y = y + dsign.y;
				if (tmax.y > maxt || y < 0 || y >= size.y) break;
				tmax.y += tdelta.y;
			}
			else
			{
				z = z + dsign.z;
				if (tmax.z > maxt || z < 0 || z >= size.z) break;
				tmax.z += tdelta.z;
			}
		}
	}

	return false;
}

Voxels::Result VoxelGrid::intersect(Ray ray, float maxt, uint8_t materialType) const
{
	Result result;
	result.normal = vec3(1.0f, 1.0f, 1.0f);
	result.hit = false;

	IntersectResult bresult = aabb.intersect(ray, 0.0f);
	if (!bresult.hit) return result;

	float	t = bresult.tmin;
	vec3 p = ray.pos + ray.dir * t;
	int x = p.x, y = p.y, z = p.z;
	x = clamp(x, 0, size.x);
	y = clamp(y, 0, size.y);
	z = clamp(z, 0, size.z);

	vec3 dsign = vec3(
		sign(ray.dir.x),
		sign(ray.dir.y),
		sign(ray.dir.z)
	);

	vec3 next;
	next.x = 0.0f + (ray.dir.x >= 0 ? x + 1 : x);
	next.y = 0.0f + (ray.dir.y >= 0 ? y + 1 : y);
	next.z = 0.0f + (ray.dir.z >= 0 ? z + 1 : z);

	vec3 tmax = (next - p) * ray.invDir, tdelta = dsign * ray.invDir;

	iters = 0;
	int side = 0;

	for (int l = 0; l < 1000; l++)
	{
		iters++;

		int cell = get(x, y, z);

		if (cell != materialType)
		{
			result.pos = vec3(x, y, z);
			result.hit = true;
			result.side = side;
			result.sign = dsign[side];
			result.normal = normals[side] * dsign;
			return result;
		}

		if (tmax.x < tmax.y)
		{
			if (tmax.x < tmax.z)
			{
				x = x + dsign.x;
				result.t = tmax.x;
				tmax.x += tdelta.x;
				side = 0;
				if (result.t > maxt || x < 0 || x >= size.x) break;
			}
			else
			{
				z = z + dsign.z;
				result.t = tmax.z;
				tmax.z += tdelta.z;
				side = 2;
				if (result.t > maxt || z < 0 || z >= size.z) break;
			}
		}
		else
		{
			if (tmax.y < tmax.z)
			{
				y = y + dsign.y;
				result.t = tmax.y;
				tmax.y += tdelta.y;
				side = 1;
				if (result.t > maxt || y < 0 || y >= size.y) break;
			}
			else
			{
				z = z + dsign.z;
				result.t = tmax.z;
				tmax.z += tdelta.z;
				side = 2;
				if (result.t > maxt || z < 0 || z >= size.z) break;
			}
		}
	}

	return result;
}

void VoxelGrid::generateOctree()
{
	space = 1;
	ptrSpace = 0;

	// Create octree for air
	doNode(0U, 0, 0, 0, size.x, 0U);

	cout << "Allocated: " << space << endl;
	cout << "Allocated Ptrs: " << ptrSpace << endl;

	// Create octree for water
	//waterNodeIndex = space;
	//doNode(waterNodeIndex, 0, 0, 0, size.x, 4U);

	//cout << "[W] Allocated: " << space << endl;
	//cout << "[W] Allocated Ptrs: " << ptrSpace << endl;

	cout << "Octree space: " << sizeof(Node)*space << " bytes" << endl;
	cout << "Pointer space: " << sizeof(uint32_t)*ptrSpace << " bytes" << endl;
	cout << "Voxels space: " << sizeof(uint8_t)*size.x*size.y*size.z << " bytes" << endl;
}

uint8_t VoxelGrid::doNode(uint32_t id, int x, int y, int z, int scale, uint8_t type)
{
	if (scale == 1)
	{
		int cell = get(x, y, z);
		return cell == type ? 0 : 2;
	}

	Node& self = nodes[id];
	self.children = 0U;
	self.leaf = 0U;
	self.child = 0U;

	uint32_t children = space - id;

	if (children > 0x7FFF)
	{
		uint16_t ptrIndex = ptrSpace++;
		ptrTable[ptrIndex] = children;

		self.children = ptrIndex | 0x8000;
	}
	else
		self.children = children;

	space += 8;

	bool empty = true;
	bool full = true;
	for (uint8_t i = 0; i < 8; i++)
	{
		const int hscale = scale >> 1;
		const int cx = x + (i & 1 ? hscale : 0);
		const int cy = y + (i & 2 ? hscale : 0);
		const int cz = z + (i & 4 ? hscale : 0);

		const uint8_t res = doNode(id + children + i, cx, cy, cz, hscale, type);
		if (res != 2) full = false;
		if (res)
		{
			empty = false;
			self.child |= (1 << i);

			if (res == 2) self.leaf |= (1 << i);
		}
	}

	if (empty || full) space -= 8;

	return full ? 2 : (empty ? 0 : 1);
}

void VoxelGrid::remake()
{
	for (int x = 0; x < size.x; x++)
		for (int y = 0; y <- size.y; y++)
			for (int z = 0; z < size.z; z++)
				get(x, y, z) = 0;

	test(0U, 0, 0, 0, size.x);
}

void VoxelGrid::test(uint32_t id, int x, int y, int z, int scale)
{	
	Node& self = nodes[id];

	for (int8_t i = 0; i < 8; i++)
	{
		const int hscale = scale >> 1;
		const int cx = x + (i & 1 ? hscale : 0);
		const int cy = y + (i & 2 ? hscale : 0);
		const int cz = z + (i & 4 ? hscale : 0);

		uint8_t mask = 1 << i;
		if ((self.child & mask) != 0)
		{
			if ((self.leaf & mask) != 0)
			{
				for (int kx = 0; kx < hscale; kx++)
					for (int ky = 0; ky < hscale; ky++)
						for (int kz = 0; kz < hscale; kz++)
							get(cx+kx, cy+ky, cz+kz) = 1;
			}
			else
				test(id + self.children + i, cx, cy, cz, hscale);
		}
	}
}