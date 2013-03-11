/*
 * VoxelGrid.h
 *
 *  Created on: 6 Nov 2011
 *      Author: dburger
 */

#ifndef VOXELGRID_H_
#define VOXELGRID_H_

#include <memory>
#include <stack>
#include <iostream>
#include <limits>
#include <glm/glm.hpp>

#include "Ray.h"
#include "AABB.h"
#include "Voxels.h"

using namespace std;
using namespace glm;

struct Node
{
	uint16_t	children;
	uint8_t		child;
	uint8_t		leaf;
};

class VoxelGrid : public Voxels
{
	public:
		VoxelGrid(ivec3 size);
		virtual ~VoxelGrid();

		inline uint8_t	get(int x, int y, int z) const {return voxels[x + (y + z * size.y) * size.x];}
		inline uint8_t&	get(int x, int y, int z) {return voxels[x + (y + z * size.y) * size.x];}

		uint8_t*	getVoxels() {return voxels;}

		bool 	intersectShadowOctree(Ray ray, float maxt) const;
		Result 	intersectOctree(Ray ray, float maxt, uint8_t materialType = 0) const;
		Result 	intersectOctreeFast(Ray ray, float maxt, uint8_t materialType = 0) const;

		bool 	intersectShadow(Ray ray, float maxt) const;
		Result 	intersect(Ray ray, float maxt, uint8_t materialType = 0) const;

		void	generateOctree();
		void	remake();

		ivec3		size;
		uint8_t*	voxels;
		AABB		aabb;

		Node*		nodes;
		uint32_t*	ptrTable;

		uint32_t	space;
		uint16_t	ptrSpace;

		uint32_t	waterNodeIndex;

		mutable int iters;

	private:
		uint8_t doNode(uint32_t id, int x, int y, int z, int scale, uint8_t type);

		void test(uint32_t id, int x, int y, int z, int scale);

		static const vec3 normals[3];
};

#endif /* VOXELGRID_H_ */
