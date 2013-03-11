/*
 * Voxels.h
 *
 *  Created on: 6 Nov 2011
 *      Author: dburger
 */

#ifndef VOXELS_H_
#define VOXELS_H_

#include <glm/glm.hpp>

using namespace glm;

class Voxels
{
	public:
		struct Result
		{
			bool 	hit;
			float	t;
			vec3	pos;
			int		side;
			float	sign;
			vec3	normal;
		};
};

#endif /* VOXELS_H_ */
