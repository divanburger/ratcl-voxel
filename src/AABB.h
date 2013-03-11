/*
 * AABB.h
 *
 *  Created on: 26 Oct 2011
 *      Author: dburger
 */

#ifndef AABB_H_
#define AABB_H_

#include <glm/glm.hpp>

#include "Ray.h"

using namespace glm;

struct IntersectResult
{
	bool hit;
	float tmin;
	float tmax;
};

class AABB
{
	public:
		AABB(vec3 min, vec3 max);
		virtual ~AABB() {}

		IntersectResult intersect(Ray ray, float start) const;

		vec3	min;
		vec3 max;
};

#endif /* AABB_H_ */
