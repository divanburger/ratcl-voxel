/*
 * Ray.h
 *
 *  Created on: 26 Oct 2011
 *      Author: dburger
 */

#ifndef RAY_H_
#define RAY_H_

#include <glm/glm.hpp>

using namespace glm;

struct Ray
{
	vec3	pos;
	vec3	dir;
	vec3 invDir;

	Ray() : pos(), dir(0.0f, 0.0f, 1.0f) {}
	Ray(vec3 pos, vec3 dir) : pos(pos), dir(dir) {calcDirInv();}

	void calcDirInv()
	{
		invDir = vec3(
			fabs(dir.x) >= 1e-6f ? (1.0f / dir.x) : 1e-6f * sign(dir.x),
			fabs(dir.y) >= 1e-6f ? (1.0f / dir.y) : 1e-6f * sign(dir.y),
			fabs(dir.z) >= 1e-6f ? (1.0f / dir.z) : 1e-6f * sign(dir.z)
		);
	}
};

#endif /* RAY_H_ */
