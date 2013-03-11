/*
 * AABB.cpp
 *
 *  Created on: 26 Oct 2011
 *      Author: dburger
 */

#include "AABB.h"

AABB::AABB(vec3 min, vec3 max) : min(min), max(max)
{

}

IntersectResult AABB::intersect(Ray ray, float start) const
{
	IntersectResult result;
	result.hit = false;

	float tmin, tmax, tmin2, tmax2;

	if (ray.invDir.x >= 0)
	{
		tmin = (min.x - ray.pos.x) * ray.invDir.x;
		tmax = (max.x - ray.pos.x) * ray.invDir.x;
	}
	else
	{
		tmin = (max.x - ray.pos.x) * ray.invDir.x;
		tmax = (min.x - ray.pos.x) * ray.invDir.x;
	}

	if (ray.invDir.y >= 0)
	{
		tmin2 = (min.y - ray.pos.y) * ray.invDir.y;
		tmax2 = (max.y - ray.pos.y) * ray.invDir.y;
	}
	else
	{
		tmin2 = (max.y - ray.pos.y) * ray.invDir.y;
		tmax2 = (min.y - ray.pos.y) * ray.invDir.y;
	}

	if ((tmin > tmax2) || (tmin2 > tmax)) return result;

	tmin = (tmin2 > tmin) ? tmin2 : tmin;
	tmax = (tmax2 < tmax) ? tmax2 : tmax;

	if (ray.invDir.z > 0)
	{
		tmin2 = (min.z - ray.pos.z) * ray.invDir.z;
		tmax2 = (max.z - ray.pos.z) * ray.invDir.z;
	}
	else
	{
		tmin2 = (max.z - ray.pos.z) * ray.invDir.z;
		tmax2 = (min.z - ray.pos.z) * ray.invDir.z;
	}

	if ((tmin > tmax2) || (tmin2 > tmax)) return result;

	tmin = (tmin2 > tmin) ? tmin2 : tmin;
	tmax = (tmax2 < tmax) ? tmax2 : tmax;

	result.hit = (tmax > start);

	if (result.hit)
	{
		result.tmin = (tmin < start) ? start : tmin;
		result.tmax = tmax;
	}

	return result;
}
