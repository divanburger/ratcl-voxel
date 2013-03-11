#ifndef _Math_h
#define _Math_h

namespace math
{
	inline float mix(float a, float b, float f)
	{
		return a + (b - a) * f;
	}
}

inline float smoothmix(float a, float b, float f)
{
	return math::mix(a, b, f * f * (3.0f - 2.0f * f));
}

#endif
