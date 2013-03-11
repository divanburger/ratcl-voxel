#ifndef _Vector_h
#define _Vector_h

#include <cmath>
#include <ostream>

using namespace std;

struct Vector2i
{
	int x, y;

	Vector2i() : x(0), y(0) {}
	Vector2i(const Vector2i& v) : x(v.x), y(v.y) {}
	explicit Vector2i(int x, int y) : x(x), y(y) {}
};

ostream& operator<< (ostream& o, Vector2i v);

struct Vector2f
{
	float x, y;

	Vector2f() : x(0.0f), y(0.0f) {}
	Vector2f(const Vector2f& v) : x(v.x), y(v.y) {}
	explicit Vector2f(float x, float y) : x(x), y(y) {}

	Vector2f operator-() {return Vector2f(-x, -y);}

	Vector2f operator+(Vector2f v) {return Vector2f(x+v.x, y+v.y);}
	Vector2f operator-(Vector2f v) {return Vector2f(x-v.x, y-v.y);}
	Vector2f operator*(Vector2f v) {return Vector2f(x*v.x, y*v.y);}

	Vector2f operator+(float f) {return Vector2f(x+f, y+f);}
	Vector2f operator-(float f) {return Vector2f(x-f, y-f);}
	Vector2f operator*(float f) {return Vector2f(x*f, y*f);}

	Vector2f operator+=(Vector2f v) {x+=v.x; y+=v.y; return *this;}
	Vector2f operator-=(Vector2f v) {x-=v.x; y-=v.y; return *this;}
	Vector2f operator*=(Vector2f v) {x*=v.x; y*=v.y; return *this;}

	float		length() const {return sqrt(x*x + y*y);}
	Vector2f	normalise() const {float ilen = 1.0f / length(); return Vector2f(x*ilen, y*ilen);}
	float		dot(Vector2f v) const {return x*v.x + y*v.y;}
	float		distance(Vector2f to) const {return (to - *this).length();}

	Vector2f	rotate(float angle) {float cosa = cos(angle), sina = sin(angle); return Vector2f(x * cosa + y * sina, y * cosa - x * sina);}

	Vector2f	perpendicular() {return Vector2f(y, -x);}
};

typedef Vector2f float2;

inline Vector2f abs(Vector2f v) {return Vector2f(fabs(v.x), fabs(v.y));}
inline Vector2f reflect(Vector2f n, Vector2f i) {return i - n * n.dot(i) * 2;}

ostream& operator<< (ostream& o, Vector2f v);

struct Vector3i
{
	union
	{
		int v[3];
		struct {int x, y, z;};
	};

	Vector3i() : x(0), y(0), z(0) {}
	Vector3i(const Vector3i& v) : x(v.x), y(v.y), z(v.z) {}
	Vector3i(const int* v) : x(v[0]), y(v[1]), z(v[2]) {}
	Vector3i(int x, int y = 0, int z = 0) : x(x), y(y), z(z) {}

	Vector3i operator-() const {return Vector3i(-x, -y, -z);}

	Vector3i operator+(Vector3i v) const {return Vector3i(x+v.x, y+v.y, z+v.z);}
	Vector3i operator-(Vector3i v) const {return Vector3i(x-v.x, y-v.y, z-v.z);}
	Vector3i operator*(Vector3i v) const {return Vector3i(x*v.x, y*v.y, z*v.z);}

	Vector3i operator+=(Vector3i v) {x+=v.x; y+=v.y; z+=v.z; return *this;}
	Vector3i operator-=(Vector3i v) {x-=v.x; y-=v.y; z-=v.z; return *this;}
	Vector3i operator*=(Vector3i v) {x*=v.x; y*=v.y; z*=v.z; return *this;}

	Vector3i operator+(int v) const {return Vector3i(x+v, y+v, z+v);}
	Vector3i operator-(int v) const {return Vector3i(x-v, y-v, z-v);}
	Vector3i operator*(int v) const {return Vector3i(x*v, y*v, z*v);}

	Vector3i operator>>(int v) const {return Vector3i(x>>v, y>>v, z>>v);}
	Vector3i operator<<(int v) const {return Vector3i(x<<v, y<<v, z<<v);}
};

struct Vector3f
{
	union
	{
		float v[3];
		struct {float x, y, z;};
		struct {float r, g, b;};
	};

	Vector3f() : x(0.0f), y(0.0f), z(0.0f) {}
	Vector3f(const Vector3i& v) : x(v.x), y(v.y), z(v.z) {}
	Vector3f(const Vector3f& v) : x(v.x), y(v.y), z(v.z) {}
	Vector3f(const float* v) : x(v[0]), y(v[1]), z(v[2]) {}
	explicit Vector3f(const Vector2f& v, float z) : x(v.x), y(v.y), z(z) {}
	Vector3f(float x, float y = 0.0f, float z = 0.0f) : x(x), y(y), z(z) {}

	operator Vector3i() {return Vector3i(x, y, z);}

	Vector3f operator-() const {return Vector3f(-x, -y, -z);}

	Vector3f operator+(Vector3f v) const {return Vector3f(x+v.x, y+v.y, z+v.z);}
	Vector3f operator-(Vector3f v) const {return Vector3f(x-v.x, y-v.y, z-v.z);}
	Vector3f operator*(Vector3f v) const {return Vector3f(x*v.x, y*v.y, z*v.z);}
	Vector3f operator/(Vector3f v) const {return Vector3f(x/v.x, y/v.y, z/v.z);}

	Vector3f operator+=(Vector3f v) {x+=v.x; y+=v.y; z+=v.z; return *this;}
	Vector3f operator-=(Vector3f v) {x-=v.x; y-=v.y; z-=v.z; return *this;}
	Vector3f operator*=(Vector3f v) {x*=v.x; y*=v.y; z*=v.z; return *this;}

	Vector3f operator+(float v) const {return Vector3f(x+v, y+v, z+v);}
	Vector3f operator-(float v) const {return Vector3f(x-v, y-v, z-v);}
	Vector3f operator*(float v) const {return Vector3f(x*v, y*v, z*v);}

	Vector3f operator+=(float v) {return Vector3f(x+=v, y+=v, z+=v);}
	Vector3f operator-=(float v) {return Vector3f(x-=v, y-=v, z-=v);}
	Vector3f operator*=(float v) {return Vector3f(x*=v, y*=v, z*=v);}

	float		length() const {return sqrt(x*x + y*y + z*z);}
	Vector3f	normalise() const {float ilen = 1.0f / length(); return Vector3f(x*ilen, y*ilen, z*ilen);}

	float		dot(Vector3f v) const {return x*v.x + y*v.y + z*v.z;}
	Vector3f	cross(Vector3f v) const {return Vector3f(y*v.z-z*v.y, z*v.x-x*v.z, x*v.y-y*v.x);}
	float		distance(Vector3f to) const {return (to - *this).length();}
};

typedef Vector3f float3;

inline Vector3f operator/(float s, Vector3f v) {return Vector3f(s/v.x, s/v.y, s/v.z);}

inline float	maxComponent(Vector3f v) {return std::max(max(v.x, v.y), v.z);}
inline float	minComponent(Vector3f v) {return min(min(v.x, v.y), v.z);}
inline Vector3f abs(Vector3f v) {return Vector3f(fabs(v.x), fabs(v.y), fabs(v.z));}
inline Vector3f floor(Vector3f v) {return Vector3f(floor(v.x), floor(v.y), floor(v.z));}
inline Vector3f pow(Vector3f v, float f) {return Vector3f(pow(v.x, f), pow(v.y, f), pow(v.z, f));}
inline Vector3f mix(Vector3f a, Vector3f b, float f) {return a + (b - a) * f;}
inline Vector3f reflect(Vector3f normal, Vector3f in) {return in - normal * normal.dot(in) * 2;}

inline Vector3f headingAndPitchToVector(float heading, float pitch)
{
	float	sinp = sin(pitch);
	return Vector3f(cos(heading) * sinp, cos(pitch), sin(heading) * sinp);
}

inline Vector3f headingToVectorY(float heading) {return Vector3f(cos(heading), 0.0f, sin(heading));}

ostream& operator<< (ostream& o, Vector3f v);

struct Vector4f
{
	union
	{
		float v[4];
		float xyz[3];
		struct {float x, y, z, w;};
		struct {float r, g, b, a;};
	};

	Vector4f() : x(0.0f), y(0.0f), z(0.0f), w(0.0f) {}
	Vector4f(const Vector4f& v) : x(v.x), y(v.y), z(v.z), w(v.w) {}
	explicit Vector4f(const Vector3f& v, float w = 1.0f) : x(v.x), y(v.y), z(v.z), w(w) {}
	Vector4f(float x, float y = 0.0f, float z = 0.0f, float w = 1.0f) : x(x), y(y), z(z), w(w) {}

	operator Vector3f() const {return Vector3f(x, y, z);}

	Vector4f operator+(Vector4f v) const {return Vector4f(x+v.x, y+v.y, z+v.z, w);}
	Vector4f operator-(Vector4f v) const {return Vector4f(x-v.x, y-v.y, z-v.z, w);}
	Vector4f operator*(Vector4f v) const {return Vector4f(x*v.x, y*v.y, z*v.z, w*v.w);}

	Vector4f operator+(float v) const {return Vector4f(x+v, y+v, z+v, w);}
	Vector4f operator-(float v) const {return Vector4f(x-v, y-v, z-v, w);}
	Vector4f operator*(float v) const {return Vector4f(x*v, y*v, z*v, w);}
};

inline Vector4f pow(Vector4f v, float f) {return Vector4f(pow(v.x, f), pow(v.y, f), pow(v.z, f), v.w);}

ostream& operator<< (ostream& o, Vector4f v);


typedef Vector2f Vector2;
typedef Vector3f Vector3;
typedef Vector4f Vector4;

#endif
