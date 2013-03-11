#ifndef _Util_h
#define _Util_h

#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <cstring>
#include <iomanip>

using namespace std;

template<typename T>
string toString(T value)
{
	std::stringstream ss;
	ss << value;
	return ss.str();
}

template<typename T>
string toString(T value, int precision)
{
	std::stringstream ss;
	ss << setprecision(precision) << value;
	return ss.str();
}

template<>
inline string toString<bool>(bool value)
{
	return value ? "true" : "false";
}

template<typename T>
T toNumber(const string& str, int fallback = 0)
{
	T result;
	return (std::stringstream(str) >> result) ? result : fallback;
}

#endif