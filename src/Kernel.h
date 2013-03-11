#ifndef _Kernel_h
#define _Kernel_h

#include <string>
#include <iostream>
#include <cstring>

#define CL_USE_DEPRECATED_OPENCL_1_1_APIS

#include <CL/cl.h>
#include <CL/cl_gl.h>

using namespace std;

class Kernel
{
	public:
		Kernel(string name, string filename, string mainFunc, cl_context context, cl_device_id device);
		~Kernel();

		bool build();

		bool	isValid() {return valid;}
		string	getName() {return name;}

	private:
		string name;
		string filename;
		string mainFunc;

		bool valid;
		cl_context context;
		cl_device_id device;
		cl_program program;

	public:
		cl_kernel kernel;
};

#endif