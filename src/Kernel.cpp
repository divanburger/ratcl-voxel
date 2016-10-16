#include "Kernel.h"

#define MAX_SOURCE_SIZE (0x100000)

Kernel::Kernel(string name, string filename, string mainFunc, cl_context context, cl_device_id device)
	: name(name), filename(filename), mainFunc(mainFunc), valid(false),
	context(context), device(device)
{
	FILE *fp;
	char *source_str;
	size_t source_size;
	fp = fopen(("src/" + filename).c_str(), "r");
	if (!fp)
	{
		cerr << "Failed to load kernel: " << filename << endl;
		exit(1);
	}

	source_str = (char*)malloc(MAX_SOURCE_SIZE);
	source_size = fread(source_str, 1, MAX_SOURCE_SIZE, fp);
	fclose(fp);

	// Create a program from the kernel source
	cl_int ret;
	program = clCreateProgramWithSource(context, 1,
					(const char **)&source_str, (const size_t *)&source_size, &ret);

	free(source_str);
}

Kernel::~Kernel()
{
	if (valid)
	{
		clReleaseKernel(kernel);
		clReleaseProgram(program);
	}
}

bool Kernel::build()
{
	// Build the program
	cl_int ret = clBuildProgram(program, 1, &device, "-cl-nv-verbose -cl-fast-relaxed-math -cl-no-signed-zeros -cl-mad-enable", NULL, NULL);

	if (ret != CL_SUCCESS)
		cerr << "Build failed: " << filename << endl;
	else
		cout << "Build success: " << filename << endl;

	size_t retValSize;
	ret = clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &retValSize);
    if (ret != CL_SUCCESS)
    {
        cerr << "Error getting build log: " << filename << endl;
        return false;
    }

	char *buildLog = (char *)malloc(retValSize + 1);
	ret = clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, retValSize, buildLog, NULL);

    if (retValSize > 2) printf("Build Log:\n\n%s", buildLog);

	if (ret == CL_SUCCESS)
	{
		valid = true;
		kernel = clCreateKernel(program, mainFunc.c_str(), &ret);
	}

	return ret == CL_SUCCESS;
}