#include <SDL/SDL.h>
#include <SDL/SDL_ttf.h>
#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/glx.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

#include <stdio.h>
#include <stdlib.h>

#define CL_USE_DEPRECATED_OPENCL_1_1_APIS

#include <CL/cl.h>
#include <CL/cl_gl.h>

#include <glm/glm.hpp>

#include "Scene.h"
#include "Raytracer.h"
#include "Noise.h"
#include "Voxels.h"
#include "VoxelGrid.h"
#include "Ray.h"
#include "Camera.h"
#include "Texture.h"
#include "GUI.h"
#include "Kernel.h"

using namespace std;

const int WIDTH = 640;
const int HEIGHT = 480;
const int HWIDTH = WIDTH >> 1;
const int HHEIGHT = HEIGHT >> 1;
const bool FULLSCREEN = false;

const float EPSILON = 1e-4;

const ivec3	voxelSize(128, 128, 128);
const int	voxelCount = voxelSize.x * voxelSize.y * voxelSize.z;
const int	seaLevel = 64;
VoxelGrid	voxels(voxelSize);
bool		voxelsDirty = true;

bool key[256] = {false};

vec3 tangents[3] = {vec3(0.0f, 1.0f, 0.0f), vec3(0.0f, 0.0f, 1.0f), vec3(1.0f, 0.0f, 0.0f)};

Camera	camera(vec3(8, 0, 8.), 7.01199, -0.292656);

vec3 velocity = vec3(0.0f, 0.0f, 0.0f);
vec3 acceleration = vec3(0.0f, 0.0f, 0.0);

float 	timeOfDay = 8.0f;
float 	exposure = 0.1f;
Scene	scene;

bool	opencl = true;
bool	grab = false;

Raytracer*	raytracer;
GLuint		outputTexture;
GLuint		screenTexture[2];
GLuint		finalTexture;

Texture*	grassTexture;

GUI*		gui;

void printCLErrorCode(string tag, cl_int code)
{
	if (code == 0) return;

	cout << tag << ": ";
	switch (code)
	{
		case CL_INVALID_VALUE: cout << "Invalid value"; break;
		case CL_INVALID_MEM_OBJECT: cout << "Invalid memory object"; break;
		case CL_INVALID_COMMAND_QUEUE: cout << "Invalid command queue"; break;
		case CL_INVALID_CONTEXT: cout << "Invalid OpenCL context"; break;
		case CL_INVALID_GL_OBJECT: cout << "Invalid OpenGL object"; break;
		case CL_INVALID_KERNEL_ARGS: cout << "Invalid Kernel args"; break;
		default: cout << "Unknown code " << code; break;
	}
	cout << endl;
}

bool setupOpenGL()
{
	GLenum glewError = glewInit();
	if (glewError != GLEW_OK)
	{
		cout << "GLEW Error: " << glewGetErrorString(glewError) << endl;
		return false;
	}

	glEnable(GL_FRAMEBUFFER_SRGB);

	return true;
}

GLuint createTexture()
{
	GLuint texture;
	glGenTextures(1, &texture);
	return texture;
}

void deleteBlock()
{
	Ray actionRay(camera.position, vec3(camera.getViewMatrix() * vec4(0.0f, 0.0f, 1.0f, 0.0f)));
	VoxelGrid::Result result = voxels.intersect(actionRay, 100.0f);
	if (result.hit)
	{
		cout << voxels.get(result.pos.x, result.pos.y, result.pos.z) << endl;
		voxels.get(result.pos.x, result.pos.y, result.pos.z) = 0;
	}
}

cl_platform_id selectPlatform()
{
	cl_uint num_platforms;
	cl_platform_id * platforms;
	clGetPlatformIDs(0, NULL, &num_platforms);
	platforms = new cl_platform_id[num_platforms];
	clGetPlatformIDs(num_platforms, platforms, NULL);

	for (unsigned int i = 0; i < num_platforms; i++)
	{
		cl_platform_id platform = platforms[i];

		size_t vendor_size;
		clGetPlatformInfo(platform, CL_PLATFORM_VENDOR, 0, NULL, &vendor_size);
		char* vendor = new char[vendor_size];
		clGetPlatformInfo(platform, CL_PLATFORM_VENDOR, vendor_size, vendor, NULL);

		size_t name_size;
		clGetPlatformInfo(platform, CL_PLATFORM_NAME, 0, NULL, &name_size);
		char* name = new char[name_size];
		clGetPlatformInfo(platform, CL_PLATFORM_NAME, name_size, name, NULL);

		size_t extensions_size;
		clGetPlatformInfo(platform, CL_PLATFORM_EXTENSIONS, 0, NULL, &extensions_size);
		char* extensions = new char[extensions_size];
		clGetPlatformInfo( platform, CL_PLATFORM_EXTENSIONS, extensions_size, extensions, NULL);

		cout << "Platform ID: " << i << endl;
		cout << "Platform vendor: " << vendor << endl;
		cout << "Platform name: " << name << endl;
		cout << "Platform extensions: " << extensions << endl;
		delete [] vendor;
		delete [] name;
		delete [] extensions;
	}

	cl_platform_id platform = platforms[0];

	delete [] platforms;

	return platform;
}

cl_device_id selectDevice(cl_platform_id platform)
{
	cl_uint num_devices;
	cl_device_id* devices;

	cl_int ret = clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, 0, NULL, &num_devices);
	printCLErrorCode("clGetDeviceIDs", ret);

	devices = new cl_device_id[num_devices];
	ret = clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, num_devices, devices, NULL);
	printCLErrorCode("clGetDeviceIDs", ret);

	for (unsigned int i = 0; i < num_devices; i++)
	{
		cl_device_id device = devices[i];

		size_t version_size;
		clGetDeviceInfo(device, CL_DEVICE_VERSION, 0, NULL, &version_size);
		char* version = new char[version_size];
		clGetDeviceInfo(device, CL_DEVICE_VERSION, version_size, version, NULL);

		cout << "Device ID: " << i << endl;
		cout << "Device version: " << version << endl;
		delete [] version;
	}

	cl_device_id device = devices[0];

	delete [] devices;

	return device;
}

void drawToScreen()
{
	// Draw to screen
	glBindTexture(GL_TEXTURE_2D, finalTexture);

	glLoadIdentity();
	glOrtho(0.0f, 1.0f, 1.0f, 0.0f, -1.0f, 1.0f);

	glEnable(GL_TEXTURE_2D);
	glColor3f(1.0f, 1.0f, 1.0f);
	glBegin(GL_TRIANGLES);
		glTexCoord2f(0.0f, 0.0f); glVertex2f(0.0f, 0.0f);
		glTexCoord2f(0.0f, 2.0f); glVertex2f(0.0f, 2.0f);
		glTexCoord2f(2.0f, 0.0f); glVertex2f(2.0f, 0.0f);
	glEnd();
	glBindTexture(GL_TEXTURE_2D, 0);
	glDisable(GL_TEXTURE_2D);

	if (grab)
	{
		float py = 0.01f;
		float px = py * HEIGHT / WIDTH;
		glColor3f(1.0f, 1.0f, 0.0f);
		glBegin(GL_QUADS);
			glVertex2f(0.5f, 0.5f);
			glVertex2f(0.5f, 0.5f+py);
			glVertex2f(0.5f+px, 0.5f+py);
			glVertex2f(0.5f+px, 0.5f);
		glEnd();
	}

	gui->render(!grab);

	SDL_GL_SwapBuffers();
}

void openCLLog(const char *msg, const void *, size_t, void *)
{
	cout << msg << endl;
	exit(1);
}

int main(void)
{
	// Setup scene details
	scene.fogDensity = 0.003f;
	scene.ambientColour = vec3(0.10f, 0.16f, 0.20f);
	scene.lightColour = vec3(0.63f, 0.55f, 0.5f);
	scene.lightDirection = normalize(vec3(1.0f, 0.25f, 0.5f));

	// Setup raytracer
	raytracer = new Raytracer(WIDTH, HEIGHT, scene, voxels, camera);

	// Setup SDL
	SDL_Init(SDL_INIT_EVERYTHING);
	SDL_SetVideoMode(WIDTH, HEIGHT, 32, SDL_OPENGL | (FULLSCREEN ? SDL_FULLSCREEN : 0));

	TTF_Init();

	// Setup screen
	setupOpenGL();

	gui = new GUI(ivec2(WIDTH, HEIGHT));

	// Get platform and device information
	cl_platform_id platform = selectPlatform();
	cl_device_id device = selectDevice(platform);

	// Create an OpenCL context
	cl_int ret;

    cl_context_properties properties[] =
    {
        CL_GL_CONTEXT_KHR, (cl_context_properties)glXGetCurrentContext(),
        CL_GLX_DISPLAY_KHR, (cl_context_properties)glXGetCurrentDisplay(),
        CL_CONTEXT_PLATFORM, (cl_context_properties)platform,
        0
    };

	cl_context context = clCreateContext(properties, 1, &device, openCLLog, NULL, &ret);
	printCLErrorCode("createContext", ret);

	// Create a command queue
	cl_command_queue command_queue = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &ret);

	// Create texture
	outputTexture = createTexture();
	screenTexture[0] = createTexture();
	screenTexture[1] = createTexture();
	finalTexture = createTexture();

	// Load textures
	grassTexture = new Texture("textures/grass.png");

	// Create pixel array
	int bytes = WIDTH * HEIGHT * sizeof(vec4);
	vec4 *pixels = (vec4*)malloc(bytes);

	glBindTexture(GL_TEXTURE_2D, outputTexture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, WIDTH, HEIGHT, 0, GL_RGBA, GL_FLOAT, NULL);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glBindTexture(GL_TEXTURE_2D, 0);

	glBindTexture(GL_TEXTURE_2D, screenTexture[0]);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, WIDTH, HEIGHT, 0, GL_RGBA, GL_FLOAT, NULL);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glBindTexture(GL_TEXTURE_2D, 0);

	glBindTexture(GL_TEXTURE_2D, screenTexture[1]);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, WIDTH, HEIGHT, 0, GL_RGBA, GL_FLOAT, NULL);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glBindTexture(GL_TEXTURE_2D, 0);

	glBindTexture(GL_TEXTURE_2D, finalTexture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, WIDTH, HEIGHT, 0, GL_RGBA, GL_FLOAT, NULL);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glBindTexture(GL_TEXTURE_2D, 0);

	// Create memory buffers on the device for each vector
	cl_mem cameraPosMemory = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(vec3), NULL, &ret);
	cl_mem cameraViewMatrixMemory = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(mat4), NULL, &ret);
	cl_mem lightDirectionMemory = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(vec3), NULL, &ret);
	
	cl_mem voxelSizeMemory = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(int) * 3, NULL, &ret); 

	printCLErrorCode("clCreateBuffer", ret);
	if (ret != 0) return 1;

	cl_mem outputMemory = clCreateFromGLTexture2D(context, CL_MEM_READ_WRITE, GL_TEXTURE_2D, 0, outputTexture, &ret);
	cl_mem screenMemory[2] = {
		clCreateFromGLTexture2D(context, CL_MEM_READ_WRITE, GL_TEXTURE_2D, 0, screenTexture[0], &ret),
		clCreateFromGLTexture2D(context, CL_MEM_READ_WRITE, GL_TEXTURE_2D, 0, screenTexture[1], &ret)
	};
	cl_mem finalMemory = clCreateFromGLTexture2D(context, CL_MEM_WRITE_ONLY, GL_TEXTURE_2D, 0, finalTexture, &ret);

	// Create the OpenCL kernel
	Kernel* kernelFast = new Kernel("Fast", "raytracer_fast.cl", "raytrace", context, device);
	Kernel* kernelNormal = new Kernel("AO", "raytracer_normal.cl", "raytrace", context, device);
	Kernel* kernelSlow = new Kernel("GI", "raytracer.cl", "raytrace", context, device);
	Kernel* kernelPath = new Kernel("Path", "raytracer_path.cl", "raytrace", context, device);

	vector<Kernel*>	raytracerKernels = {
		kernelFast,
		kernelNormal,
		kernelSlow,
		kernelPath
	};

	for (Kernel* kernel : raytracerKernels) kernel->build();

	Kernel postprocess("postprocess", "postprocess.cl", "main", context, device);
	postprocess.build();

	Kernel tonemapper("tonemapper", "tonemap.cl", "main", context, device);	
	tonemapper.build();

	bool	reloadKernel = true;
	int		kernelType = 0;
	Kernel* kernel = kernelFast;

	// Create voxels
	const int AIR = 0, GRASS = 1, GROUND = 2, STONE = 3, WATER = 4, WOOD = 5, LEAVES = 6;

	#pragma omp parallel for schedule(dynamic)
	for (int z = 0; z < voxelSize.z; z++)
		for (int y = 0; y < voxelSize.y; y++)
			for (int x = 0; x < voxelSize.x; x++)
			{
				float h = (perlinNoise(0.6f, 6, x / 128.0f, y / 128.0f, z / 128.0f, 0) - (0.5f - (float)y / voxelSize.y)) * 3.0f;

				int v = AIR;

				if (h >= 1.2f)
					v = STONE;
				else if (h >= 1.0f)
					v = GROUND;

				if (y > 0 && v == GROUND && voxels.get(x, y - 1, z) == AIR)
					v = GRASS;

				if (v == 0 && y > voxelSize.y - seaLevel)
					v = WATER;

				voxels.get(x, y, z) = v;
			}

	int trees = voxelSize.x*voxelSize.z*0.01f;

	int plantedTrees = 0;
	for (int i = 0; i < trees*trees && plantedTrees < trees; i++)
	{
		int x = noise(i, 10) * (voxelSize.x - 8) + 4;
		int y = 0;
		int z = noise(i, 20) * (voxelSize.z - 8) + 4;

		bool valid = true;
		for (y = 0; y < voxelSize.y; y++)
		{
			int type = voxels.get(x, y, z);
			valid = (type == GRASS);
			if (type != AIR) break;
		}

		if (!valid || y < 10) continue;

		plantedTrees++;

		for (int yt = y - 6; yt < y; yt++)
			voxels.get(x, yt, z) = WOOD;

		for (int yt = y - 8; yt < y - 2; yt++)
		{
			for (int xt = -4; xt <= 4; xt++)
				for (int zt = -4; zt <= 4; zt++)
				{
					int yft = y - yt - 6;
					int dist = zt*zt + xt*xt + yft*yft + perlinNoise(0.5, 3, (x+xt)/2.0f, yt, (z+zt)/2.0f, 13) * 8.0f;

					if (dist <= 12)
						voxels.get(x+xt, yt, z+zt) = LEAVES;
				}
		}
	}

	voxels.generateOctree();
	voxelsDirty = false;
	//voxels.remake();

	cout << "Node: " << sizeof(Node) << endl;

	// Create OpenCL data buffers
	int voxelsSpace = voxels.space * 1.2;
	int voxelsPtrSpace = voxels.ptrSpace * 1.5;
	size_t gpuMemoryEstimate = sizeof(uint8_t) * voxelCount + sizeof(Node) * voxelsSpace + sizeof(uint32_t) * voxelsPtrSpace;
	cout << "GPU Memory Estimate: " << ((gpuMemoryEstimate / 1024.0) / 1024) << "MB" << endl;

	cl_mem voxelMemory, octreeMemory, ptrTableMemory;
	voxelMemory = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint8_t) * voxelCount, NULL, &ret);
	octreeMemory = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(Node) * voxelsSpace, NULL, &ret);
	ptrTableMemory = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint32_t) * voxelsPtrSpace, NULL, &ret);

	// Main loop
	SDL_Event event;
	float delta = -1.0f;
	long last = SDL_GetTicks();
	long lastTimePhysics = last;
	long nextReport = last;
	bool running = true;

	bool postprocessEnabled = true;
	int	framesLeftToRender = -1;

	int forceMixOff = 2;
	int succesiveFrames = 1;

	while (running)
	{
		bool onceBomb = false;

		while (SDL_PollEvent(&event))
			switch (event.type)
			{
				case SDL_QUIT:
					running = false;
					break;
				case SDL_KEYUP:
					if (!grab) gui->keyPress(event.key.keysym.sym, false);

					key[event.key.keysym.sym] = false;

					if (event.key.keysym.sym == SDLK_ESCAPE)
					{
						running = false;
					}
					else if (event.key.keysym.sym == SDLK_p)
					{
						postprocessEnabled = !postprocessEnabled;
						forceMixOff = 1;
					}
					else if (event.key.keysym.sym == SDLK_t)
					{
						kernelType = (kernelType + 1) % raytracerKernels.size();
						kernel = raytracerKernels[kernelType];

						cout << "Kernel Type: " << kernelType << endl;
						cout << "Kernel: " << kernel->getName() << endl;

						reloadKernel = true;
					}
					else if (event.key.keysym.sym == SDLK_y)
					{
						if (framesLeftToRender == -1)
						{
							framesLeftToRender = 3;

							kernel = kernelPath;
							kernelType = 3;
						}
						else
						{
							framesLeftToRender = -1;

							kernel = kernelFast;
							kernelType = 0;
						}

						reloadKernel = true;
						delta = 0; // Reset delta smoothing
					}
					else if (event.key.keysym.sym == SDLK_l)
					{
						opencl = !opencl;
						delta = 0; // Reset delta smoothing
					}
					else if (event.key.keysym.sym == SDLK_g)
					{
						grab = !grab;

						if (grab)
						{
							SDL_ShowCursor(0);
							SDL_WM_GrabInput(SDL_GRAB_ON);
						}
						else
						{
							SDL_ShowCursor(1);
							SDL_WM_GrabInput(SDL_GRAB_OFF);
						}
					}
					else if (event.key.keysym.sym == SDLK_i)
					{
						cout << camera.position.x << ", " << camera.position.y << ", " << camera.position.z << endl;
					}
					else if (event.key.keysym.sym == SDLK_o)
					{
						ofstream out("camera.txt");
						out << camera.position.x << " " << camera.position.y << " " << camera.position.z;
						out << " " << camera.getYaw() << " " << camera.getPitch() << endl;
						out.close();
						cout << "Camera state saved" << endl;
					}
					else if (event.key.keysym.sym == SDLK_p)
					{
						ifstream in("camera.txt");
						float yaw, pitch;
						in >> camera.position.x >> camera.position.y >> camera.position.z; 
						in >> yaw >> pitch;
						in.close();
						camera.setYaw(yaw);
						camera.setPitch(pitch);
						cout << "Camera state restored" << endl;
					}

					//if (event.key.keysym.sym == SDLK_ESCAPE) running = false;
					break;
				case SDL_KEYDOWN:
					if (!grab) gui->keyPress(event.key.keysym.sym, true);

					if (event.key.keysym.sym == SDLK_n)
						onceBomb = true;

					key[event.key.keysym.sym] = true;
					break;
				case SDL_MOUSEBUTTONDOWN:
					if (!grab) 
					{
						gui->mouseButton(ivec2(event.button.x, event.button.y), true);
						gui->mouseMove(ivec2(event.button.x, event.button.y));
					}
					else if (event.button.button == 1) 
						deleteBlock();
					break;
				case SDL_MOUSEBUTTONUP:
					if (!grab) 
					{
						gui->mouseButton(ivec2(event.button.x, event.button.y), false);
						gui->mouseMove(ivec2(event.button.x, event.button.y));
					}
					break;
				case SDL_MOUSEMOTION:
					if (grab)
					{
						if (event.motion.xrel != 0 || event.motion.yrel != 0)
							forceMixOff = 1;

						camera.changeYaw(0.012f * event.motion.xrel);
						camera.changePitch(-0.015f * event.motion.yrel);
					}
					else
						gui->mouseMove(ivec2(event.motion.x, event.motion.y));
					break;
			}

		if (!running) break;

		if (framesLeftToRender > 0)
		{
			if (--framesLeftToRender == 0)
				cout << "Paused" << endl;
		}
		else if (framesLeftToRender == 0)
			continue;

		if (gui->getOption(1, 0).choice != kernelType)
		{
			kernelType = gui->getOption(1, 0).choice;
			kernel = raytracerKernels[kernelType];

			cout << "Kernel Type: " << kernelType << endl;
			cout << "Kernel: " << kernel->getName() << endl;

			reloadKernel = true;
		}

		{
			int newState = gui->getOption(2, 1).choice == 1;
			if (postprocessEnabled != newState)
			{
				postprocessEnabled = newState;
				forceMixOff = 1;
			}
		}

		{
			timeOfDay = gui->getOption(0, 0).value;

			//timeOfDay = (timeOfDay + delta * 0.025f);
			if (timeOfDay > 24.0f) timeOfDay -= 24.0f;
			float sunAngle = timeOfDay * M_PI / 12.0f;
			scene.lightDirection = normalize(vec3(sin(sunAngle), -cos(sunAngle), 0.25f));
		}

		{
			exposure = gui->getOption(2, 0).value;
			exposure *= exposure;
		}

		{
			float speed = 8.0f;

			//acceleration.y = 9.8f * 4.0f;
			velocity.x *= delta * 0.1f;
			velocity.y *= delta * 0.1f;
			velocity.z *= delta * 0.1f;

			vec3 up = vec3(0.0f, 1.0f, 0.0f);
			vec3 forward = camera.getWalkDirection();
			vec3 right = cross(camera.getWalkDirection(), up);

			if (key[SDLK_a])
			{
				velocity += right * speed;
				forceMixOff = 1;
			}
			if (key[SDLK_d])
			{
				velocity -= right * speed;
				forceMixOff = 1;
			}
			if (key[SDLK_w])
			{
				velocity += forward * speed;
				forceMixOff = 1;
			}
			if (key[SDLK_s])
			{
				velocity -= forward * speed;
				forceMixOff = 1;
			}
			if (key[SDLK_q])
			{
				velocity -= up * speed;
				forceMixOff = 1;
			}
			if (key[SDLK_z])
			{
				velocity += up * speed;
				forceMixOff = 1;
			}

			if (key[SDLK_b] || onceBomb)
			{
				Ray actionRay = Ray(camera.position, vec3(camera.getViewMatrix() * vec4(0.0f, 0.0f, 1.0f, 0.0f)));
				VoxelGrid::Result result = voxels.intersect(actionRay, 100.0f);

				if (result.hit)
				{
					voxelsDirty = true;

					for (int x = -4; x <= 4; x++)
						for (int y = -4; y <= 4; y++)
							for (int z = -4; z <= 4; z++)
								if (x*x + y*y + z*z < 16)
								{
									int xp = x + result.pos.x, yp = y + result.pos.y, zp = z + result.pos.z;
									if (xp < 0 || xp >= voxelSize.x || yp < 0 || yp >= voxelSize.y || zp < 0 || zp >= voxelSize.z) continue;
									voxels.get(xp, yp, zp) = 0;
								}
				}
			}
			else if (key[SDLK_DELETE])
			{
				deleteBlock();
			}
			else if (key[SDLK_RETURN])
			{
				Ray actionRay(camera.position, vec3(camera.getViewMatrix() * vec4(0.0f, 0.0f, 1.0f, 0.0f)));
				VoxelGrid::Result result = voxels.intersect(actionRay, 100.0f);

				if (result.hit)
					voxels.get(result.pos.x, result.pos.y, result.pos.z) = 3;
			}

			long now = SDL_GetTicks();
			bool firstPhysics = true;
			for (;lastTimePhysics < now; lastTimePhysics += 10)
			{
				float pdelta = 10.0f / 1000.0f;

				// Integrate position
				camera.position += velocity * pdelta + acceleration * 0.5f * pdelta * pdelta;
				velocity += acceleration * pdelta;

				// Move out of terrain
				Ray groundRay(camera.position, vec3(0.0f, 1.0f, 0.0f));
				VoxelGrid::Result groundResult = voxels.intersect(groundRay, 100.0f);

				if (groundResult.hit && groundResult.t < 2.0f)
				{
					camera.position.y -= 2.0f - groundResult.t;
					if (velocity.y > 0.0f) velocity.y = 0.0f;
					if (key[SDLK_SPACE] && velocity.y > -0.1f && firstPhysics) velocity.y = -20.0f;
				}

				firstPhysics = false;
			}
		}

		if (opencl)
		{
			static int bufferIndex = 1;
			bufferIndex = 1 - bufferIndex;

			if (reloadKernel)
			{
				reloadKernel = false;
				ret |= clSetKernelArg(kernel->kernel, 1, sizeof(cl_mem), (void *)&outputMemory);
				ret |= clSetKernelArg(kernel->kernel, 2, sizeof(cl_mem), (void *)&cameraPosMemory);
				ret |= clSetKernelArg(kernel->kernel, 3, sizeof(cl_mem), (void *)&cameraViewMatrixMemory);
				ret |= clSetKernelArg(kernel->kernel, 4, sizeof(cl_mem), (void *)&lightDirectionMemory);
				ret |= clSetKernelArg(kernel->kernel, 6, sizeof(cl_mem), (void *)&voxelSizeMemory);

				ret |= clEnqueueWriteBuffer(command_queue, voxelMemory, CL_TRUE, 0, sizeof(uint8_t) * voxelCount, voxels.getVoxels(), 0, NULL, NULL);
				ret |= clEnqueueWriteBuffer(command_queue, octreeMemory, CL_TRUE, 0, sizeof(Node) * voxels.space, voxels.nodes, 0, NULL, NULL);
				ret |= clEnqueueWriteBuffer(command_queue, ptrTableMemory, CL_TRUE, 0, sizeof(uint32_t) * voxels.ptrSpace, voxels.ptrTable, 0, NULL, NULL);

				ret |= clEnqueueWriteBuffer(command_queue, voxelSizeMemory, CL_TRUE, 0, sizeof(int) * 3, &voxelSize, 0, NULL, NULL);

				forceMixOff = 2;
			}

			if (forceMixOff > 0)
			{
				forceMixOff--;
				succesiveFrames = 1;
			}
			else if (succesiveFrames < 64)
				succesiveFrames++;

			float mixFactor = 1.0f - 1.0f / succesiveFrames;

			static int k = 0;
			k++;

			vec3 cameraPos = camera.position;
			mat4 cameraViewMatrix = transpose(camera.getViewMatrix());

			ret = clEnqueueWriteBuffer(command_queue, cameraPosMemory, CL_FALSE, 0, sizeof(vec3), &cameraPos,
					0, NULL, NULL);
			ret = clEnqueueWriteBuffer(command_queue, cameraViewMatrixMemory, CL_FALSE, 0, sizeof(mat4), &cameraViewMatrix,
					0, NULL, NULL);
			ret = clEnqueueWriteBuffer(command_queue, lightDirectionMemory, CL_FALSE, 0, sizeof(vec3), &scene.lightDirection,
					0, NULL, NULL);

			if (voxelsDirty)
			{
				voxelsDirty = false;
				voxels.generateOctree();

				ret |= clEnqueueWriteBuffer(command_queue, voxelMemory, CL_FALSE, 0, sizeof(uint8_t) * voxelCount, voxels.getVoxels(), 0, NULL, NULL);
				ret |= clEnqueueWriteBuffer(command_queue, octreeMemory, CL_FALSE, 0, sizeof(Node) * voxels.space, voxels.nodes, 0, NULL, NULL);
				ret |= clEnqueueWriteBuffer(command_queue, ptrTableMemory, CL_FALSE, 0, sizeof(uint32_t) * voxels.ptrSpace, voxels.ptrTable, 0, NULL, NULL);
			}

			// Execute the OpenCL kernel on the list
			ret = clSetKernelArg(kernel->kernel, 0, sizeof(int), (void *)&k);
			ret = clSetKernelArg(kernel->kernel, 5, sizeof(cl_mem), (void *)&voxelMemory);
			ret = clSetKernelArg(kernel->kernel, 7, sizeof(cl_mem), (void *)&octreeMemory);
			ret = clSetKernelArg(kernel->kernel, 8, sizeof(cl_mem), (void *)&ptrTableMemory);

			glFinish();
			ret = clEnqueueAcquireGLObjects(command_queue, 1, &outputMemory, 0, NULL, NULL);
			ret = clEnqueueAcquireGLObjects(command_queue, 1, &screenMemory[0], 0, NULL, NULL);
			ret = clEnqueueAcquireGLObjects(command_queue, 1, &screenMemory[1], 0, NULL, NULL);
			ret = clEnqueueAcquireGLObjects(command_queue, 1, &finalMemory, 0, NULL, NULL);

			size_t global_item_size[2];
			size_t local_item_size[2] = {8, 8};

			global_item_size[0] = WIDTH;//((WIDTH + 7) >> 3) << 3;
			global_item_size[1] = HEIGHT;//((HEIGHT + 15) >> 4) << 4;

			clEnqueueNDRangeKernel(command_queue, kernel->kernel, 2, NULL, global_item_size, local_item_size, 0, NULL, NULL);

			if (postprocessEnabled)
			{
				clSetKernelArg(postprocess.kernel, 0, sizeof(int), (void *)&k);
				clSetKernelArg(postprocess.kernel, 1, sizeof(float), (void *)&mixFactor);
				clSetKernelArg(postprocess.kernel, 2, sizeof(cl_mem), (void *)&outputMemory);
				clSetKernelArg(postprocess.kernel, 3, sizeof(cl_mem), (void *)&screenMemory[1-bufferIndex]);
				clSetKernelArg(postprocess.kernel, 4, sizeof(cl_mem), (void *)&screenMemory[bufferIndex]);

				clEnqueueNDRangeKernel(command_queue, postprocess.kernel, 2, NULL, global_item_size, local_item_size, 0, NULL, NULL);
			}
			
			if (postprocessEnabled)
				clSetKernelArg(tonemapper.kernel, 0, sizeof(cl_mem), (void *)&screenMemory[bufferIndex]);
			else
				clSetKernelArg(tonemapper.kernel, 0, sizeof(cl_mem), (void *)&outputMemory);

			clSetKernelArg(tonemapper.kernel, 1, sizeof(cl_mem), (void *)&finalMemory);
			clSetKernelArg(tonemapper.kernel, 2, sizeof(float), (void *)&exposure);

			clEnqueueNDRangeKernel(command_queue, tonemapper.kernel, 2, NULL, global_item_size, local_item_size, 0, NULL, NULL);

			clEnqueueReleaseGLObjects(command_queue, 1, &outputMemory, 0, NULL, NULL);
			clEnqueueReleaseGLObjects(command_queue, 1, &screenMemory[0], 0, NULL, NULL);
			clEnqueueReleaseGLObjects(command_queue, 1, &screenMemory[1], 0, NULL, NULL);
			clEnqueueReleaseGLObjects(command_queue, 1, &finalMemory, 0, NULL, NULL);
			clFinish(command_queue);

			drawToScreen();
		}
		else
		{
			static int k = 0;
			#pragma omp parallel for schedule(dynamic) lastprivate(k)
			for (int y = 0; y < HEIGHT; y++)
				for (int x = 0; x < WIDTH; x++)
					pixels[x + y * WIDTH] = raytracer->trace(x, y, x + y * WIDTH + k++);

			glBindTexture(GL_TEXTURE_2D, screenTexture[0]);
			glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, WIDTH, HEIGHT, GL_RGBA, GL_FLOAT, (GLvoid*)pixels);
			glBindTexture(GL_TEXTURE_2D, 0);

			drawToScreen();
		}

		long now = SDL_GetTicks();

		//if (now - start > 10000) break;

		float cdelta = (now - last) * 0.001f;
		last = now;

		delta = (delta > 0.0f) ? (delta * 2.0f + cdelta) * 0.333f : cdelta;

		if (nextReport < now)
		{
			cout << fixed << setprecision(2) << (delta * 1000) << "ms - FPS: " << (1.0f / delta) << endl;
			nextReport = now + 500;
		}

		gui->update(delta, !grab);
	}

	glDeleteTextures(2, screenTexture);

	// Clean up
	for (Kernel* kernel : raytracerKernels) delete kernel;

	ret = clFinish(command_queue);

	ret = clReleaseMemObject(cameraPosMemory);
	ret = clReleaseMemObject(cameraViewMatrixMemory);
	ret = clReleaseMemObject(lightDirectionMemory);
	ret = clReleaseMemObject(outputMemory);
	ret = clReleaseMemObject(screenMemory[0]);
	ret = clReleaseMemObject(screenMemory[1]);
	ret = clReleaseMemObject(voxelSizeMemory);

	ret = clReleaseMemObject(voxelMemory);
	ret = clReleaseMemObject(octreeMemory);
	ret = clReleaseMemObject(ptrTableMemory);

	ret = clReleaseCommandQueue(command_queue);
	ret = clReleaseContext(context);
	free(pixels);

	if (grab)
	{
		SDL_ShowCursor(1);
		SDL_WM_GrabInput(SDL_GRAB_OFF);
	}

	TTF_Quit();
	SDL_Quit();

	return 0;
}
