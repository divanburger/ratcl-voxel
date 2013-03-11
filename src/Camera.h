/*
 * Camera.h
 *
 *  Created on: 31 Oct 2011
 *      Author: dburger
 */

#ifndef CAMERA_H_
#define CAMERA_H_

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

using namespace glm;

class Camera
{
	public:
		Camera(vec3 position, float yaw = 0.0f, float pitch = 0.0f);

		void moveForward(float amount);
		void changeYaw(float angle);
		void changePitch(float angle);

		mat4 getViewMatrix();
		vec3 getWalkDirection();

		float getYaw() const {return yaw;}
		float getPitch() const {return pitch;}

		void setYaw(float yaw) {this->yaw = yaw; dirtyMatrix = true;}
		void setPitch(float pitch) {this->pitch = pitch; dirtyMatrix = true;}

		vec3	position;

	private:
		bool	dirtyMatrix;

		mat4	matrix;
		float	yaw;
		float	pitch;
};

#endif /* CAMERA_H_ */
