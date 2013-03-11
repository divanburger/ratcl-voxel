/*
 * Camera.cpp
 *
 *  Created on: 31 Oct 2011
 *      Author: dburger
 */

#include "Camera.h"

Camera::Camera(vec3 position, float yaw, float pitch) : position(position), dirtyMatrix(true), yaw(yaw), pitch(pitch)
{
}

mat4 Camera::getViewMatrix()
{
	if (!dirtyMatrix) return matrix;
	dirtyMatrix = false;

	matrix = mat4(1.0f);
	matrix = rotate(matrix, degrees(yaw), vec3(0.0f, 1.0f, 0.0f));
	matrix = rotate(matrix, degrees(pitch), vec3(1.0f, 0.0f, 0.0f));
	return matrix;
}

vec3 Camera::getWalkDirection()
{
	return vec3(sin(yaw), 0.0f, cos(yaw));
}

void Camera::moveForward(float amount)
{
	position += vec3(getViewMatrix() * vec4(amount, 0.0f, 0.0f, 0.0f));
}

void Camera::changeYaw(float angle)
{
	dirtyMatrix = true;
	yaw += angle;
}

void Camera::changePitch(float angle)
{
	dirtyMatrix = true;
	pitch += angle;
	if (pitch < -M_PI*0.499f) pitch = -M_PI*0.499f;
	if (pitch > M_PI*0.499f) pitch = M_PI*0.499f;
}
