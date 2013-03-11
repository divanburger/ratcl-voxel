/*
 * Texture.h
 *
 *  Created on: 8 Nov 2011
 *      Author: dburger
 */

#ifndef TEXTURE_H_
#define TEXTURE_H_

#include <SDL/SDL.h>
#include <SDL/SDL_image.h>
#include <GL/glew.h>
#include <glm/glm.hpp>

#include <iostream>
#include <string>

using namespace std;
using namespace glm;

class Texture
{
	public:
		Texture(string filename);
		virtual ~Texture();

		void load(string filename);

		int getWidth() const {return width;}
		int	getHeight() const {return height;}
		ivec2 getSize() const {return ivec2(width, height);}

		GLuint	internalTexture;

	private:
		int		width;
		int		height;
		bool	alpha;
};

#endif /* TEXTURE_H_ */
