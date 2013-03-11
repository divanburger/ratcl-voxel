#include "Texture.h"

Texture::Texture(string filename)
{
	load(filename);
}

void Texture::load(string filename)
{
	SDL_Surface* data = IMG_Load(filename.c_str());

	if (!data)
	{
		cerr << "Could not load image: " << filename << endl;
		return;
	}

	width = data->w;
	height = data->h;

	glGenTextures(1, &internalTexture);
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, internalTexture);

	alpha = (data->format->BitsPerPixel == 32);
	GLenum format = alpha ? GL_RGBA : GL_RGB;

	glTexImage2D(GL_TEXTURE_2D, 0, GL_SRGB_ALPHA, width, height, 0, format, GL_UNSIGNED_BYTE, data->pixels);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

	glBindTexture(GL_TEXTURE_2D, 0);
	glDisable(GL_TEXTURE_2D);

	SDL_FreeSurface(data);
}

Texture::~Texture()
{
	glDeleteTextures(1, &internalTexture);
}
