#ifndef _Font_h
#define _Font_h

#include <SDL/SDL.h>
#include <SDL/SDL_ttf.h>

#include <iostream>
#include <string>

using namespace std;

class FontFace
{
	public:
		FontFace(string filename, int ptsSize);
		virtual ~FontFace();

		TTF_Font*	internalFont;
};

#endif
