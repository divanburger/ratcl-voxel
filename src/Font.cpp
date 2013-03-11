#include "Font.h"

FontFace::FontFace(string filename, int ptsSize)
{
	internalFont = TTF_OpenFont(filename.c_str(), ptsSize);
	if (!internalFont) 
		cerr << "Could not load font: " << filename << endl;
}

FontFace::~FontFace()
{
	//TTF_CloseFont(internalFont);
}
