#ifndef _GUI_h
#define _GUI_h

#include <SDL/SDL.h>
#include <SDL/SDL_image.h>
#include <SDL/SDL_ttf.h>
#include <GL/glew.h>
#include <glm/glm.hpp>

#include <iostream>
#include <string>
#include <vector>

#include "Font.h"
#include "Util.h"

using namespace std;
using namespace glm;

class GUI
{
	public:
		struct Option
		{
			enum Type {VALUE, CHOICE};

			Option() : type(VALUE) {}
			
			string	name;
			Type	type;

			vector<string> choices;

			union
			{
				struct {float maxValue, minValue, value;};
				struct {int choice;};
			};
		};

		struct Menu
		{
			string			name;
			vector<Option>	options;
		};

		GUI(ivec2 screenSize);
		~GUI();

		void update(float delta, bool editMode);
		void render(bool editMode);

		void keyPress(int key, bool down);
		void mouseMove(ivec2 pos);
		void mouseButton(ivec2 pos, bool down);

		Option&	getOption(int menu, int option);

	private:
		void 	rectFill(ivec2 pos, ivec2 size, vec3 colour);
		void 	rectFill(ivec2 pos, ivec2 size, vec4 colour);

		void 	renderText(ivec2 pos, string text, FontFace* font, vec3 colour);
		int		getTextWidth(string text, FontFace* font);

		ivec2 		screenSize;

		FontFace*	defaultFont;

		float		frames, intervalTime;

		float		fps, frameTime;

		bool			dragging;
		int				selectedMenu;
		vector<Menu>	menus;

		static const vec3	WHITE;
		static const vec3	BLACK;
		static const vec3	BLUE;
		static const vec4	BLACK_TRANS;
};

#endif