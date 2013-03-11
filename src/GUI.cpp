#include "GUI.h"

const vec3 GUI::WHITE = vec3(1.0f, 1.0f, 1.0f);
const vec3 GUI::BLACK = vec3(0.0f, 0.0f, 0.0f);
const vec3 GUI::BLUE = vec3(0.1f, 0.3f, 0.6f);
const vec4 GUI::BLACK_TRANS = vec4(0.0f, 0.0f, 0.0f, 0.9f);

GUI::GUI(ivec2 screenSize) : screenSize(screenSize), frames(0), intervalTime(0), dragging(false), selectedMenu(-1)
{
	defaultFont = new FontFace("fonts/DejaVuSans.ttf", 16);

	Option option;

	Menu lightMenu;
	lightMenu.name = "Light";
		option.name = "Time of Day";
		option.minValue = 0.0f;
		option.maxValue = 24.0f;
		option.value = 8.0f;
		lightMenu.options.push_back(option);

		option.name = "Sun Brightness";
		option.minValue = 1.0f;
		option.maxValue = 20.0f;
		option.value = 2.5f;
		lightMenu.options.push_back(option);
	menus.push_back(lightMenu);

	Menu renderMenu;
	renderMenu.name = "Render";
		Option kernelOption;
		kernelOption.name = "Type";
		kernelOption.type = Option::CHOICE;
		kernelOption.choice = 0;
		kernelOption.choices.push_back("Fast");
		kernelOption.choices.push_back("AO");
		kernelOption.choices.push_back("GI");
		kernelOption.choices.push_back("Path");
		renderMenu.options.push_back(kernelOption);	
	menus.push_back(renderMenu);

	Menu postMenu;
	postMenu.name = "Post";
		option.name = "Exposure";
		option.minValue = 0.0f;
		option.maxValue = 1.5f;
		option.value = 0.35f;
		postMenu.options.push_back(option);

		Option persistenceOption;
		persistenceOption.name = "Persistence";
		persistenceOption.type = Option::CHOICE;
		persistenceOption.choice = 1;
		persistenceOption.choices.push_back("Off");
		persistenceOption.choices.push_back("On");
		postMenu.options.push_back(persistenceOption);
	menus.push_back(postMenu);	
}

GUI::~GUI()
{

}

void GUI::update(float delta, bool editMode)
{
	frames++;
	intervalTime += delta;

	if (intervalTime > 0.5f)
	{
		fps = frames / intervalTime;
		intervalTime = 0.0f;
		frames = 0;
	}
}

void GUI::keyPress(int key, bool down)
{
	if (key == SDLK_ESCAPE && down) selectedMenu = -1;
}

void GUI::mouseMove(ivec2 pos)
{
	if (dragging && pos.y > 40 && pos.x >= screenSize.x - 300 && selectedMenu != -1)
	{
		int leftSide = screenSize.x - 300;

		int row = (pos.y - 40) / 60;

		Menu& menu = menus[selectedMenu];
		if (row >= 0 && row < menu.options.size())
		{
			Option& option = menu.options[row];

			if (option.type == Option::VALUE)
			{
				float ratio = glm::clamp((pos.x - leftSide - 10) / 280.0f, 0.0f, 1.0f);
				option.value = option.minValue + ratio * (option.maxValue - option.minValue);
			}
			else
			{
				option.choice = glm::clamp(
					int((pos.x - leftSide - 10) * option.choices.size() / 280), 0, (int)option.choices.size());
			}
		}
	}
}

void GUI::mouseButton(ivec2 pos, bool down)
{
	dragging = down;

	if (pos.y < 40)
	{
		selectedMenu = -1;

		for (size_t i = 0; i < menus.size(); i++)
		{
			int x = screenSize.x + (i - menus.size()) * 80;

			if (pos.x >= x && pos.x < x + 80) selectedMenu = i;
		}
	}
}

GUI::Option& GUI::getOption(int menu, int option)
{
	return menus[menu].options[option];
}

void GUI::render(bool editMode)
{
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glLoadIdentity();
	glOrtho(0.0f, screenSize.x, screenSize.y, 0.0f, -1.0f, 1.0f);

	if (editMode)
	{
		rectFill(ivec2(0, 0), ivec2(screenSize.x, 40), BLUE);

		for (size_t i = 0; i < menus.size(); i++)
		{
			int x = screenSize.x + (i - menus.size()) * 80;

			if (i == selectedMenu)
			{
				rectFill(ivec2(x, 0), ivec2(80, 40), WHITE);
				renderText(ivec2(10+x, 10), menus[i].name, defaultFont, BLACK);
			}
			else
				renderText(ivec2(10+x, 10), menus[i].name, defaultFont, WHITE);
		}

		if (selectedMenu != -1)
		{
			int leftSide = screenSize.x - 300;
			rectFill(ivec2(leftSide, 40), ivec2(300, screenSize.y - 40), BLACK_TRANS);

			Menu& menu = menus[selectedMenu];
			for (size_t i = 0; i < menu.options.size(); i++)
			{
				Option& option = menu.options[i];
				float ratio = (option.value - option.minValue) / (option.maxValue - option.minValue);

				ivec2 start = ivec2(leftSide, 40 + 60 * i);
				renderText(start + ivec2(10, 10), option.name, defaultFont, WHITE);

				string valueText = "<unknown>";
				if (option.type == Option::VALUE)
					valueText = toString(option.value, 3);
				else
					valueText = option.choices[option.choice];

				int valueTextWidth = getTextWidth(valueText, defaultFont);
				renderText(ivec2(screenSize.x - 10 - valueTextWidth, start.y + 10), valueText, defaultFont, WHITE);

				if (option.type == Option::VALUE)
				{
					rectFill(start + ivec2(10, 41), ivec2(280, 3), BLUE);
					rectFill(start + ivec2(10, 38) + ivec2(ratio * (280 - 9), 0), ivec2(9, 9), WHITE);
				}
				else
				{
					int choiceWidth = 290 / option.choices.size() - 10;

					for (int i = 0; i < option.choices.size(); i++)
					{
						int x = (290 * i) / option.choices.size();
						rectFill(start + ivec2(10 + x, 35), ivec2(choiceWidth, 25), (i == option.choice) ? WHITE : BLUE);
						renderText(start + ivec2(15 + x, 40), option.choices[i], defaultFont, (i == option.choice) ? BLACK : WHITE);
					}
				}
			}
		}
	}

	renderText(ivec2(10, 10), "FPS: " + toString(fps, 4), defaultFont, WHITE);

	glDisable(GL_BLEND);
}

void GUI::rectFill(ivec2 pos, ivec2 size, vec3 colour)
{
	glColor3f(colour.r, colour.g, colour.b);
	glBegin(GL_QUADS);
		glVertex2f(pos.x, pos.y);
		glVertex2f(pos.x+size.x, pos.y);
		glVertex2f(pos.x+size.x, pos.y+size.y);
		glVertex2f(pos.x, pos.y+size.y);
	glEnd();
}

void GUI::rectFill(ivec2 pos, ivec2 size, vec4 colour)
{
	glColor4f(colour.r, colour.g, colour.b, colour.a);
	glBegin(GL_QUADS);
		glVertex2f(pos.x, pos.y);
		glVertex2f(pos.x+size.x, pos.y);
		glVertex2f(pos.x+size.x, pos.y+size.y);
		glVertex2f(pos.x, pos.y+size.y);
	glEnd();
}

void GUI::renderText(ivec2 pos, string text, FontFace* font, vec3 colour)
{
	SDL_Surface* data = TTF_RenderText_Blended(font->internalFont, text.c_str(), SDL_Colour { 255, 255, 255, 255 });

	GLuint texture;
	glGenTextures(1, &texture);
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, texture);

	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, data->w, data->h, 0, GL_RGBA, GL_UNSIGNED_BYTE, data->pixels);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

	glColor3f(colour.r, colour.g, colour.b);

	glBegin(GL_QUADS);
	glTexCoord2f(0.0f, 0.0f);
	glVertex2i(pos.x, pos.y);
	glTexCoord2f(1.0f, 0.0f);
	glVertex2i(pos.x + data->w, pos.y);
	glTexCoord2f(1.0f, 1.0f);
	glVertex2i(pos.x + data->w, pos.y + data->h);
	glTexCoord2f(0.0f, 1.0f);
	glVertex2i(pos.x, pos.y + data->h);
	glEnd();

	glBindTexture(GL_TEXTURE_2D, 0);
	glDisable(GL_TEXTURE_2D);

	glDeleteTextures(1, &texture);

	SDL_FreeSurface(data);
}

int GUI::getTextWidth(string text, FontFace* font)
{
	int width, height;
	TTF_SizeText(font->internalFont, text.c_str(), &width, &height);
	return width;
}