CC=ccache g++

DCFLAGS=-pipe -Wall -std=c++0x -g
PCFLAGS=-pipe -Wall -std=c++0x -g -pg
OCFLAGS=-pipe -Wall -std=c++0x -fno-rtti -Ofast -ffast-math -march=native -fopenmp -flto
CFLAGS=-Ilibs $(OCFLAGS)
LFLAGS=
LIBS=-lGL -lGLEW -lSDL -lSDL_image -lSDL_ttf -lOpenCL

BIN_SRC=$(shell find src -iname main.cpp)
BIN=RatCL

SRC=$(shell find src -iname *.cpp)
SOURCES=$(SRC:src/%.cpp=%.cpp)
OBJECTS=$(SOURCES:%.cpp=objs/%.o)
DEPS=$(SOURCES:%.cpp=deps/%.d)

all: release

release: $(OBJECTS)
	@echo Making release: $(BIN)
	$(CC) $(CFLAGS) $(LFLAGS) -o $(BIN) $(OBJECTS) $(LIBS)

-include $(DEPS)

objs/%.o: src/%.cpp deps/%.d
	@echo Compiling: $*
	@mkdir -p $(@D)
	@$(CC) $(CFLAGS) -c $< -o $@

depend: $(DEPS)

deps/%.d: src/%.cpp
	@echo Generating dependencies: $*
	@mkdir -p $(@D)
	@$(CC) $(CFLAGS) -MM -MT 'objs/$*.o' $< > deps/$*.d

clean:
	rm -rf deps
	rm -rf objs
	rm -ff $(BIN)
