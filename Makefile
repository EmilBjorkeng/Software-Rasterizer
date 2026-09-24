CC = g++
CXXFLAGS = -Isrc/include -std=c++26 -Wall -Wextra
PKG_CFLAGS := $(shell pkg-config --cflags sdl3 sdl3-ttf)
PKG_LDFLAGS := $(shell pkg-config --libs sdl3 sdl3-ttf)

TARGET = main
SRC = $(wildcard src/*.cpp)
OBJ = $(SRC:src/%.cpp=%.o)

ifeq ($(OS),Windows_NT)
	RM = del /Q
	EXE = .exe
	RUN_CMD = .\$(TARGET)$(EXE)
else
	RM = rm -f
	EXE =
	RUN_CMD = ./$(TARGET)
endif

.PHONY: all clean run debug

all: $(TARGET)$(EXE)

%.o: src/%.cpp
	$(CC) $(CXXFLAGS) $(PKG_CFLAGS) -c $< -o $@

$(TARGET)$(EXE): $(OBJ)
	$(CC) $(CXXFLAGS) $^ -o $@ $(PKG_LDFLAGS)

clean:
	-$(RM) $(TARGET)$(EXE)
	-$(RM) *.o

run: all
	$(RUN_CMD)

debug: clean
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -g"
