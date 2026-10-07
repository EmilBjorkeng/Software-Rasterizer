CC = g++
CXXFLAGS = -Isrc/include -std=c++26 -Wall -Wextra -O2
PKG_CFLAGS := $(shell pkg-config --cflags sdl3 sdl3-ttf)
PKG_LDFLAGS := $(shell pkg-config --libs sdl3 sdl3-ttf)

TARGET = main
SRC = $(wildcard src/*.cpp)
OBJ = $(SRC:src/%.cpp=%.o)

.PHONY: all clean run debug

all: $(TARGET)$(EXE)

%.o: src/%.cpp
	$(CC) $(CXXFLAGS) $(PKG_CFLAGS) -c $< -o $@

$(TARGET)$(EXE): $(OBJ)
	$(CC) $(CXXFLAGS) $^ -o $@ $(PKG_LDFLAGS)

clean:
	-rm -f $(TARGET)
	-rm -f *.o

run:
	$(MAKE) -j$(NPROC) all
	./$(TARGET)

debug: clean
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -g"
