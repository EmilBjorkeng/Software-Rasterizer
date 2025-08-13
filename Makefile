CC = g++
CFLAGS = -Isrc/include -std=c++26 -Wall
LDFLAGS = -Lsrc/lib
LDLIBS = -lmingw32 -lSDL3 -lSDL3_ttf

TARGET = main
SRC = $(filter-out src/$(TARGET).cpp, $(wildcard src/*.cpp))
OBJ = $(SRC:src/%.cpp=%.o)

.PHONY: all clean run debug

all: $(TARGET)

%.o: src/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) $^ src/$(TARGET).cpp $(LDFLAGS) $(LDLIBS) -o $@

clean:
	-del $(TARGET).exe 2>nul || true
	-del *.o 2>nul || true

run: clean all
	.\$(TARGET).exe
	$(MAKE) clean

debug: CFLAGS += -g
debug: clean all