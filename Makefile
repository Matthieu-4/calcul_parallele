CC = mpic++

CFLAGS = -I./include
CFLAGS += -std=c++11
CFLAGS +=  -g -O0

#LDFLAGS = -fopenmp

SRC_FILES = $(notdir $(wildcard src/*.cpp))
OBJ_FILES = $(addprefix obj/,$(SRC_FILES:.cpp=.o))

all: bin/main.exe


bin/main.exe: $(OBJ_FILES)
	$(CC) -o $@ $^ $(LDFLAGS)

obj/%.o: src/%.cpp
	$(CC) $< -o $@ -c $(CFLAGS)



.PHONY: install clean

install:
	mkdir -p bin obj


clean:
	rm -f bin/* obj/*
