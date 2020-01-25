NB_PROC_MPI=4

CC = mpic++
CCF = mpif90

CFLAGS = -I./include
CFLAGS += -std=c++11 -Wall -Werror
#CFLAGS +=  -g -O0
CFLAGS += -O3

#LDFLAGS = -fopenmp

SRC_FILES = $(notdir $(wildcard src2/*.cpp))
OBJ_FILES = $(addprefix obj/,$(SRC_FILES:.cpp=.o))

SRCF_FILES = $(wildcard src/*.f90) #$(notdir $(wildcard src/*.f90))
#OBJF_FILES = $(addprefix obj/f90/,$(SRCF_FILES:.f90=.o))


all: bin/main.exe# bin/mainf90.exe


bin/main.exe: $(OBJ_FILES)
	$(CC) -o $@ $^ $(LDFLAGS)

bin/mainf90.exe: $(SRCF_FILES)
	$(CCF) $^ $(LDFLAGSF)
	mv a.out bin/mainf90
	rm matrice.mod
	rm function.mod

obj/%.o: src2/%.cpp
	$(CC) $< -o $@ -c $(CFLAGS)

#obj/f90/%.o: src/%.f90
#	$(CCF) $< -o $@ -c $(CFLAGSF)

.PHONY: install run clean clean_result

install:
	mkdir -p bin obj Result #obj/f90

run: all
	mpirun -n $(NB_PROC_MPI) --mca pml ob1 ./bin/main.exe data.txt

clean:
	rm -f bin/* obj/*.o obj/f90/*.o

clean_result:
	rm -f Result/*
