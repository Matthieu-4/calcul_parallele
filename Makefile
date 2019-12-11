CC = mpic++
# CCF = mpif90

CFLAGS = -I./include
CFLAGS += -std=c++11
CFLAGS +=  -g -O0


#LDFLAGS = -fopenmp

SRC_FILES = $(notdir $(wildcard src/*.cpp))
OBJ_FILES = $(addprefix obj/,$(SRC_FILES:.cpp=.o))

SRCF_FILES = $(wildcard src/*.f90) #$(notdir $(wildcard src/*.f90))
#OBJF_FILES = $(addprefix obj/f90/,$(SRCF_FILES:.f90=.o))


all: bin/main.exe


bin/main.exe: $(OBJ_FILES)
	$(CC) -o $@ $^ $(LDFLAGS)

# bin/mainf90.exe: $(SRCF_FILES)
# 	$(CCF) $^ $(LDFLAGSF)
# 	mv a.out bin/mainf90
# 	rm matrice.mod
# 	rm function.mod

obj/%.o: src/%.cpp
	$(CC) $< -o $@ -c $(CFLAGS)

#obj/f90/%.o: src/%.f90
#	$(CCF) $< -o $@ -c $(CFLAGSF)

.PHONY: install run clean

install:
	mkdir -p bin obj Result #obj/f90

run: all
	mpirun -n 2 --mca pml ob1 ./bin/main.exe data.txt

clean:
	rm -f bin/* obj/*.o obj/f90/*.o

plot:
	gnuplot Result/sol.plot
