CC = mpic++

CFLAGS = -I./include
CFLAGS += -std=c++11
#CFLAGS +=  -g -O0

CFLAGS += -O3


SRC_FILES = $(notdir $(wildcard src/*.cpp))
OBJ_FILES = $(addprefix obj/,$(SRC_FILES:.cpp=.o))

all: bin/main.exe


bin/main.exe: $(OBJ_FILES)
	$(CC) -o $@ $^ $(LDFLAGS)

obj/%.o: src/%.cpp
	$(CC) $< -o $@ -c $(CFLAGS)


.PHONY: install run clean clean_result

install:
	mkdir -p bin obj Result pdf

# Use default settings for overlap and method
H =
V =
run: all
	mpirun -n 8 --mca pml ob1 ./bin/main.exe data.txt $(V) $(H)

clean:
	rm -f bin/* obj/*.o

plot:
	gnuplot Result/sol.plot

NB_PORC = 2 4 8
VERIONS = 0 1 2
OUT = Result/data_speedup.csv

data_speedup: bin/main.exe
	@echo "v,n,t,e" > $(OUT)
	@for v in $(VERIONS); do \
			for p in $(NB_PORC); do \
				echo -n "$$v,$$p," >> $(OUT) ; \
				mpirun -n $$p --mca pml ob1 ./bin/main.exe data.txt $$v >> $(OUT); \
			done ; \
	done

OUT2 = Result/data_h.csv
SIZE_H = 1 2 3 4 5 6 7 8 9 10
data_h: bin/main.exe
	@echo "v,h,t,e" > $(OUT2)
	@for v in $(VERIONS); do \
			for h in $(SIZE_H); do \
				echo -n "$$v,$$h," >> $(OUT2) ; \
				mpirun -n 8 --mca pml ob1 ./bin/main.exe data.txt $$v >> $(OUT2); \
			done ; \
	done

h: data_h
	R CMD BATCH graph/h.R

speedup: data_speedup
	R CMD BATCH graph/speedup.R

error: data_speedup
	R CMD BATCH graph/error.R

graph: h speedup error
