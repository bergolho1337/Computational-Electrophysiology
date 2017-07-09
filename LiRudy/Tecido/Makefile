CPP = g++
CPP_FLAGS = -Wall -O2
LINKER = -lm

all: prd

prd: main.c prd.c
	$(CPP) $(CPP_FLAGS) -o $@ $^ $(LINKER)

clear:
	rm -r prd

plot:
	python plot.py
