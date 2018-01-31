CC = g++
CFLAGS = -Wall -Wno-write-strings -std=c++11 -O2
PROGRAM_NAME = noble

all: $(PROGRAM_NAME)

$(PROGRAM_NAME): $(PROGRAM_NAME).cpp
	$(CC) $(CFLAGS) -o $@ $^ -lm

clean:
	rm -f $(PROGRAM_NAME)

clcResults:
	rm -f solution.dat

remade:
	$(MAKE) clean
	$(MAKE)

plot:
	python plotter.py

print-%  : ; @echo $* = $($*)
