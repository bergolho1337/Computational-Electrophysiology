INCLUDE_DIR = include
SOURCE_DIR = src
OBJECT_DIR = obj
CC = g++
CFLAGS = -Wall -O2 -std=c++11 -I$(INCLUDE_DIR)
PFLAGS = -DOUTPUT
CC_SOURCES = $(wildcard src/*)
STRING_OBJ_AUX = $(CC_SOURCES:.cpp=.o)
STRING_OBJ = $(subst src/,,$(STRING_OBJ_AUX))
CC_OBJ = $(patsubst %,$(OBJECT_DIR)/%,$(STRING_OBJ))
PROGRAM_NAME = celular_solver

all: $(PROGRAM_NAME)

$(PROGRAM_NAME): $(CC_OBJ)
	$(CC) $(CFLAGS) $(PFLAGS) -o $@ $^

$(OBJECT_DIR)/%.o: $(SOURCE_DIR)/%.cpp
	$(CC) $(CFLAGS) $(PFLAGS) -o $@ -c $<

clean:
	rm -f $(OBJECT_DIR)/*.o Output/*.dat Output/*.pdf Output/*.log $(PROGRAM_NAME)

plot:
	python Output/plotter.py

print-%  : ; @echo $* = $($*)
