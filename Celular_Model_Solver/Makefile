INCLUDE_DIR = include
SOURCE_DIR = src
OBJECT_DIR = obj
CC = g++
CFLAGS = -Wall -O2 -I$(INCLUDE_DIR)
CC_SOURCES = $(wildcard src/*)
STRING_OBJ_AUX = $(CC_SOURCES:.cpp=.o)
STRING_OBJ = $(subst src/,,$(STRING_OBJ_AUX))
CC_OBJ = $(patsubst %,$(OBJECT_DIR)/%,$(STRING_OBJ))
PROGRAM_NAME = celular_solver

all: $(PROGRAM_NAME)

$(PROGRAM_NAME): $(CC_OBJ)
	$(CC) $(CFLAGS) -o $@ $^

$(OBJECT_DIR)/%.o: $(SOURCE_DIR)/%.cpp
	$(CC) $(CFLAGS) -o $@ -c $<

clean:
	rm $(OBJECT_DIR)/*.o *.dat $(PROGRAM_NAME)

print-%  : ; @echo $* = $($*)
