# indicate where the object files are to be created
CC         := g++
LINKER     := $(CC)
#CFLAGS     := -O3 -fopenmp -Wall -Isrc/ -Isrc/DLA/ -Isrc/tensors/ -Isrc/LLDLA
#CFLAGS     := -O3 -Wall
CFLAGS	   := -g -Wall -Isrc/ -Isrc/DLA/ -Isrc/tensors -Isrc/LLDLA
#CFLAGS	   := -g -Wall
#CFLAGS	   := -pg -Wall

HEADERS :=  $(shell find src -type f -name '*.h')
SOURCES :=  $(shell find src -type f -name '*.cpp')
OBJS := $(patsubst src/%.cpp, obj/%.o, $(SOURCES))

all: dxter.x

obj/%.o: src/%.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

dxter.x: $(OBJS) $(HEADERS)
	$(LINKER) $(CFLAGS) $(OBJS) -o $@

clean:
	rm -f obj/*.o obj/DLA/*.o obj/tensors/*.o src/*~ src/DLA/*~ src/LLDLA/*~ src/LLDLA/*.o *.x *~

open:
	emacs -nw src/*cpp src/*h src/DLA/*cpp src/DLA/*h src/tensors/*cpp src/tensors/*h src/LLDLA/*cpp src/LLDLA/*h makefile

opencpp:
	emacs -nw src/*cpp src/DLA/*cpp src/tensors/*cpp src/LLDLA/*cpp

openh:
	emacs -nw src/*h src/DLA/*h src/tensors/*h src/LLDLA/*h