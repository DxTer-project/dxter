# indicate where the object files are to be created
CC         := mpicxx
LINKER     := $(CC)
#CFLAGS     := -O3 -fopenmp -Wall -Isrc/ -Isrc/DLA/ -Isrc/tensors/ -Isrc/LLDLA
#CFLAGS     := -O3 -Wall
CFLAGS	   := -g -Wall -I/Users/bamarker/work/tensormental/build/include
#CFLAGS	   := -g -Wall
#CFLAGS	   := -pg -Wall

TENSORMENTALDIR := /Users/bamarker/work/tensormental/build
TENSORMENTALLIB := $(TENSORMENTALDIR)/libtensormental.a

HEADERS :=  $(shell find src -type f -name '*.h')
SOURCES :=  $(shell find src -type f -name '*.cpp')
OBJS := $(patsubst %.cpp, %.o, $(SOURCES))

all: Test1.x

Test1.o: Test1.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

Test1.x: Test1.o $(HEADERS)
	$(LINKER) $(CFLAGS) Test1.o $(TENSORMENTALLIB) -o $@

clean:
	rm -f *.o *.x *~