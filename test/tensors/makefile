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

BLIS_INSTALL_PREFIX := /Users/bamarker/work/blisInstall
BLIS_LIB_PATH  := $(BLIS_INSTALL_PREFIX)/lib
BLIS_LIB       := $(BLIS_LIB_PATH)/libblis.a

#CFLAGS += -fno-underscoring -fno-second-underscore 


all: Test1.x Test2.x Test3.x MP3.x

Test1.o: Test1.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

Test1.x: Test1.o $(HEADERS) $(TENSORMENTALLIB)
	$(LINKER) $(CFLAGS) Test1.o $(TENSORMENTALLIB) -framework Accelerate -o $@


Test2.o: Test2.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

Test2.x: Test2.o $(HEADERS) $(TENSORMENTALLIB)
	$(LINKER) $(CFLAGS) Test2.o $(TENSORMENTALLIB) -framework Accelerate -o $@



Test3.o: Test3.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

Test3.x: Test3.o $(HEADERS) $(TENSORMENTALLIB)
	$(LINKER) $(CFLAGS) Test3.o $(TENSORMENTALLIB) -framework Accelerate -o $@

MP3.o: MP3.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

MP3.x: MP3.o $(HEADERS) $(TENSORMENTALLIB)
	$(LINKER) $(CFLAGS) MP3.o $(TENSORMENTALLIB) -framework Accelerate -o $@

clean:
	rm -f *.o *.x *~