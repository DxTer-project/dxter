CC         := g++
LINKER     := $(CC)
CFLAGS	   := -O3 -g -Wall -mfma -march=native -mavx -std=c++0x

HEADERS :=  $(shell find ./ -type f -name '*.h')
SOURCES :=  $(shell find ./ -type f -name '*.cpp')
OBJS := $(patsubst %.cpp, obj/%.o, $(SOURCES))

all: benchmark.x

obj/%.o: %.cpp $(HEADERS)
	@mkdir  -p obj
	$(CC) $(CFLAGS) -c $< -o $@

benchmark.x: $(OBJS) $(HEADERS)
	$(LINKER) $(CFLAGS) $(OBJS) -o $@

clean:
	rm -rf obj/
	rm -f *~
	rm -f *.x
