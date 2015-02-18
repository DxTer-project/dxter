CC         := g++
LINKER     := $(CC)
CFLAGS	   := -O3 -g -Wall -std=c++0x -Isrc/ -Isrc/DLA/ -Isrc/tensors -Isrc/LLDLA

HEADERS :=  $(shell find src -type f -name '*.h')
SOURCES :=  $(shell find src -type f -name '*.cpp')
OBJS := $(patsubst %.cpp, %.o, $(SOURCES))
DEPS := $(patsubst %.o, %.dxt_deps, $(OBJS))

all: dxter.x

dxter.x: $(notdir $(OBJS))
	$(LINKER) $(CFLAGS) $(notdir $(OBJS)) -o $@

include $(DEPS)

#obj/%.o: src/%.cpp $(HEADERS)
#	@mkdir  -p obj
#	@mkdir  -p obj/linearization
#	@mkdir  -p obj/DLA
#	@mkdir  -p obj/LLDLA
#	@mkdir  -p obj/tensors
#	$(CC) $(CFLAGS) -c $< -o $@

%.dxt_deps: %.cpp
	bash depends.sh $(CFLAGS) src/ $*.cpp > $@

clean:
	find . -type f -name '*.dxt_deps' -delete
	find . -type f -name '*.o' -delete
	find . -type f -name '*~' -delete
	find . -type f -name '*.x' -delete

open:
	emacs -nw src/*cpp src/*h src/linearization/*cpp src/linearization/*h src/DLA/*cpp src/DLA/*h src/tensors/*cpp src/tensors/*h src/LLDLA/*cpp src/LLDLA/*h makefile

opencpp:
	emacs -nw src/*cpp src/linearization/*cpp src/DLA/*cpp src/tensors/*cpp src/LLDLA/*cpp

openh:
	emacs -nw src/*h src/linearization/*h src/DLA/*h src/tensors/*h src/LLDLA/*h
