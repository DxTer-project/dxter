# indicate where the object files are to be created
CC         := gcc
LINKER     := $(CC)
CFLAGS	   := -O3 -msse3 -mfpmath=sse -Wall -I../kernels

HEADERS = row_stride_tests.h ../kernels/utils.h gen_stride_tests.h
SOURCES = all_tests.c row_stride_tests.c ../kernels/utils.c gen_stride_tests.c
#HEADERS +=  $(shell find ../kernels -type f -name '*.h')
#SOURCES +=  $(shell find ../kernels -type f -name '*.c')

#HEADERS :=  utils.h ../kernels/gen_stride_lldla_primitives.h gen_stride_tests.h row_stride_tests.h row_stride_lldla_primitives.h
#SOURCES :=  utils.c all_tests.c gen_stride_lldla_primitives.c gen_stride_tests.c row_stride_tests.c row_stride_lldla_primitives.c
OBJS := $(patsubst %.c, %.o, $(SOURCES))

all_tests : $(OBJS)
	$(CC) $(CFLAGS) -I../kernels $(OBJS) -o all_tests 

clean :
	rm -f $(OBJS) all_tests