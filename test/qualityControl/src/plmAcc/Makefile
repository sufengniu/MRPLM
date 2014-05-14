dbg	= 1
LA_HOME = /home/sniu/lapack-3.5.0
#nvcc_loc = $(shell which nvcc)
#CUDA_INSTALL_PATH := $($(patsubst %/bin/nvcc, %, $(nvcc_loc)))

CC 	= gcc
NVCC	= nvcc
CFLAGS	:= -m64 -Wall 
NFLAGS 	:= -m64 
LINK	= 
INC	= -I. -I${JAVA_HOME}/include/ -I${JAVA_HOME}/include/linux/ -I$(LA_HOME)/lapacke/include/ -I${CUDA_INSTALL_PATH}/include/ -I/home/sniu/CBLAS/include/
TARGET	= jniWrapper 

LIBS    := -L$(LA_HOME) -L/home/sniu/CBLAS/lib -L${CUDA_INSTALL_PATH}/lib64
LIBS	+= -lm -llapacke -llapack -lcblas -lblas -lgfortran -lrt
LIBS	+= -lcudart -lcublas

ifeq ($(dbg),1)
	CFLAGS += -g
else
	CFLAGS += -O3
endif

OBJS    = wls_acc.o \
	  utils.o \
	  cholesky_inv.o \
	  jniWrapper.o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(shell source ~/hadoop-dist/bin/setenv.sh)
	$(CC) $(CFLAGS) $(LINK) $^ -o $@ $(LIBS) 

jniWrapper.o: jniWrapper.c jniWrapper.h define_t.h
	$(CC) $(CFLAGS) -c $< -o $@ $(INC)

cholesky_inv.o: cholesky_inv.c cholesky_inv.h
	$(CC) $(CFLAGS) -c $< -o $@ $(INC)

wls_acc.o: wls_acc.cu wls_acc.h utils.h define_t.h 
	$(NVCC) $(NFLAGS) -c $< -o $@ $(INC)

utils.o: utils.c utils.h define_t.h
	$(CC) $(CFLAGS) -c $< -o $@ $(INC)	

clean:
	rm -rf *.o *.so
