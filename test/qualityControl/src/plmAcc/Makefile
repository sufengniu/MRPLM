LA_HOME = /home/sniu/lapack-3.5.0
CUDA_HOME = /usr/local/cuda-5.5
CC 	= gcc
NVCC	= nvcc
CLFAGS	:= -Wall -fPIC
LINK	:= -shared -Wl
INC	:= -I${JAVA_HOME}/include -I${JAVA_HOME}/include/linux -I$(LA_HOME)/lapacke/include -I$(CUDA_HOME)/include
TARGET	= jniWrapper 

LIBS	:= -lm -llapacke -llapack -lblas -lgfortran
LIBS	+= -L(LA_HOME) -L/usr/lib -L$(CUDA_HOME)/lib64

ifeq ($(dbg),1)
	CFLAGS += -g
else
	CFLAGS += -O3
endif

OBJS    = cholesky_inv.o \
	  utils.o \
	  wls_acc.o

all:
	$(CC) $(CFLAGS) -o libjniWrapper.so $(INC) *.o

jniWrapper.o: jniWrapper.c jniWrapper.h define_t.h
	$(CC) $(CFLAGS) -c $< -o $@ $(INC)

cholesky_inv.o:cholesky_inv.c cholesky_inv.h
	$(CC) $(CFLAGS) -c $< -o $(INC)

wls_acc.o: wls_acc.cu wls_acc.h define_t.h
	$(NVCC) $(CFLAGS) -c $< -o $@ $(INC)

utils.o: utils.cpp utils.h define_t.h
	$(CC) $(CFLAGS) -c $< -o $@ $(INC)	

clean:
	rm *.o *.so
