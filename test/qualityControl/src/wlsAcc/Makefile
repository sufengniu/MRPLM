dbg	= 1

LA_HOME = /home/sniu/lapack-3.5.0
#nvcc_loc = $(shell which nvcc)
#CUDA_INSTALL_PATH := $($(patsubst %/bin/nvcc, %, $(nvcc_loc)))

CC 	= gcc
NVCC	= nvcc

CFLAGS	:= -m64 -fPIC 
NFLAGS 	:= -m64 -arch=sm_20 -Xcompiler '-fPIC' 
LINK	= -shared

INC	= -I. -I${JAVA_HOME}/include/ -I${JAVA_HOME}/include/linux/ -I$(LA_HOME)/lapacke/include/ -I${CUDA_INSTALL_PATH}/include/ -I/home/sniu/CBLAS/include/
TARGET	= libjniWrapper.so

LIBS    := -L$(LA_HOME) -L/home/sniu/CBLAS/lib -L${CUDA_INSTALL_PATH}/lib64
LIBS	+= -lm -llapacke -llapack -lblas -lgfortran -lrt
LIBS	+= -lcudart -lcublas

ifeq ($(dbg),1)
	CFLAGS += -g
	NFLAGS += -g -G
else
	CFLAGS += -O3
	NFLAGS +=
endif

OBJS    = jniWrapper.o \
	  cholesky_inv.o \
	  wls.o \
	  wls_acc.o \
	  wls_link.o \
	  utils.o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(LINK) $^ -o $@ $(LIBS) 

jniWrapper.o: affy_qualityControl_jniWrapper.c affy_qualityControl_jniWrapper.h define_t.h
	$(CC) $(CFLAGS) -c $< -o $@ $(INC)

cholesky_inv.o: cholesky_inv.c cholesky_inv.h
	$(CC) $(CFLAGS) -c $< -o $@ $(INC)

wls_acc.o: wls_acc.cu wls_acc.h utils.h define_t.h 
	$(NVCC) $(NFLAGS) -rdc=true -dc $< -o $@ $(INC)
	$(NVCC) $(NFLAGS) -dlink $@ -o wls_link.o 

wls.o: wls.c define_t.h
	$(CC) $(CFLAGS) -c $< -o $@ $(INC)

utils.o: utils.cpp utils.h define_t.h
	$(CC) $(CFLAGS) -c $< -o $@ $(INC)	

cp:
	cp $(TARGET) ~/hadoop-dist/app/lib
	cp $(TARGET) ~/hadoop-dist/hadoop/lib/

clean:
	rm -rf *.o *.so *~ test
	rm ~/hadoop-dist/app/lib/$(TARGET) ~/hadoop-dist/hadoop/lib/$(TARGET)
