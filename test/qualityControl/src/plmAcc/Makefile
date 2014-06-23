#dbg	= 1

LA_HOME = /home/sniu/lapack-3.5.0
CUDA_HOME = /opt/cuda-toolkit/5.5.22

CC 	= g++
NVCC	= nvcc

CFLAGS	:= -m64 -fPIC -Wall
NFLAGS 	:= -m64 -arch=sm_35 -Xcompiler '-fPIC'
LINK	= -shared

INC	= -I. -I${JAVA_HOME}/include/ -I${JAVA_HOME}/include/linux -I$(LA_HOME)/lapacke/include/ -I$(CUDA_HOME)/include/
TARGET	= libjniWrapper.so

LIBS    := -L$(LA_HOME) -L$(CUDA_HOME)/lib64
LIBS	+= -lm -llapacke -llapack -lblas -lgfortran -lrt
LIBS	+= -lcudart -lcublas

ifeq ($(dbg),1)
	CFLAGS += -g
	NFLAGS += -g -G
else
	CFLAGS += -O3
	NFLAGS +=
endif

OBJS  = jniWrapper.o \
	matrix_functions.o \
	PLMsummarize_utils.o \
	PLMsummarize_cpu.o \
	rlm_fit_anova_gpu.o \
	rlm_compute_se_anova_gpu.o \
	rlm_fit.o \
	rlm_se.o \
	gpu_utils.o 
		
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(LINK) $^ -o $@ $(LIBS)

jniWrapper.o: affy_qualityControl_jniWrapper.cpp affy_qualityControl_jniWrapper.h define_t.h
	$(CC) $(CFLAGS) -c $< -o $@ $(INC)

matrix_functions.o: matrix_functions.cpp matrix_functions.h define_t.h
	$(CC) $(CFLAGS) -c $< -o $@ $(INC)

PLMsummarize_utils.o: PLMsummarize_utils.cpp PLMsummarize_utils.h define_t.h
	$(CC) $(CFLAGS) -c $< -o $@ $(INC)

PLMsummarize_cpu.o: PLMsummarize_cpu.cpp PLMsummarize_cpu.h define_t.h PLMsummarize_utils.h
	$(CC) $(CFLAGS) -c $< -o $@ $(INC)

rlm_fit_anova_gpu.o: rlm_fit_anova_gpu.cu rlm_fit_anova_gpu.h PLMsummarize_utils.h gpu_utils.h
	$(NVCC) $(NFLAGS) -rdc=true -dc $< -o $@ $(INC)
	$(NVCC) $(NFLAGS) -dlink $@ -o rlm_fit.o

rlm_compute_se_anova_gpu.o: rlm_compute_se_anova_gpu.cu rlm_compute_se_anova_gpu.h PLMsummarize_utils.h gpu_utils.h
	$(NVCC) $(NFLAGS) -rdc=true -dc $< -o $@ $(INC)
	$(NVCC) $(NFLAGS) -dlink $@ -o rlm_se.o

gpu_utils.o: gpu_utils.cpp gpu_utils.h define_t.h
	$(CC) $(CFLAGS) -c $< -o $@ $(INC)

cp:
	cp $(TARGET) ~/hadoop-dist/app/lib
	cp $(TARGET) ~/hadoop-dist/hadoop/lib/
	
clean:
	rm -rf *.o *.so *~ test
