INCLUDES=-I/home/sniu/lapack-3.5.0/lapacke/include -I/usr/local/cuda-5.5/include

LIB=-L/home/sniu/lapack-3.5.0/
LIB+=/usr/lib/libgfortran.so
LIB+=/usr/local/cuda-5.5/lib64/libcudart.so.5.5.22 /usr/local/cuda-5.5/lib64/libcublas.so.5.5.22

CC = gcc

test: test.o matrix_functions.o test_cuda.o utils.o
	$(CC) -o test test.o matrix_functions.o test_cuda.o utils.o -lm -llapacke -llapack -lblas -lgfortran $(LIB)

test.o: test.c test.h define_t.h
	$(CC) -O3 -c test.c -o test.o

matrix_functions.o: matrix_functions.c matrix_functions.h
	$(CC) $(INCLUDES) -O3 -c matrix_functions.c -o matrix_functions.o

test_cuda.o: test_cuda.cu test_cuda.h define_t.h
	nvcc  $(INCLUDES) -c test_cuda.cu -o test_cuda.o

utils.o: utils.cpp utils.h define_t.h
	$(CC) -O3 -c utils.cpp -o utils.o

clean: 
	rm -f *.o test
