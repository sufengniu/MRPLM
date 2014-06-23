
#ifndef DEFINE_T_H_
#define DEFINE_T_H_

#include <stdio.h>
#include <stdlib.h>
#include <new>
#include <math.h>
#include <sys/time.h>
#include "string.h"

#define Y_ROWS 12
#define Y_COLS 10000

#define RANDOM_MAX 1
#define RANDOM_MIN -1

#define BLOCK_WIDTH 256
#define TILE_WIDTH 16

//P_BLOCKS_NUM is the number of blocks that matrix P is divided into
//such that it fits in the GPU
//Each block is of size P_BLOCK_SIZE * P_BLOCK_SIZE
//Use minimum possible value of P_BLOCKS_NUM such that each block fits in GPU

//P_BLOCKS_NUM should be divisble by Y_COLS
//This is done to keep the code simple
#define P_BLOCKS_NUM 10	//should be divisble by Y_COLS
//#define P_BLOCK_SIZE (Y_COLS/P_BLOCKS_NUM)

typedef long int int_t;

#define CPU_COMPUTE

#endif
