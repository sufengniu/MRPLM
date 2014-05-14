
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "define_t.h"

void gpu_handle(int_t y_cols, int_t y_rows, double* wts, double* y, double* out_beta_gpu);

int Choleski_inverse(double *X, double *Xinv, double *work, int n, int upperonly);

