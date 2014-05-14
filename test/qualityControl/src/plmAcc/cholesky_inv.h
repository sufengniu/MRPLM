#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>

int Choleski_decompose(double *X, double *L, int n, int lapack);
int Choleski_2_inverse(double *X, double *Xinv, int n,int upperonly, int lapack);
int Choleski_inverse(double *X, double *Xinv, double *work, int n, int upperonly);

