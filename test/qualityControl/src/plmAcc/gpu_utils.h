#ifndef _GPU_UTILS_H_
#define _GPU_UTILS_H_

#include "define_t.h"
#include "cublas_v2.h"

const int_t CUSTOM_VAL = 256; //this value is manually tuned

void host_D_from_Dtemp(double* h_D, double* h_D_temp, int_t y_rows, int_t y_cols);

void mat_transpose(double* mat, double* mat_trans, int_t y_cols, int_t y_rows);

void mat_subtract(double* A, double* B, double* C, int_t y_cols, int_t y_rows);

void mat_negate(double* mat, double* mat_neg, int_t y_cols, int_t y_rows);

void mat_inverse(double* mat, double* mat_inv, int_t y_rows);

void mat_mult(const double *A, const double *B, double *C, const int_t m, const int_t p, const int_t n);

void print_mat(double* A, int_t y_cols, int_t y_rows);

void gpu_blas_mmul(const double *A, const double *B, double *C, const int m, const int k, const int n); 

void gpu_blas_mtrans(const double *A, double *A_trans, const int m, const int n);

const char* cublasGetErrorString(cublasStatus_t status);

int Choleski_inverse(double *X, double *Xinv, double *work, int n, int upperonly);

#endif
