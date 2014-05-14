#ifndef	WLS_ACC_H_
#define WLS_ACC_H_

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <cublas_v2.h>

#include "define_t.h"
#include "utils.h"

void wls_acc(int_t y_cols, int_t y_rows, double* wts, double* y, double* out_beta_gpu);

void gpu_blas_mmul(const double *A, const double *B, double *C, const int m, const int k, const int n, double alf, double bet) ;

void gpu_blas_mtrans(const double *A, double *A_trans, const int m, const int n, double alf, double bet); 

void alloc_host_mem(int_t y_cols, int_t y_rows);

void free_host_mem(void);

void alloc_device_mem(int_t y_cols, int_t y_rows);

void free_device_mem(void);

const char* cublasGetErrorString(cublasStatus_t status);

extern void gpu_XTWY(int_t y_rows, int_t y_cols, double *wts,double *y, double *xtwy);

extern void mat_transpose(double* mat, double* mat_trans, int_t y_cols, int_t y_rows);

extern void mat_subtract(double* A, double* B, double* C, int_t y_cols, int_t y_rows);

extern void mat_negate(double* mat, double* mat_neg, int_t y_cols, int_t y_rows);

extern void mat_inverse(double* mat, double* mat_inv, int_t y_rows);

extern void mat_mult(const double *A, const double *B, double *C, const int_t m, const int_t p, const int_t n);

extern void host_D_from_Dtemp(double* h_D, double* h_D_temp, int_t y_cols, int_t y_rows);

extern void print_mat(double* A, int_t y_cols, int_t y_rows);

#endif
