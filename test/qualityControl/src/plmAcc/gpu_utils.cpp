/*********************************************************************************
 **
 ** This file contains various functions used by the cuda files to calculate 
 **
 ** rlm_fit_anova and rlm_compute_se_anova
 **
 *********************************************************************************/

#include "gpu_utils.h"

/*
   Calcuate on Host -> D from D_temp
   The matrix D calculates from weights is done in two steps.
   Dtemp is calculated in GPU.
   This function calculates D from Dtemp
   */
void host_D_from_Dtemp(double* h_D, double* h_D_temp, int_t y_rows, int_t y_cols){

	int_t i,j;

	double val_temp = 0.0; //h_D_temp[i][row-1]
	memset(h_D, 0, (y_rows-1)*(y_rows-1)*sizeof(double));

	for(i=0; i<CUSTOM_VAL; i++){
		for(j=0; j<(y_rows-1);j++){
			h_D[j*(y_rows-1)+j] += h_D_temp[i*y_rows+j];
		}
		val_temp += h_D_temp[i*y_rows+y_rows-1];
	}

	for(i=0; i<(y_rows-1); i++){
		for(j=0; j<(y_rows-1); j++){
			h_D[i*(y_rows-1)+j] += val_temp;
		}
	}
}

//Calculate transpose of matrix
//mat is of size y_rows*y_cols stored in column major format
void mat_transpose(double* mat, double* mat_trans, int_t y_rows, int_t y_cols){
	int_t i,j;
	for(j=0; j<y_cols; j++){
		for(i=0; i<y_rows; i++){
			mat_trans[i*y_cols+j] = mat[j*y_rows+i];
		}
	}
}

//C = A-B
//matrices are of size y_rows*y_cols stored in column major format
void mat_subtract(double* A, double* B, double* C, int_t y_rows, int_t y_cols){
	int_t i,j;
	for(j=0; j<y_cols; j++){
		for(i=0; i<y_rows; i++){
			C[j*y_rows +i] = A[j*y_rows +i] - B[j*y_rows +i];
		}
	}
}

//mat_neg = -mat
//mat is of size y_rows*y_cols stored in column major format
void mat_negate(double* mat, double* mat_neg, int_t y_rows, int_t y_cols){
	int_t i,j;
	for(j=0; j<y_cols; j++){
		for(i=0; i<y_rows; i++){
			mat_neg[j*y_rows +i]= -mat[j*y_rows +i];
		}
	}
}

//Calculates inverse of a positive symmetric matrix using Cholesky decomposition
//Uses lapack
void mat_inverse(double* mat, double* mat_inv, int_t y_rows)
{
	double* work = (double*)calloc(y_rows*y_rows,sizeof(double));
	int err = Choleski_inverse(mat, mat_inv, work, (int)y_rows, 0);

	if(err) {
		printf("\n Choleski error code = %d\n", err);
		print_mat(mat, y_rows, y_rows);
		//exit(1);
	}

	free(work);
}

//C(m,n) = A(m,p) * B(p,n)
//Assumes matrices are stored in column major format
void mat_mult(const double *A, const double *B, double *C, const int_t m, const int_t p, const int_t n){

	int_t i,j,k;
	for (j = 0; j<n; j++){
		for (i = 0; i<m; i++){
			C[j*m+i] = 0.0;
			for (k = 0; k<p; k++){
				C[j*m+i]+=A[k*m+i]*B[j*p+k];
			}
		}
	}
}

//print matrix of size y_rows*y_cols stored in column major format
void print_mat(double* A, int_t y_rows, int_t y_cols){
	printf("\n");
	printf("[");
	int_t i,j;
	for (j=0; j<y_cols; j++){
		for(i=0;i<y_rows;i++){
			printf("%e ",A[j*y_rows+i]);
		}
		printf(";");
	}
	printf("]");
	printf("\n");
}

// Multiply the arrays A and B on GPU and save the result in C
// C(m,n) = A(m,k) * B(k,n)
//Assumes matrices are stored in column major format
void gpu_blas_mmul(const double *A, const double *B, double *C, const int m, const int k, const int n) {

	const double alf = 1;
	const double bet = 0;

	int lda=m,ldb=k,ldc=m;

	const double *alpha = &alf;
	const double *beta = &bet;

	// Create a handle for CUBLAS
	cublasHandle_t handle;
	cublasStatus_t cublas_status;	

	cublas_status = cublasCreate(&handle);
	if(cublas_status!=CUBLAS_STATUS_SUCCESS)
		printf("\ngpu_blas_mmul cublasCreate Status = %s\n", cublasGetErrorString(cublas_status));

	// Do the actual multiplication
	cublas_status = cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);

	if(cublas_status!=CUBLAS_STATUS_SUCCESS)
		printf("\ngpu_blas_mmul cublasDgemm Status = %s\n", cublasGetErrorString(cublas_status));

	// Destroy the handle
	cublasDestroy(handle);
}

//Transposes matrix A using cublas
void gpu_blas_mtrans(const double *A, double *A_trans, const int m, const int n) {

	const double alf = 1;
	const double bet = 0;

	int lda=m,ldc=n;
	const double *alpha = &alf;
	const double *beta = &bet;

	// Create a handle for CUBLAS
	cublasHandle_t handle;
	cublasStatus_t cublas_status;

	cublas_status = cublasCreate(&handle);
	if(cublas_status!=CUBLAS_STATUS_SUCCESS)
		printf("\ngpu_blas_mtrans cublasCreate Status = %s\n", cublasGetErrorString(cublas_status));

	// Do the actual transpose
	cublas_status = cublasDgeam(handle, CUBLAS_OP_T, CUBLAS_OP_T, n, m, alpha, A, lda, beta, A, lda, A_trans, ldc);

	if(cublas_status!=CUBLAS_STATUS_SUCCESS)
		printf("\ngpu_blas_mtrans cublasDgeam Status = %s\n", cublasGetErrorString(cublas_status));

	// Destroy the handle
	cublasDestroy(handle);
}

//Cublas Error Messages
const char* cublasGetErrorString(cublasStatus_t status)
{
	switch(status)
	{
		case CUBLAS_STATUS_SUCCESS: 		return "CUBLAS_STATUS_SUCCESS";
		case CUBLAS_STATUS_NOT_INITIALIZED: return "CUBLAS_STATUS_NOT_INITIALIZED";
		case CUBLAS_STATUS_ALLOC_FAILED: 	return "CUBLAS_STATUS_ALLOC_FAILED";
		case CUBLAS_STATUS_INVALID_VALUE: 	return "CUBLAS_STATUS_INVALID_VALUE"; 
		case CUBLAS_STATUS_ARCH_MISMATCH: 	return "CUBLAS_STATUS_ARCH_MISMATCH"; 
		case CUBLAS_STATUS_MAPPING_ERROR: 	return "CUBLAS_STATUS_MAPPING_ERROR";
		case CUBLAS_STATUS_EXECUTION_FAILED:return "CUBLAS_STATUS_EXECUTION_FAILED"; 
		case CUBLAS_STATUS_INTERNAL_ERROR: 	return "CUBLAS_STATUS_INTERNAL_ERROR"; 
	}
	return "unknown error";
}

