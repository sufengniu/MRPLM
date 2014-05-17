#include "wls_acc.h"
#include "utils.h"

//device memory
double *d_wts, *d_xtwy, *d_out_beta_gpu, *d_Ainv, *d_B, *d_C, *d_D_temp, *d_AinvB, *d_CAinvB, *d_CAinv, *d_Q, *d_P_block, *d_P_trans_block;

//host memory
double *h_Ainv, *h_D, *h_D_temp, *h_CAinvB, *h_Q, *h_S;
double *xtwy;

//Kernel to calculate A, Ainv and B
__global__ void Kernel_Ainv_B(double *d_wts, double* d_Ainv, double *d_B, int_t y_rows, int_t y_cols)
{
	int_t col = blockIdx.x * blockDim.x + threadIdx.x; 
	int_t i;
	double wts_yrow_1, d_wts_elem, d_Ainv_col;

	if(col < y_cols){		
		wts_yrow_1 = d_wts[col*y_rows+y_rows-1];	
		d_Ainv_col = wts_yrow_1;

		for (i=0; i<(y_rows-1); i++){

			d_wts_elem = d_wts[col*y_rows+i];
			d_Ainv_col += d_wts_elem;
			d_B[i*y_cols+col] = d_wts_elem - wts_yrow_1;
		}
		d_Ainv[col] = 1.0/d_Ainv_col;
	}
}

//Kernel to Calculate D_temp
__global__ void Kernel_Dtemp(double *d_wts, double *d_D_temp, int_t y_rows, int_t y_cols)
{
	int_t col = blockIdx.x * blockDim.x + threadIdx.x; 
	int_t row = blockIdx.y * blockDim.y + threadIdx.y; 
	int_t i;

	//
	int_t x_range = (int_t)(y_cols/CUSTOM_VAL) + 1;
	double d_D_temp_val = 0.0;

	for(i=0 ;i<x_range; i++){
		if( ((col*x_range+i)<y_cols) && (row<y_rows) )
			d_D_temp_val+= d_wts[(col*x_range+i)*y_rows+row];
	}	
	if(row < y_rows){
		d_D_temp[col*y_rows+row] = d_D_temp_val;
	}		
}

__global__ void Kernel_AinvB (double *d_Ainv, double* d_B, double* d_AinvB, int_t y_rows, int_t y_cols)
{
	int_t col = blockIdx.x * blockDim.x + threadIdx.x; 
	int_t i;

	double d_Ainv_elem, d_B_elem;

	if(col < y_cols){

		d_Ainv_elem = d_Ainv[col];

		for (i=0; i<(y_rows-1); i++)
		{
			d_B_elem = d_B[i*(y_cols)+col];
			d_AinvB[i*(y_cols)+col] = d_B_elem*d_Ainv_elem;
		}
	}
}

__global__ void Kernel_outbeta_P (double *d_P_block, double* d_xtwy, double* d_out_beta, int_t size, 
		int_t out_index_row, int_t out_index_col)
{
	int_t row = blockIdx.x * blockDim.x + threadIdx.x; 
	int_t i;

	double d_out_beta_elem;

	if(row<size){
		d_out_beta_elem = d_out_beta[out_index_row+row];
		for(i=0;i<size; i++){
			d_out_beta_elem -= d_P_block[i*size+row]*d_xtwy[out_index_col+i];
		}	
		d_out_beta[out_index_row+row] = d_out_beta_elem;
	}
}

__global__ void Kernel_outbeta_P_trans (double *d_P_block, double* d_xtwy, double* d_out_beta, int_t size, 
		int_t out_index_row, int_t out_index_col)
{
	int_t col = blockIdx.x * blockDim.x + threadIdx.x;
	int_t i;

	double d_out_beta_elem;

	if(col<size){
		d_out_beta_elem = d_out_beta[out_index_col+col];
		for(i=0;i<size; i++){
			d_out_beta_elem -= d_P_block[col*size+i]*d_xtwy[out_index_row+i];
		}		
		d_out_beta[out_index_col+col] = d_out_beta_elem;
	}
}

extern "C" void wls_gpu(int_t y_cols, int_t y_rows, double* wts, double* y, double* out_beta_gpu)
{
	struct timeval start, end;
	long utime;	

	alloc_host_mem(y_cols,y_rows);
	alloc_device_mem(y_cols,y_rows);

	gpu_XTWY(y_rows, y_cols, wts, y, xtwy);

	//Copy wts from Host to device
	cudaMemcpy(d_wts, wts, y_cols*y_rows*sizeof(double), cudaMemcpyHostToDevice);

	{
		dim3 dimGrid(ceil((y_cols)/((float)BLOCK_WIDTH)), 1, 1);
		dim3 dimBlock(BLOCK_WIDTH, 1, 1);	

		//Launch Kernel to calculate Ainv and B
		Kernel_Ainv_B <<<dimGrid, dimBlock>>> (d_wts, d_Ainv, d_B, y_rows, y_cols);
	}

	//Copy Ainv from Device to Host
	cudaMemcpy(h_Ainv, d_Ainv, y_cols*sizeof(double), cudaMemcpyDeviceToHost);

	//sqrt(70,000) = 265
	//launch grid with 256*16
	//this is hard-coded need to improve
	{
		dim3 dimGrid((CUSTOM_VAL/TILE_WIDTH), 1, 1);
		dim3 dimBlock(TILE_WIDTH, TILE_WIDTH, 1);

		//Launch Kernel to calculate D_temp
		Kernel_Dtemp <<<dimGrid, dimBlock>>> (d_wts, d_D_temp, y_rows, y_cols);
	}

	//Copy D_temp from Device to Host
	cudaMemcpy(h_D_temp, d_D_temp, CUSTOM_VAL*y_rows*sizeof(double), cudaMemcpyDeviceToHost);

	//Calcuate on Host -> D from D_temp
	host_D_from_Dtemp(h_D, h_D_temp, y_rows, y_cols);

	//Inverse start
	gpu_blas_mtrans(d_B, d_C, (int)y_cols, (int)(y_rows-1), 1.0, 0.0);

	{
		dim3 dimGrid(ceil((y_cols)/((float)BLOCK_WIDTH)), 1, 1);
		dim3 dimBlock(BLOCK_WIDTH, 1, 1);	

		//Launch Kernel to calculate AinvB
		Kernel_AinvB <<<dimGrid, dimBlock>>> (d_Ainv, d_B, d_AinvB, y_rows, y_cols);
	}

	//Calculate in Device CAinvB(y_rows-1,y_rows-1) -> C(y_rows-1, y_cols) * AinvB(y_cols, y_rows-1)	
	// C(m,n) = A(m,k) * B(k,n)
	//gpu_blas_mmul(const double *A, const double *B, double *C, const int m, const int k, const int n, double alf, double bet) 
	gpu_blas_mmul(d_C, d_AinvB, d_CAinvB, (int)(y_rows-1), (int)y_cols, (int)(y_rows-1), 1.0, 0.0);

	//Copy CAinvB from Device to Host
	cudaMemcpy(h_CAinvB, d_CAinvB, (y_rows-1)*(y_rows-1)*sizeof(double), cudaMemcpyDeviceToHost);

	//CAinvB = (D-CAinvB)
	mat_subtract(h_D, h_CAinvB, h_CAinvB, (y_rows-1), (y_rows-1));

	//Calculate (D-CAinvB)inv using Choleski inverse
	//Q = (D-CAinvB)inv; D  -> (y_rows-1)x(y_rows-1)
	mat_inverse(h_CAinvB, h_S, (y_rows-1));

	//h_CAinvB = -(D-CAinvB)inv
	mat_negate(h_S, h_CAinvB, (y_rows-1), (y_rows-1));

	//Copy -(D-CAinvB)inv from Host to Device
	cudaMemcpy(d_CAinvB, h_CAinvB, (y_rows-1)*(y_rows-1)*sizeof(double), cudaMemcpyHostToDevice);

	//Calculate on Device Q (y_cols, y_rows-1) = AinvB (y_cols, y_rows-1) *[- (D-CAinvB)inv] (y_rows-1, y_rows-1)	
	// C(m,n) = A(m,k) * B(k,n)
	//gpu_blas_mmul(const double *A, const double *B, double *C, const int m, const int k, const int n, double alf, double bet) 
	gpu_blas_mmul(d_AinvB, d_CAinvB, d_Q, (int)y_cols, (int)(y_rows-1), (int)(y_rows-1), 1.0, 0.0);

	//Copy d_B [-AinvC *(D-BAinvC)inv] from Device to Host
	cudaMemcpy(h_Q, d_Q, y_cols*(y_rows-1)*sizeof(double), cudaMemcpyDeviceToHost);

	//Calculate on Host -> CAinv = trans(AinvB) 
	gpu_blas_mtrans(d_AinvB, d_CAinv, y_cols, y_rows-1, 1.0, 0.0);

	//Copy out_beta_gpu from Host to Device
	cudaMemcpy(d_out_beta_gpu, out_beta_gpu, (y_rows+y_cols)*sizeof(double), cudaMemcpyHostToDevice);

	//Copy xtwy from Host to Device
	cudaMemcpy(d_xtwy, xtwy, (y_rows+y_cols)*sizeof(double), cudaMemcpyHostToDevice);

	int_t i,j,k;

	double* temp;
	temp = (double*)calloc(y_cols*(y_rows-1), sizeof(double));	
	double** X_dptr = (double**)calloc(BLOCKS_NUM, sizeof(double*));
	for(i=0;i<BLOCKS_NUM;i++){
		X_dptr[i] = temp+i*(BLOCK_SIZE*(y_rows-1));
	}

	for(k=0; k<BLOCKS_NUM; k++){
		for(j=0 ;j<(y_rows-1); j++){
			memcpy(X_dptr[k]+j*BLOCK_SIZE, h_Q+j*y_cols+k*BLOCK_SIZE, BLOCK_SIZE*sizeof(double));
		}
	}

	double **d_X_dptr;
	double **d_Y_dptr;		

	d_X_dptr = (double**)calloc(BLOCKS_NUM, sizeof(double*));
	d_Y_dptr = (double**)calloc(BLOCKS_NUM, sizeof(double*));

	cudaMalloc((void**)&temp, y_cols*(y_rows-1)*sizeof(double));	
	for(i=0;i<BLOCKS_NUM;i++){
		d_X_dptr[i] = temp+i*(BLOCK_SIZE*(y_rows-1));
	}

	for(i=0;i<BLOCKS_NUM;i++){
		d_Y_dptr[i] = d_CAinv+i*((y_rows-1)*BLOCK_SIZE);
	}

	cudaMemcpy(d_X_dptr[0], X_dptr[0], y_cols*(y_rows-1)*sizeof(double), cudaMemcpyHostToDevice);

	//non diag blocks
	dim3 dimGrid(ceil((BLOCK_SIZE)/((float)BLOCK_WIDTH)), 1, 1);
	dim3 dimBlock(BLOCK_WIDTH, 1, 1);

#if 0		
	for(j=0; j<BLOCKS_NUM; j++){
		for(i=0; i<j; i++){	

			gettimeofday(&start, NULL);
			// C(m,n) = A(m,k) * B(k,n)
			//gpu_blas_mmul(const double *A, const double *B, double *C, const int m, const int k, const int n) 
			gpu_blas_mmul(d_X_dptr[i], d_Y_dptr[j], d_P_block, (int)BLOCK_SIZE, (int)(y_rows-1), (int)BLOCK_SIZE, 1.0, 0.0);

			//__global__ void Kernel_outbeta_P (double *d_P_block, double* d_xtwy, double* d_out_beta, int_t size, 
			//					int_t out_index_row, int_t out_index_col)					

			cudaDeviceSynchronize();

			gettimeofday(&end, NULL);
			seconds  = end.tv_sec  - start.tv_sec;
			useconds = end.tv_usec - start.tv_usec;
			utime = ((seconds) * 1000000 + useconds);
			printf("\nTime GPU gpu_blas_mmul = %ld us\n",utime);

			gettimeofday(&start, NULL);
			//Launch Kernel to calculate out_beta for P
			Kernel_outbeta_P <<<dimGrid, dimBlock>>> (d_P_block, d_xtwy, d_out_beta_gpu, (int_t)BLOCK_SIZE,
					i*BLOCK_SIZE, j*BLOCK_SIZE);
			cudaDeviceSynchronize();

			gettimeofday(&end, NULL);
			seconds  = end.tv_sec  - start.tv_sec;
			useconds = end.tv_usec - start.tv_usec;
			utime = ((seconds) * 1000000 + useconds);
			printf("\nTime GPU Kernel_outbeta_P = %ld us\n",utime);

			gettimeofday(&start, NULL);
			//Launch Kernel to calculate out_beta for P trans
			Kernel_outbeta_P_trans<<<dimGrid, dimBlock>>> (d_P_block, d_xtwy, d_out_beta_gpu, (int_t)BLOCK_SIZE,
					i*BLOCK_SIZE, j*BLOCK_SIZE);
			cudaDeviceSynchronize();

			gettimeofday(&end, NULL);
			seconds  = end.tv_sec  - start.tv_sec;
			useconds = end.tv_usec - start.tv_usec;
			utime = ((seconds) * 1000000 + useconds);
			printf("\nTime GPU Kernel_outbeta_P_trans = %ld us\n",utime);

		}
	}

	//diag blocks
	for(j=0,i=0; j<BLOCKS_NUM; j++,i++){

		// C(m,n) = A(m,k) * B(k,n)
		//gpu_blas_mmul(const double *A, const double *B, double *C, const int m, const int k, const int n)
		gpu_blas_mmul(d_X_dptr[i], d_Y_dptr[j], d_P_block, (int)BLOCK_SIZE, (int)(y_rows-1), (int)BLOCK_SIZE, 1.0, 0.0);	

		Kernel_outbeta_P <<<dimGrid, dimBlock>>> (d_P_block, d_xtwy, d_out_beta_gpu, (int_t)BLOCK_SIZE,
				i*BLOCK_SIZE, j*BLOCK_SIZE);
	}

#else
	for(j=0; j<BLOCKS_NUM; j++){
		for(i=0; i<BLOCKS_NUM; i++){ 

			// C(m,n) = A(m,k) * B(k,n)
			//gpu_blas_mmul(const double *A, const double *B, double *C, const int m, const int k, const int n) 
			gpu_blas_mmul(d_X_dptr[i], d_Y_dptr[j], d_P_block, (int)BLOCK_SIZE, (int)(y_rows-1), (int)BLOCK_SIZE, 1.0, 0.0);

			//gpu_blas_mtrans(d_P_block, d_P_trans_block, (int)BLOCK_SIZE, (int)BLOCK_SIZE, 1.0, 0.0);

			//Launch Kernel to calculate out_beta for P
			//__global__ void Kernel_outbeta_P (double *d_P_block, double* d_xtwy, double* d_out_beta, int_t size, 
			//											int_t out_index_row, int_t out_index_col)					
			Kernel_outbeta_P <<<dimGrid, dimBlock>>> (d_P_block, d_xtwy, d_out_beta_gpu, (int_t)BLOCK_SIZE,
					i*BLOCK_SIZE, j*BLOCK_SIZE);			
			//Kernel_outbeta_P <<<dimGrid, dimBlock>>> (d_P_trans_block, d_xtwy, d_out_beta_gpu, (int_t)BLOCK_SIZE,
			//					j*BLOCK_SIZE, i*BLOCK_SIZE);			
		}
	}

	//diag blocks
	/*
	   for(j=0,i=0; j<BLOCKS_NUM; j++,i++){

// C(m,n) = A(m,k) * B(k,n)
/gpu_blas_mmul(const double *A, const double *B, double *C, const int m, const int k, const int n)
gpu_blas_mmul(d_X_dptr[i], d_Y_dptr[j], d_P_block, (int)BLOCK_SIZE, (int)(y_rows-1), (int)BLOCK_SIZE, 1.0, 0.0);	

Kernel_outbeta_P <<<dimGrid, dimBlock>>> (d_P_block, d_xtwy, d_out_beta_gpu, (int_t)BLOCK_SIZE,
i*BLOCK_SIZE, j*BLOCK_SIZE);
}
	 */
#endif

	//Copy d_out_beta_gpu from Device to Host
	cudaMemcpy(out_beta_gpu, d_out_beta_gpu, (y_rows+y_cols)*sizeof(double), cudaMemcpyDeviceToHost);

	for(j=0; j<y_cols; j++){
		out_beta_gpu[j] += h_Ainv[j]*xtwy[j];
	}

	for(j=0; j<(y_rows-1); j++){
		for(i=0; i<y_cols; i++){
			out_beta_gpu[i] += h_Q[j*y_cols+i]*xtwy[y_cols+j];
			out_beta_gpu[y_cols+j] += h_Q[j*y_cols+i]*xtwy[i];
		}
	}

	for(j=0; j<(y_rows-1); j++){
		for(i=0; i<(y_rows-1); i++){
			out_beta_gpu[y_cols+i] += h_S[j*(y_rows-1)+i]*xtwy[y_cols+j];
		}
	}

	//done

	//free GPU memory
	cudaFree(*d_X_dptr);

	//Free Host Memory
	free(*X_dptr);		
	free(X_dptr);	
	free(d_X_dptr);
	free(d_Y_dptr);

	free_device_mem();
	free_host_mem();
}

// Multiply the arrays A and B on GPU and save the result in C
// C(m,n) = A(m,k) * B(k,n)
void gpu_blas_mmul(const double *A, const double *B, double *C, const int m, const int k, const int n,
		double alf, double bet) {

	//const double alf = 1;
	//const double bet = 0;

	int lda=m,ldb=k,ldc=m;
	const double *alpha = &alf;
	const double *beta = &bet;

	// Create a handle for CUBLAS
	cublasHandle_t handle;
	cublasStatus_t cublas_status;

	cublas_status = cublasCreate(&handle);
	if(cublas_status!=CUBLAS_STATUS_SUCCESS)
		printf("\nCublas Status = %s\n", cublasGetErrorString(cublas_status));

	// Do the actual multiplication
	cublas_status = cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);

	if(cublas_status!=CUBLAS_STATUS_SUCCESS)
		printf("\nCublas Status = %s\n", cublasGetErrorString(cublas_status));

	// Destroy the handle
	cublasDestroy(handle);
}


void gpu_blas_mtrans(const double *A, double *A_trans, const int m, const int n,	double alf, double bet) {

	//const double alf = 1;
	//const double bet = 0;

	int lda=m,ldc=n;
	const double *alpha = &alf;
	const double *beta = &bet;

	// Create a handle for CUBLAS
	cublasHandle_t handle;
	cublasStatus_t cublas_status;

	cublas_status = cublasCreate(&handle);
	if(cublas_status!=CUBLAS_STATUS_SUCCESS)
		printf("\nCublas Status = %s\n", cublasGetErrorString(cublas_status));

	// Do the actual transpose
	//cublasStatus_t cublasDgeam(cublasHandle_t handle, cublasOperation_t transa, cublasOperation_t transb, int m, int n, const double *alpha, 
	// const double *A, int lda, const double *beta, const double *B, int ldb, double *C, int ldc)
	cublas_status = cublasDgeam(handle, CUBLAS_OP_T, CUBLAS_OP_T, n, m, alpha, A, lda, beta, A, lda, A_trans, ldc);

	if(cublas_status!=CUBLAS_STATUS_SUCCESS)
		printf("\nCublas Status = %s\n", cublasGetErrorString(cublas_status));

	// Destroy the handle
	cublasDestroy(handle);
}

void alloc_host_mem(int_t y_cols, int_t y_rows)
{
	//Allocate host memory for matrices
#if 0	
	h_Ainv = 	(double*)malloc(y_cols*sizeof(double));
	h_D = 		(double*)malloc((y_rows-1)*(y_rows-1)*sizeof(double));
	h_D_temp = 	(double*)malloc(CUSTOM_VAL*y_rows*sizeof(double));
	h_CAinvB = 	(double*)malloc((y_rows-1)*(y_rows-1)*sizeof(double));
	h_Q = 		(double*)malloc(y_cols*(y_rows-1)*sizeof(double));
	h_S = 		(double*)malloc((y_rows-1)*(y_rows-1)*sizeof(double));	
	xtwy = 		(double*)malloc((y_rows+y_cols)*sizeof(double));	
#else
	cudaHostAlloc((void**)&h_Ainv, 		(y_cols*sizeof(double)), 				cudaHostAllocDefault);
	cudaHostAlloc((void**)&h_D, 		((y_rows-1)*(y_rows-1)*sizeof(double)), cudaHostAllocDefault);
	cudaHostAlloc((void**)&h_D_temp, 	(CUSTOM_VAL*y_rows*sizeof(double)), 	cudaHostAllocDefault);
	cudaHostAlloc((void**)&h_CAinvB, 	((y_rows-1)*(y_rows-1)*sizeof(double)), cudaHostAllocDefault);
	cudaHostAlloc((void**)&h_Q, 		(y_cols*(y_rows-1)*sizeof(double)), 	cudaHostAllocDefault);
	cudaHostAlloc((void**)&h_S, 		((y_rows-1)*(y_rows-1)*sizeof(double)), cudaHostAllocDefault);
	cudaHostAlloc((void**)&xtwy, 		((y_rows+y_cols)*sizeof(double)), 		cudaHostAllocDefault);
#endif
}

void free_host_mem(void)
{		
#if 0
	free(h_Ainv);
	free(h_D);	
	free(h_D_temp);
	free(h_CAinvB);
	free(h_Q);
	free(h_S);
	free(xtwy);
#else
	cudaFreeHost(h_Ainv);
	cudaFreeHost(h_D);	
	cudaFreeHost(h_D_temp);
	cudaFreeHost(h_CAinvB);
	cudaFreeHost(h_Q);
	cudaFreeHost(h_S);
	cudaFreeHost(xtwy);
#endif
}

void alloc_device_mem(int_t y_cols, int_t y_rows)
{

	//Allocate device memory for matrices
	cudaMalloc((void**)&d_wts, 			y_cols*y_rows*sizeof(double));
	cudaMalloc((void**)&d_Ainv, 		y_cols*sizeof(double));
	cudaMalloc((void**)&d_B, 			y_cols*(y_rows-1)*sizeof(double));
	cudaMalloc((void**)&d_C, 			(y_rows-1)*y_cols*sizeof(double));	
	cudaMalloc((void**)&d_D_temp, 		CUSTOM_VAL*y_rows*sizeof(double));
	cudaMalloc((void**)&d_AinvB, 		y_cols*(y_rows-1)*sizeof(double));
	cudaMalloc((void**)&d_CAinvB, 		(y_rows-1)*(y_rows-1)*sizeof(double));
	cudaMalloc((void**)&d_CAinv, 		(y_rows-1)*y_cols*sizeof(double));
	cudaMalloc((void**)&d_Q, 			y_cols*(y_rows-1)*sizeof(double));
	cudaMalloc((void**)&d_out_beta_gpu, (y_rows+y_cols)*sizeof(double));
	cudaMalloc((void**)&d_xtwy, 		(y_rows+y_cols)*sizeof(double));
	cudaMalloc((void**)&d_P_block, 		BLOCK_SIZE*BLOCK_SIZE*sizeof(double));
	//cudaMalloc((void**)&d_P_trans_block,BLOCK_SIZE*BLOCK_SIZE*sizeof(double));
}

void free_device_mem(void)
{
	cudaFree(d_wts);
	cudaFree(d_Ainv);
	cudaFree(d_B);	
	cudaFree(d_C);	
	cudaFree(d_D_temp);
	cudaFree(d_AinvB);
	cudaFree(d_CAinvB);
	cudaFree(d_CAinv);
	cudaFree(d_Q);
	cudaFree(d_out_beta_gpu);
	cudaFree(d_xtwy);
	cudaFree(d_P_block);
	//cudaFree(d_P_trans_block);
}

const char* cublasGetErrorString(cublasStatus_t status)
{
	switch(status)
	{
		case CUBLAS_STATUS_SUCCESS: return "CUBLAS_STATUS_SUCCESS";
		case CUBLAS_STATUS_NOT_INITIALIZED: return "CUBLAS_STATUS_NOT_INITIALIZED";
		case CUBLAS_STATUS_ALLOC_FAILED: return "CUBLAS_STATUS_ALLOC_FAILED";
		case CUBLAS_STATUS_INVALID_VALUE: return "CUBLAS_STATUS_INVALID_VALUE"; 
		case CUBLAS_STATUS_ARCH_MISMATCH: return "CUBLAS_STATUS_ARCH_MISMATCH"; 
		case CUBLAS_STATUS_MAPPING_ERROR: return "CUBLAS_STATUS_MAPPING_ERROR";
		case CUBLAS_STATUS_EXECUTION_FAILED: return "CUBLAS_STATUS_EXECUTION_FAILED"; 
		case CUBLAS_STATUS_INTERNAL_ERROR: return "CUBLAS_STATUS_INTERNAL_ERROR"; 
	}
	return "unknown error";
}
