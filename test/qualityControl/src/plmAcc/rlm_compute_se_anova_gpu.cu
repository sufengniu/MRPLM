/******************************************************************************

  This file contains functions to calculate "rlm_compute_se_anova" by offloading compute intensive 
  parts on the GPU

  Specifically, the parts offloaded are the calculations of XTWX and XTWXinv

  void gpu_XTX(int_t y_cols, int_t y_rows, double* wts, double* XTX)

  This function takes wts as input and generates XTX which the diagonal of XTWXinv.
  wts  - y_rows * y_cols
  All matrices are stored and operated in column major order

  The code is optimized for y_cols >> y_rows

  From wts the symmetric matrix (xtwx) is calculated. xtwx is represented as

  A B
  C D

A: diagonal matrix of size y_cols * y_cols
B: y_cols * (y_rows-1)
C: (y_rows-1) * y_cols
D: (y_rows-1) *(y_rows-1) 

The inverse of the above matrix (xtwx_inv) is

P Q
R S

P: Inv(A) + Inv(A) * B * Inv(D - C * Inv(A) * B) * C * Inv(A)
Q: -Inv(A) * B * Inv(D - C * Inv(A) * B) 
R: -Inv(D - C * Inv(A) * B) * C * Inv(A)
S: Inv(D - C * Inv(A) * B) 

The output of the function is XTX which is the diagonal of xtwx_inv. Hence only the diagonal is calculated for the inverse function

 ******************************************************************************/


#include "rlm_compute_se_anova_gpu.h"

//device memory
static double *d_wts, *d_Ainv, *d_B, *d_C, *d_D_temp, *d_AinvB, *d_CAinvB, *d_CAinv, *d_Q, *d_XTX;

//host memory
static double *h_D, *h_D_temp, *h_CAinvB, *h_S;

/******************************************************************************

  Input - wts [y_rows * y_cols]
  Output - Inv(A) and B [y_cols * (y_rows-1)]
  Since A is a diagonal matrix, it is stored as a vector of size [y_cols]

 ******************************************************************************/
__global__ static void Kernel_Ainv_B(double *d_wts, double* d_Ainv, double *d_B, int_t y_rows, int_t y_cols)
{
	int_t x = blockIdx.x * blockDim.x + threadIdx.x; 
	int_t i;
	double wts_yrow_1, d_wts_elem, d_Ainv_col;

	if(x < y_cols){		
		wts_yrow_1 = d_wts[x*y_rows+y_rows-1];	
		d_Ainv_col = wts_yrow_1;

		for (i=0; i<(y_rows-1); i++){

			d_wts_elem = d_wts[x*y_rows+i];
			d_Ainv_col += d_wts_elem;
			d_B[i*y_cols+x] = d_wts_elem - wts_yrow_1;
		}
		d_Ainv[x] = 1.0/d_Ainv_col;
	}
}

/******************************************************************************

  Input - wts [y_rows * y_cols]
  Output - Dtemp [CUSTOM_VAL * y_rows]

  The calculation of D from wts is divided into two stages to exploit better parallelism.
  This first stage is Kernel_Dtemp which is computed on GPU
  The second stage is computed on the host using the function host_D_from_Dtemp

  CUSTOM_VAL is a manually set parameter

 ******************************************************************************/
__global__ static void Kernel_Dtemp(double *d_wts, double *d_D_temp, int_t y_rows, int_t y_cols)
{
	int_t x = blockIdx.x * blockDim.x + threadIdx.x; 
	int_t y = blockIdx.y * blockDim.y + threadIdx.y; 
	int_t i;

	int_t range = ceil(y_cols/(float)CUSTOM_VAL);
	double d_D_temp_val = 0.0;

	for(i=0 ;i<range; i++){
		if( ((y*range+i)<y_cols) && (x<y_rows) )
			d_D_temp_val+= d_wts[(y*range+i)*y_rows+x];
	}	
	if((x < y_rows) && (y<CUSTOM_VAL)){
		d_D_temp[y*y_rows+x] = d_D_temp_val;
	}		
}



/******************************************************************************

  Input - Ainv and B [y_cols * (y_rows-1)]
  Output - AinvB[y_cols * (y_rows-1)]

  AinvB = (Ainv * B) 
  The multiplication uses the fact that A is a diagonal matrix

 ******************************************************************************/

__global__ static void Kernel_AinvB (double *d_Ainv, double* d_B, double* d_AinvB, int_t y_rows, int_t y_cols)
{
	int_t x = blockIdx.x * blockDim.x + threadIdx.x; 
	int_t i;

	double d_Ainv_elem, d_B_elem;

	if(x < y_cols){

		d_Ainv_elem = d_Ainv[x];

		for (i=0; i<(y_rows-1); i++)
		{
			d_B_elem = d_B[i*(y_cols)+x];
			d_AinvB[i*(y_cols)+x] = d_B_elem*d_Ainv_elem;
		}
	}
}

/******************************************************************************

  Input - A [m * p], B [p * m], Ainv
  Output - d_XTX [m]

  This function multiplies two matrices A [m * p] and B[p * m].
  The diagonal of the resultant matrix is added to Ainv. 
  Ainv is also a diagonal matrix stored as vector.

 ******************************************************************************/

__global__ static void Kernel_XTX (double *A, double* B, double* d_XTX, int_t m, int_t p, double* d_Ainv)
{
	int_t x = blockIdx.x * blockDim.x + threadIdx.x;

	int_t k;
	if(x<m){
		double temp = d_Ainv[x];
		for(k=0; k<p ; k++){
			temp+=A[k*m+x]*B[x*p+k];
		}
		d_XTX[x]=temp;
	}
}

/******************************************************************************

  Input - wts [y_rows * y_cols]
  Output - XTX [y_cols + y_rows - 1]

  XTX is actually a diagonal matrix of size (y_cols + y_rows - 1) * (y_cols + y_rows - 1)

 ******************************************************************************/

void gpu_XTX(int_t y_cols, int_t y_rows, double* wts, double* XTX)
{
	cudaError_t cudaError;

	//Copy wts from Host to device
	cudaError = cudaMemcpy(d_wts, wts, y_cols*y_rows*sizeof(double), cudaMemcpyHostToDevice);
	if(cudaError != cudaSuccess)	printf("\n cudaMemcpyHostToDevice error = %s\n",cudaGetErrorString (cudaError));

	//Launch Kernel to calculate Ainv and B
	{
		dim3 dimGrid(ceil((y_cols)/((float)BLOCK_WIDTH)), 1, 1);
		dim3 dimBlock(BLOCK_WIDTH, 1, 1);	
		Kernel_Ainv_B <<<dimGrid, dimBlock>>> (d_wts, d_Ainv, d_B, y_rows, y_cols);
	}

	//Launch Kernel to calculate D_temp
	//CUSTOM_VAL is manually set
	{
		dim3 dimGrid(ceil(y_rows/(float)TILE_WIDTH), ceil(CUSTOM_VAL/(float)TILE_WIDTH), 1);
		dim3 dimBlock(TILE_WIDTH, TILE_WIDTH, 1);
		Kernel_Dtemp <<<dimGrid, dimBlock>>> (d_wts, d_D_temp, y_rows, y_cols);
	}

	//Copy D_temp from Device to Host
	cudaError = cudaMemcpy(h_D_temp, d_D_temp, CUSTOM_VAL*y_rows*sizeof(double), cudaMemcpyDeviceToHost);
	if(cudaError != cudaSuccess)	printf("\n cudaMemcpyDeviceToHost error = %s\n",cudaGetErrorString (cudaError));

	//Calcuate on Host -> D from D_temp
	host_D_from_Dtemp(h_D, h_D_temp, y_rows, y_cols);

	//Calculate C = trans(B)
	gpu_blas_mtrans(d_B, d_C, (int)y_cols, (int)(y_rows-1));

	//Launch Kernel to calculate AinvB
	{
		dim3 dimGrid(ceil((y_cols)/((float)BLOCK_WIDTH)), 1, 1);
		dim3 dimBlock(BLOCK_WIDTH, 1, 1);	
		Kernel_AinvB <<<dimGrid, dimBlock>>> (d_Ainv, d_B, d_AinvB, y_rows, y_cols);
	}

	//Calculate in Device CAinvB [(y_rows-1) * (y_rows-1)] = C [(y_rows-1) * y_cols] * AinvB [y_cols * (y_rows-1)]	
	gpu_blas_mmul(d_C, d_AinvB, d_CAinvB, (int)(y_rows-1), (int)y_cols, (int)(y_rows-1));

	//Copy CAinvB from Device to Host
	cudaError = cudaMemcpy(h_CAinvB, d_CAinvB, (y_rows-1)*(y_rows-1)*sizeof(double), cudaMemcpyDeviceToHost);
	if(cudaError != cudaSuccess)	printf("\n cudaMemcpyDeviceToHost error = %s\n",cudaGetErrorString (cudaError));

	//CAinvB = (D-C* Inv(A) * B)
	mat_subtract(h_D, h_CAinvB, h_CAinvB, (y_rows-1), (y_rows-1));

	//Calculate Inv(D-C* Inv(A) * B) using Choleski inverse
	//h_CAinvB = (D-C* Inv(A) * B)
	//h_S = Inv(h_CAinvB) = Inv(D-C* Inv(A) * B)
	mat_inverse(h_CAinvB, h_S, (y_rows-1));

	//Copy Inv(D-C* Inv(A) * B) from Host to Device
	cudaError = cudaMemcpy(d_CAinvB, h_S, (y_rows-1)*(y_rows-1)*sizeof(double), cudaMemcpyHostToDevice);
	if(cudaError != cudaSuccess)	printf("\n cudaMemcpyHostToDevice error = %s\n",cudaGetErrorString (cudaError));

	//Calculate on Device Q = AinvB * CAinvB
	//AinvB [y_cols * (y_rows-1)]= A * Inv(B)
	//CAinvB [(y_rows-1) * (y_rows-1)] = Inv(D-C* Inv(A) * B)
	//Q [y_cols * (y_rows-1)] = A * Inv(B) * Inv(D-C* Inv(A) * B)

	gpu_blas_mmul(d_AinvB, d_CAinvB, d_Q, (int)y_cols, (int)(y_rows-1), (int)(y_rows-1));

	//Calculate on Host -> CAinv = trans(AinvB) 
	//AinvB [y_cols * (y_rows-1)] = Inv(A) * B
	//CAinv [(y_rows-1) * y_cols]= C * Inv(A)
	gpu_blas_mtrans(d_AinvB, d_CAinv, y_cols, y_rows-1);

	//Launch Kernel to Calculate XTX
	//Q [y_cols * (y_rows-1)] = A * Inv(B) * Inv(D-C* Inv(A) * B)
	//CAinv [(y_rows-1) * y_cols]= C * Inv(A)
	{
		dim3 dimGrid(ceil((y_cols)/((float)BLOCK_WIDTH)), 1, 1);
		dim3 dimBlock(BLOCK_WIDTH, 1, 1);	
		Kernel_XTX<<<dimGrid, dimBlock>>> (d_Q, d_CAinv, d_XTX, y_cols, (y_rows-1), d_Ainv);
	}

	//Calculate XTX from S
	for(int_t i = 0; i<(y_rows-1); i++){
		XTX[i+y_cols] = h_S[i*(y_rows-1)+i];
	}

	//Copy XTX from Device to Host
	cudaError = cudaMemcpy(XTX, d_XTX, y_cols*sizeof(double), cudaMemcpyDeviceToHost);
	if(cudaError != cudaSuccess)	printf("\n cudaMemcpyDeviceToHost error = %s\n",cudaGetErrorString (cudaError));

	cudaError = cudaGetLastError();
	if(cudaError != cudaSuccess)	printf("\n cudaGetLastError = %s\n",cudaGetErrorString (cudaError));
}

void rlm_compute_se_anova_gpu(int_t y_rows, int_t y_cols, double* resids, double* weights, double* se_estimates, double* varcov, double psi_k){

	//double k1 = psi_k;   /*  was 1.345; */
	//double sumpsi2=0.0;  /* sum of psi(r_i)^2 */
	/*  double sumpsi=0.0; */
	//double sumderivpsi=0.0; /* sum of psi'(r_i) */
	//double Kappa=0.0;      /* A correction factor */
	//double scale=0.0;

	int_t n = y_rows * y_cols;
	int_t p = y_rows + y_cols - 1;

	//double* XTX = new double[p * p];
	double* XTX = new double[p]();

	//double [] W = new double [p*p];
	//double [] work = new double [p*p];

	double RMSEw = 0.0;

	//double vs=0.0,m,varderivpsi=0.0;
	//double [] W_tmp= new double [n];

	alloc_host_mem(y_cols,y_rows);
	alloc_device_mem(y_cols,y_rows);

	for (int_t i = 0; i < n; i++) {
		RMSEw += weights[i] * resids[i] * resids[i];
	}

	RMSEw = sqrt(RMSEw / (double) (n - p));
	//residSE[0] =  RMSEw;

	/***************** GPU offload start ***************/		

#if 0	
	XTWX(y_rows, y_cols, weights, XTX);
	if (y_rows > 1) {
		XTWXinv(y_rows, y_cols, XTX);
	} else {
		for (int_t i = 0; i < p; i++) {
			XTX[i * p + i] = 1.0 / XTX[i * p + i];
		}
	}
#else	
	gpu_XTX(y_cols, y_rows, weights, XTX);
#endif	

	/***************** GPU offload end ***************/

	for (int_t i = 0; i < p; i++) {
		//se_estimates[i] = RMSEw * sqrt(XTX[i * p + i]);
		se_estimates[i] = RMSEw * sqrt(XTX[i]);
	}

	//assuming that (varcov == NULL)
#if 0	
	if (varcov != NULL)
		for (int_t i = 0; i < p; i++)
			for (int_t j = i; j < p; j++)
				varcov[j * p + i] = RMSEw * RMSEw * XTX[j * p + i];
#endif

	delete[] XTX;

	free_device_mem();
	free_host_mem();
}

//Allocate host memory for matrices
static void alloc_host_mem(int_t y_cols, int_t y_rows)
{

#if 0	
	h_D = 		(double*)malloc((y_rows-1)*(y_rows-1)*sizeof(double));
	h_D_temp = 	(double*)malloc(CUSTOM_VAL*y_rows*sizeof(double));
	h_CAinvB = 	(double*)malloc((y_rows-1)*(y_rows-1)*sizeof(double));
	h_S = 		(double*)malloc((y_rows-1)*(y_rows-1)*sizeof(double));	
#else
	cudaError_t cudaError;

	cudaError = cudaHostAlloc((void**)&h_D, 		((y_rows-1)*(y_rows-1)*sizeof(double)), cudaHostAllocDefault);
	if(cudaError != cudaSuccess)	printf("\n cudaHostAlloc h_D error = %s\n",cudaGetErrorString (cudaError));

	cudaError = cudaHostAlloc((void**)&h_D_temp, 	(CUSTOM_VAL*y_rows*sizeof(double)), 	cudaHostAllocDefault);
	if(cudaError != cudaSuccess)	printf("\n cudaHostAlloc h_D_temp error = %s\n",cudaGetErrorString (cudaError));

	cudaError = cudaHostAlloc((void**)&h_CAinvB, 	((y_rows-1)*(y_rows-1)*sizeof(double)), cudaHostAllocDefault);
	if(cudaError != cudaSuccess)	printf("\n cudaHostAlloc h_CAinvB error = %s\n",cudaGetErrorString (cudaError));

	cudaError = cudaHostAlloc((void**)&h_S, 		((y_rows-1)*(y_rows-1)*sizeof(double)), cudaHostAllocDefault);
	if(cudaError != cudaSuccess)	printf("\n cudaHostAlloc h_S error = %s\n",cudaGetErrorString (cudaError));

#endif
}

//Free allocated host memory
static void free_host_mem(void)
{		
#if 0
	free(h_D);	
	free(h_D_temp);
	free(h_CAinvB);
	free(h_S);
#else
	cudaError_t cudaError;

	cudaError = cudaFreeHost(h_D);
	if(cudaError != cudaSuccess)	printf("\n cudaFreeHost h_D error = %s\n",cudaGetErrorString (cudaError));

	cudaError = cudaFreeHost(h_D_temp);
	if(cudaError != cudaSuccess)	printf("\n cudaFreeHost h_D_temp error = %s\n",cudaGetErrorString (cudaError));

	cudaError = cudaFreeHost(h_CAinvB);
	if(cudaError != cudaSuccess)	printf("\n cudaFreeHost h_CAinvB error = %s\n",cudaGetErrorString (cudaError));

	cudaError = cudaFreeHost(h_S);
	if(cudaError != cudaSuccess)	printf("\n cudaFreeHost h_S error = %s\n",cudaGetErrorString (cudaError));

#endif
}

//Allocate device memory for matrices
static void alloc_device_mem(int_t y_cols, int_t y_rows)
{

	cudaError_t cudaError;

	cudaError = cudaMalloc((void**)&d_wts, 			y_cols*y_rows*sizeof(double));
	if(cudaError != cudaSuccess)	printf("\n cudaMalloc d_wts error = %s\n",cudaGetErrorString (cudaError));

	cudaError = cudaMalloc((void**)&d_Ainv, 		y_cols*sizeof(double));
	if(cudaError != cudaSuccess)	printf("\n cudaMalloc d_Ainv error = %s\n",cudaGetErrorString (cudaError));

	cudaError = cudaMalloc((void**)&d_B, 			y_cols*(y_rows-1)*sizeof(double));
	if(cudaError != cudaSuccess)	printf("\n cudaMalloc d_B error = %s\n",cudaGetErrorString (cudaError));

	cudaError = cudaMalloc((void**)&d_C, 			(y_rows-1)*y_cols*sizeof(double));	
	if(cudaError != cudaSuccess)	printf("\n cudaMalloc d_C error = %s\n",cudaGetErrorString (cudaError));

	cudaError = cudaMalloc((void**)&d_D_temp, 		CUSTOM_VAL*y_rows*sizeof(double));
	if(cudaError != cudaSuccess)	printf("\n cudaMalloc d_D_temp error = %s\n",cudaGetErrorString (cudaError));

	cudaError = cudaMalloc((void**)&d_AinvB, 		y_cols*(y_rows-1)*sizeof(double));
	if(cudaError != cudaSuccess)	printf("\n cudaMalloc d_AinvB error = %s\n",cudaGetErrorString (cudaError));

	cudaError = cudaMalloc((void**)&d_CAinvB, 		(y_rows-1)*(y_rows-1)*sizeof(double));
	if(cudaError != cudaSuccess)	printf("\n cudaMalloc d_CAinvB error = %s\n",cudaGetErrorString (cudaError));

	cudaError = cudaMalloc((void**)&d_CAinv, 		(y_rows-1)*y_cols*sizeof(double));
	if(cudaError != cudaSuccess)	printf("\n cudaMalloc d_CAinv error = %s\n",cudaGetErrorString (cudaError));

	cudaError = cudaMalloc((void**)&d_Q, 			y_cols*(y_rows-1)*sizeof(double));
	if(cudaError != cudaSuccess)	printf("\n cudaMalloc d_Q error = %s\n",cudaGetErrorString (cudaError));

	cudaError = cudaMalloc((void**)&d_XTX, 			y_cols*sizeof(double));
	if(cudaError != cudaSuccess)	printf("\n cudaMalloc d_XTX error = %s\n",cudaGetErrorString (cudaError));
}

//Free alocated device memory
static void free_device_mem(void)
{
	cudaError_t cudaError;

	cudaError = cudaFree(d_wts);
	if(cudaError != cudaSuccess)	printf("\n cudaFree d_wts error = %s\n",cudaGetErrorString (cudaError));

	cudaError = cudaFree(d_Ainv);
	if(cudaError != cudaSuccess)	printf("\n cudaFree d_Ainv error = %s\n",cudaGetErrorString (cudaError));

	cudaError = cudaFree(d_B);	
	if(cudaError != cudaSuccess)	printf("\n cudaFree d_B error = %s\n",cudaGetErrorString (cudaError));

	cudaError = cudaFree(d_C);	
	if(cudaError != cudaSuccess)	printf("\n cudaFree d_C error = %s\n",cudaGetErrorString (cudaError));

	cudaError = cudaFree(d_D_temp);
	if(cudaError != cudaSuccess)	printf("\n cudaFree d_D_temp error = %s\n",cudaGetErrorString (cudaError));

	cudaError = cudaFree(d_AinvB);
	if(cudaError != cudaSuccess)	printf("\n cudaFree d_AinvB error = %s\n",cudaGetErrorString (cudaError));

	cudaError = cudaFree(d_CAinvB);
	if(cudaError != cudaSuccess)	printf("\n cudaFree d_CAinvB error = %s\n",cudaGetErrorString (cudaError));

	cudaError = cudaFree(d_CAinv);
	if(cudaError != cudaSuccess)	printf("\n cudaFree d_CAinv error = %s\n",cudaGetErrorString (cudaError));

	cudaError = cudaFree(d_Q);
	if(cudaError != cudaSuccess)	printf("\n cudaFree d_Q error = %s\n",cudaGetErrorString (cudaError));

	cudaError = cudaFree(d_XTX);
	if(cudaError != cudaSuccess)	printf("\n cudaFree d_XTX error = %s\n",cudaGetErrorString (cudaError));
}

