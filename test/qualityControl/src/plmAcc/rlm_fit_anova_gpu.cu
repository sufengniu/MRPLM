/******************************************************************************

This file contains functions to calculate "rlm_fit_anova" by offloading compute intensive 
parts on the GPU

Specifically, the parts offloaded are the calculations of XTWX and XTWXinv and 
subsequent calculation of out_beta

gpu_out_beta(int_t y_cols, int_t y_rows, double* wts, double* xtwy, double* out_beta_gpu)

This function takes wts and xtwy as inputs. All matrices are stored and operated in column major order
wts  - y_rows * y_cols
xtwy - vector of length y_cols + y_rows - 1

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

The inverse matrix xtwx_inv and xtwy are then used to calculate out_beta

******************************************************************************/

#include "rlm_fit_anova_gpu.h"

//device memory
static double *d_wts, *d_xtwy, *d_out_beta_gpu, *d_Ainv, *d_B, *d_C, *d_D_temp, *d_AinvB, *d_CAinvB, *d_CAinv, *d_Q, *d_P_block;

//host memory
static double *h_Ainv, *h_D, *h_D_temp, *h_CAinvB, *h_Q, *h_S;

int P_BLOCK_SIZE;

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

Input - P_block[size*size], xtwy, size, out_index_row, out_index_col
Output - out_beta

The inverse matrix P is divided into parts since the entire matrix does not fit in the GPU
P_block is one part of the matrix

This function takes P_block and xtwy as input and performs the corresponding calculations
to evaluate out_beta

P_block is of dimensions size*size

out_index_row and out_index_col mark the relative position of P_block in the actual P matrix.
This also determines which part of out_beta is evaluated

******************************************************************************/

__global__ static void Kernel_outbeta_P (double *d_P_block, double* d_xtwy, double* d_out_beta, int_t size, int_t out_index_row, int_t out_index_col)
{
	int_t row = blockIdx.x * blockDim.x + threadIdx.x; 
	int_t i;

	double d_out_beta_elem;
	
	if(row<size){
		d_out_beta_elem = d_out_beta[out_index_row+row];
		for(i=0; i<size; i++){
			d_out_beta_elem -= d_P_block[i*size+row]*d_xtwy[out_index_col+i];
		}	
		d_out_beta[out_index_row+row] = d_out_beta_elem;
	}
}

/******************************************************************************

Input - wts [y_rows * y_cols], xtwy [y_cols + y_rows ]
Output - out_beta_gpu [y_cols + y_rows ]

******************************************************************************/

void gpu_out_beta(int_t y_cols, int_t y_rows, double* wts, double* xtwy, double* out_beta_gpu)
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

	//Copy Ainv from Device to Host
	cudaError = cudaMemcpy(h_Ainv, d_Ainv, y_cols*sizeof(double), cudaMemcpyDeviceToHost);
	if(cudaError != cudaSuccess)	printf("\n cudaMemcpyDeviceToHost error = %s\n",cudaGetErrorString (cudaError));

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

	//Calculate C = Trans(B)
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

	//h_CAinvB = -Inv(D-C* Inv(A) * B)
	mat_negate(h_S, h_CAinvB, (y_rows-1), (y_rows-1));

	//Copy (-Inv(D-C* Inv(A) * B)) from Host to Device
	cudaError = cudaMemcpy(d_CAinvB, h_CAinvB, (y_rows-1)*(y_rows-1)*sizeof(double), cudaMemcpyHostToDevice);
	if(cudaError != cudaSuccess)	printf("\n cudaMemcpyHostToDevice error = %s\n",cudaGetErrorString (cudaError));

	//Calculate on Device Q = AinvB * CAinvB
	//AinvB [y_cols * (y_rows-1)]= A * Inv(B)
	//CAinvB [(y_rows-1) * (y_rows-1)] = (-Inv(D-C* Inv(A) * B))
	//Q [y_cols * (y_rows-1)] = (-A * Inv(B) * Inv(D-C* Inv(A) * B))	

	gpu_blas_mmul(d_AinvB, d_CAinvB, d_Q, (int)y_cols, (int)(y_rows-1), (int)(y_rows-1));

	//Copy Q from Device to Host
	//Q [y_cols * (y_rows-1)] = (-A * Inv(B) * Inv(D-C* Inv(A) * B))	

	cudaError = cudaMemcpy(h_Q, d_Q, y_cols*(y_rows-1)*sizeof(double), cudaMemcpyDeviceToHost);
	if(cudaError != cudaSuccess)	printf("\n cudaMemcpyDeviceToHost error = %s\n",cudaGetErrorString (cudaError));

	//Calculate on Host -> CAinv = trans(AinvB) 
	//AinvB [y_cols * (y_rows-1)] = Inv(A) * B
	//CAinv [(y_rows-1) * y_cols]= C * Inv(A)

	gpu_blas_mtrans(d_AinvB, d_CAinv, y_cols, y_rows-1);

	//Copy out_beta_gpu from Host to Device
	cudaError = cudaMemcpy(d_out_beta_gpu, out_beta_gpu, (y_rows+y_cols)*sizeof(double), cudaMemcpyHostToDevice);
	if(cudaError != cudaSuccess)	printf("\n cudaMemcpyHostToDevice error = %s\n",cudaGetErrorString (cudaError));
	
	//Copy xtwy from Host to Device
	cudaError = cudaMemcpy(d_xtwy, xtwy, (y_rows+y_cols)*sizeof(double), cudaMemcpyHostToDevice);
	if(cudaError != cudaSuccess)	printf("\n cudaMemcpyHostToDevice error = %s\n",cudaGetErrorString (cudaError));

	//Now we proceed to calculating P by dividing it into blocks (P_block)
	//P_block and xtwy are used to calculate corresponding values of out_beta
	//Things might get a bit difficult to follow from here on
	
	int_t i,j,k;

	//X_dptr is a double pointer which points to arrays that store Q divided into P_BLOCKS_NUM blocks
	//Q [y_cols * (y_rows-1)] = (-A * Inv(B) * Inv(D-C* Inv(A) * B))	

	//Similarly Y_dptr is a double pointer which points to arrays that store CAinv divided into P_BLOCKS_NUM blocks
	//CAinv [(y_rows-1) * y_cols]= C * Inv(A)
	
	double* temp;
	temp = (double*)calloc(y_cols*(y_rows-1), sizeof(double));
	double** X_dptr = (double**)calloc(P_BLOCKS_NUM, sizeof(double*));
	for(i=0;i<P_BLOCKS_NUM;i++){
		X_dptr[i] = temp+i*(P_BLOCK_SIZE*(y_rows-1));
	}
	
	//Fill X_dptr arrays
	for(k=0; k<P_BLOCKS_NUM; k++){
		for(j=0 ;j<(y_rows-1); j++){
			memcpy(X_dptr[k]+j*P_BLOCK_SIZE, h_Q+j*y_cols+k*P_BLOCK_SIZE, P_BLOCK_SIZE*sizeof(double));
		}
	}
	
	
	double **d_X_dptr;
	double **d_Y_dptr;		

	//d_X_dptr and d_Y_dptr are corresponding device pointers for X_dptr and Y_dptr
	d_X_dptr = (double**)calloc(P_BLOCKS_NUM, sizeof(double*));
	d_Y_dptr = (double**)calloc(P_BLOCKS_NUM, sizeof(double*));

	cudaError = cudaMalloc((void**)&temp, y_cols*(y_rows-1)*sizeof(double));	
	if(cudaError != cudaSuccess)	printf("\n cudaMalloc error = %s\n",cudaGetErrorString (cudaError));
	
	for(i=0;i<P_BLOCKS_NUM;i++){
		d_X_dptr[i] = temp+i*(P_BLOCK_SIZE*(y_rows-1));
	}

	for(i=0;i<P_BLOCKS_NUM;i++){
		d_Y_dptr[i] = d_CAinv+i*((y_rows-1)*P_BLOCK_SIZE);
	}

	//Copy data corresping to X_dptr i.e. Q from host to device
	//Data corresping to Y_dptr i.e. CAinv is already on device. Hence need not copy
	cudaError = cudaMemcpy(d_X_dptr[0], X_dptr[0], y_cols*(y_rows-1)*sizeof(double), cudaMemcpyHostToDevice);
	if(cudaError != cudaSuccess)	printf("\n cudaMemcpyHostToDevice error = %s\n",cudaGetErrorString (cudaError));

	//Perform multiplications by looping through all blocks
	//Calculate out_beta. Only out_beta is copied to host. Thus HostToDevice communication is signifcantly reduced.
	//Also host memory requirement is reduced
	for(j=0; j<P_BLOCKS_NUM; j++){
		for(i=0; i<P_BLOCKS_NUM; i++){ 
			
			gpu_blas_mmul(d_X_dptr[i], d_Y_dptr[j], d_P_block, (int)P_BLOCK_SIZE, (int)(y_rows-1), (int)P_BLOCK_SIZE);
			
			//Launch Kernel to calculate out_beta from P
			{
				dim3 dimGrid(ceil((P_BLOCK_SIZE)/((float)BLOCK_WIDTH)), 1, 1);
				dim3 dimBlock(BLOCK_WIDTH, 1, 1);
				Kernel_outbeta_P <<<dimGrid, dimBlock>>> (d_P_block, d_xtwy, d_out_beta_gpu, (int_t)P_BLOCK_SIZE,
									i*P_BLOCK_SIZE, j*P_BLOCK_SIZE);								
			}							
		}
	}

	//Copy out_beta from Device to Host
	cudaError = cudaMemcpy(out_beta_gpu, d_out_beta_gpu, (y_rows+y_cols)*sizeof(double), cudaMemcpyDeviceToHost);
	if(cudaError != cudaSuccess)	printf("\n cudaMemcpyDeviceToHost error = %s\n",cudaGetErrorString (cudaError));

	//Calculate part of out_beta from Ainv
	for(j=0; j<y_cols; j++){
		out_beta_gpu[j] += h_Ainv[j]*xtwy[j];
	}

	//Calculate part of out_beta from Q
	for(j=0; j<(y_rows-1); j++){
		for(i=0; i<y_cols; i++){
			out_beta_gpu[i] += h_Q[j*y_cols+i]*xtwy[y_cols+j];
			out_beta_gpu[y_cols+j] += h_Q[j*y_cols+i]*xtwy[i];
		}
	}
	
	//Calculate part of out_beta from S	
	for(j=0; j<(y_rows-1); j++){
		for(i=0; i<(y_rows-1); i++){
			out_beta_gpu[y_cols+i] += h_S[j*(y_rows-1)+i]*xtwy[y_cols+j];
		}
	}
	
	//done
	
	//free GPU memory
	cudaError = cudaFree(*d_X_dptr);
	if(cudaError != cudaSuccess)	printf("\n cudaFree error = %s\n",cudaGetErrorString (cudaError));

	//Free Host Memory
	free(*X_dptr);		
	free(X_dptr);	
	free(d_X_dptr);
	free(d_Y_dptr);

	cudaError = cudaGetLastError();
	if(cudaError != cudaSuccess)	printf("\n cudaGetLastError = %s\n",cudaGetErrorString (cudaError));

}

//rlm_fit_anova function with compute intensive code offloaded to gpu
void rlm_fit_anova_gpu(double* y, int_t y_rows, int_t y_cols, double* out_beta, double* resids, double* weights, double psi_k, int_t max_iter) {

	double acc = 1e-4;
	double scale = 0.0;

	double* old_resids = new double[y_rows * y_cols]();
	double* rowmeans = new double[y_rows]();

	P_BLOCK_SIZE = y_cols/P_BLOCKS_NUM;

	//double* xtwx = new double[(y_rows + y_cols - 1) * (y_rows + y_cols - 1)];
	double* xtwy = new double[y_rows + y_cols]();

	alloc_host_mem(y_cols,y_rows);
	alloc_device_mem(y_cols,y_rows);

	int_t rows = y_rows * y_cols;

	/* intially use equal weights */
	for (int_t i = 0; i < rows; i++) {
		weights[i] = 1.0;
	}


	/* starting matrix */
#if 0			
	for (int_t i = 0; i < y_rows; i++) {
		for (int_t j = 0; j < y_cols; j++) {
			resids[j * y_rows + i] = y[j * y_rows + i];
		}
	}
#endif

	/* sweep columns (ie chip effects) */
	for (int_t j = 0; j < y_cols; j++) {
		out_beta[j] = 0.0;
		double sumweights = 0.0;
		for (int_t i = 0; i < y_rows; i++) {
			out_beta[j] += weights[j * y_rows + i] * resids[j * y_rows + i];
			sumweights += weights[j * y_rows + i];
		}
		out_beta[j] /= sumweights;
		for (int_t i = 0; i < y_rows; i++) {
			resids[j * y_rows + i] = resids[j * y_rows + i] - out_beta[j];
		}
	}

	/* sweep rows  (ie probe effects) */
	for (int_t i = 0; i < y_rows; i++) {
		rowmeans[i] = 0.0;
		double sumweights = 0.0;
		for (int_t j = 0; j < y_cols; j++) {
			rowmeans[i] += weights[j * y_rows + i] * resids[j * y_rows + i];
			sumweights += weights[j * y_rows + i];
		}
		rowmeans[i] /= sumweights;
		for (int_t j = 0; j < y_cols; j++) {
			resids[j * y_rows + i] = resids[j * y_rows + i] - rowmeans[i];
		}
	}
	for (int_t i = 0; i < y_rows - 1; i++) {
		out_beta[i + y_cols] = rowmeans[i];
	}

	for (int_t iter = 0; iter < max_iter; iter++) {

		scale = med_abs(resids, rows) / 0.6745;

		if (fabs(scale) < 1e-10) {
			/*printf("Scale too small \n"); */
			break;
		}

		for (int_t i = 0; i < rows; i++) {
			old_resids[i] = resids[i];
		}

		for (int_t i = 0; i < rows; i++) {
			//weights[i] = PsiFunction.huber(resids[i] / scale, psi_k, 0);  /* psi_huber(resids[i]/scale,k,0); */
			weights[i] =  psi_huber(resids[i]/scale,psi_k,0);
		}

		/* weighted least squares */
#if 0
		int_t xtwx_length = (y_rows + y_cols - 1) * (y_rows + y_cols - 1);
		for (int_t i = 0; i < xtwx_length; i++) {
			xtwx[i] = 0.0;
		}
#endif

		/***************** GPU offload start ***************/

#if 0
		XTWX(y_rows, y_cols, weights, xtwx);
		XTWXinv(y_rows, y_cols, xtwx);
		XTWY(y_rows, y_cols, weights, y, xtwy);

		for (int_t i = 0; i < y_rows + y_cols - 1; i++) {
			out_beta[i] = 0.0;
			for (int_t j = 0; j < y_rows + y_cols - 1; j++) {
				out_beta[i] += xtwx[j * (y_rows + y_cols - 1) + i] * xtwy[j];
			}
		}
#else
		XTWY(y_rows, y_cols, weights, y, xtwy);
		for (int_t i=0;i < (y_rows+y_cols-1); i++){
			out_beta[i] = 0.0;
		}  
		gpu_out_beta(y_cols, y_rows, weights, xtwy, out_beta);		
#endif
		/***************** GPU offload end ***************/

		/* residuals */
		for (int_t i = 0; i < y_rows - 1; i++) {
			for (int_t j = 0; j < y_cols; j++) {
				resids[j * y_rows + i] = y[j * y_rows + i] - (out_beta[j] + out_beta[i + y_cols]);
			}
		}

		for (int_t j = 0; j < y_cols; j++) {
			double endprobe = 0.0;
			for (int_t i = 0; i < y_rows - 1; i++) {
				endprobe += out_beta[i + y_cols];
			}
			resids[j * y_rows + y_rows - 1] = y[j * y_rows + y_rows - 1] - (out_beta[j] - endprobe);
		}

		/*check convergence  based on residuals */
		double conv = irls_delta(old_resids, resids, rows);

		if (conv < acc) {
			/*    printf("Converged \n");*/
			break;

		}
	}

	delete[] old_resids;
	delete[] rowmeans;
	//delete[] xtwx;
	delete[] xtwy;

	free_device_mem();
	free_host_mem();
}

//Allocate host memory for matrices
static void alloc_host_mem(int_t y_cols, int_t y_rows)
{

#if 0
	h_Ainv = 	(double*)malloc(y_cols*sizeof(double));
	h_D = 		(double*)malloc((y_rows-1)*(y_rows-1)*sizeof(double));
	h_D_temp = 	(double*)malloc(CUSTOM_VAL*y_rows*sizeof(double));
	h_CAinvB = 	(double*)malloc((y_rows-1)*(y_rows-1)*sizeof(double));
	h_Q = 		(double*)malloc(y_cols*(y_rows-1)*sizeof(double));
	h_S = 		(double*)malloc((y_rows-1)*(y_rows-1)*sizeof(double));	
#else
	cudaError_t cudaError;

	cudaError = cudaHostAlloc((void**)&h_Ainv, 		(y_cols*sizeof(double)), 		cudaHostAllocDefault);
	if(cudaError != cudaSuccess)	printf("\n cudaHostAlloc h_Ainv error = %s\n",cudaGetErrorString (cudaError));

	cudaError = cudaHostAlloc((void**)&h_D, 		((y_rows-1)*(y_rows-1)*sizeof(double)), cudaHostAllocDefault);
	if(cudaError != cudaSuccess)	printf("\n cudaHostAlloc h_D error = %s\n",cudaGetErrorString (cudaError));

	cudaError = cudaHostAlloc((void**)&h_D_temp, 	(CUSTOM_VAL*y_rows*sizeof(double)), 	cudaHostAllocDefault);
	if(cudaError != cudaSuccess)	printf("\n cudaHostAlloc h_D_temp error = %s\n",cudaGetErrorString (cudaError));

	cudaError = cudaHostAlloc((void**)&h_CAinvB, 	((y_rows-1)*(y_rows-1)*sizeof(double)), cudaHostAllocDefault);
	if(cudaError != cudaSuccess)	printf("\n cudaHostAlloc h_CAinvB error = %s\n",cudaGetErrorString (cudaError));

	cudaError = cudaHostAlloc((void**)&h_Q, 		(y_cols*(y_rows-1)*sizeof(double)), 	cudaHostAllocDefault);
	if(cudaError != cudaSuccess)	printf("\n cudaHostAlloc h_Q error = %s\n",cudaGetErrorString (cudaError));

	cudaError = cudaHostAlloc((void**)&h_S, 		((y_rows-1)*(y_rows-1)*sizeof(double)), cudaHostAllocDefault);
	if(cudaError != cudaSuccess)	printf("\n cudaHostAlloc h_S error = %s\n",cudaGetErrorString (cudaError));

#endif
}

//Free allocated host memory
static void free_host_mem(void)
{

#if 0
	free(h_Ainv);
	free(h_D);	
	free(h_D_temp);
	free(h_CAinvB);
	free(h_Q);
	free(h_S);
#else
	cudaError_t cudaError;

	cudaError = cudaFreeHost(h_Ainv);
	if(cudaError != cudaSuccess)	printf("\n cudaFreeHost h_Ainv error = %s\n",cudaGetErrorString (cudaError));

	cudaError = cudaFreeHost(h_D);	
	if(cudaError != cudaSuccess)	printf("\n cudaFreeHost h_D error = %s\n",cudaGetErrorString (cudaError));

	cudaError = cudaFreeHost(h_D_temp);
	if(cudaError != cudaSuccess)	printf("\n cudaFreeHost h_D_temp error = %s\n",cudaGetErrorString (cudaError));

	cudaError = cudaFreeHost(h_CAinvB);
	if(cudaError != cudaSuccess)	printf("\n cudaFreeHost h_CAinvB error = %s\n",cudaGetErrorString (cudaError));

	cudaError = cudaFreeHost(h_Q);
	if(cudaError != cudaSuccess)	printf("\n cudaFreeHost h_Q error = %s\n",cudaGetErrorString (cudaError));

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

	cudaError = cudaMalloc((void**)&d_out_beta_gpu, (y_rows+y_cols)*sizeof(double));
	if(cudaError != cudaSuccess)	printf("\n cudaMalloc d_out_beta_gpu error = %s\n",cudaGetErrorString (cudaError));

	cudaError = cudaMalloc((void**)&d_xtwy, 		(y_rows+y_cols)*sizeof(double));
	if(cudaError != cudaSuccess)	printf("\n cudaMalloc d_xtwy error = %s\n",cudaGetErrorString (cudaError));

	cudaError = cudaMalloc((void**)&d_P_block, 		P_BLOCK_SIZE*P_BLOCK_SIZE*sizeof(double));
	if(cudaError != cudaSuccess)	printf("\n cudaMalloc d_P_block error = %s\n",cudaGetErrorString (cudaError));

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

	cudaError = cudaFree(d_out_beta_gpu);
	if(cudaError != cudaSuccess)	printf("\n cudaFree d_out_beta_gpu error = %s\n",cudaGetErrorString (cudaError));

	cudaError = cudaFree(d_xtwy);
	if(cudaError != cudaSuccess)	printf("\n cudaFree d_xtwy error = %s\n",cudaGetErrorString (cudaError));

	cudaError = cudaFree(d_P_block);
	if(cudaError != cudaSuccess)	printf("\n cudaFree d_P_block error = %s\n",cudaGetErrorString (cudaError));

}
