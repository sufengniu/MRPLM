#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "define_t.h"

void XTWX(int_t y_rows, int_t y_cols, double *wts, double *xtwx);
void XTWY(int_t y_rows, int_t y_cols, double *wts,double *y, double *xtwy);
void XTWXinv(int_t y_rows, int_t y_cols, double *xtwx);


void wls_cpu(int_t y_cols, int_t y_rows, double* wts, double* y, double* out_beta_cpu){

	double *xtwx = (double *)calloc((y_rows+y_cols-1)*(y_rows+y_cols-1),sizeof(double));
	double *xtwy = (double *)calloc((y_rows+y_cols),sizeof(double));

	//XTWX start 
	XTWX(y_rows, y_cols, wts, xtwx); 
//	printf("\nXTWX done\n");	fflush(stdout);                                                                                  

	//XTWXinv 
	XTWXinv(y_rows, y_cols, xtwx); 
//	printf("\nXTWXinv done\n");	fflush(stdout);                                                                           

	//XTWY     
	XTWY(y_rows, y_cols, wts, y, xtwy);  
//	printf("\nXTWY done\n");	fflush(stdout);	

	int i, j;
	//OUT_BETA                                                                                                                
	for (i=0; i < y_rows+y_cols-1; i++){  
		out_beta_cpu[i] = 0.0;       
		for (j=0; j < y_rows+y_cols-1; j++){
			out_beta_cpu[i] += xtwx[j*(y_rows+y_cols -1)+i]*xtwy[j];	
		}        
	}

	free(xtwx);
	free(xtwy);

	return 0;
}

void XTWY(int_t y_rows, int_t y_cols, double *wts,double *y, double *xtwy){

	int_t i,j;

	/* sweep columns (ie chip effects) */

	for (j=0; j < y_cols; j++){
		xtwy[j] = 0.0;
		for (i=0; i < y_rows; i++){
			xtwy[j] += wts[j*y_rows + i]* y[j*y_rows + i];
		}
	}

	/* sweep rows  (ie probe effects) */

	for (i=0; i < y_rows; i++){
		xtwy[i+y_cols] = 0.0;
		for (j=0; j < y_cols; j++){
			xtwy[i+y_cols] += wts[j*y_rows + i]* y[j*y_rows + i]; 
		}
	}

	for (i=0; i < y_rows-1; i++){
		xtwy[i+y_cols] = xtwy[i+y_cols] - xtwy[y_cols+y_rows-1];
	}

}


void XTWX(int_t y_rows, int_t y_cols, double *wts, double *xtwx){

	int_t Msize = y_cols +y_rows-1;
	int_t i,j,k;

	//struct timeval start, end;
	//long utime, seconds, useconds;

	/* diagonal elements of first part of matrix ie upper partition */
	for (j =0; j < y_cols;j++){
		for (i=0; i < y_rows; i++){
			xtwx[j*Msize + j]+=wts[j*y_rows + i];
		}
	}

	/* diagonal portion of lower partition matrix: diagonal elements*/ 
	for (j =0; j < y_cols;j++){
		for (i = 0; i < y_rows-1;i++){
			xtwx[(y_cols +i)*Msize + (y_cols +i)]+= wts[j*y_rows + i];
		}
	}

	/* diagonal portion of lower partition matrix: off diagonal elements*/ 
	for (j =0; j < y_cols;j++){
		for (i = 0; i < y_rows-1;i++){
			for (k=i ; k <  y_rows-1;k++){
				xtwx[(y_cols +k)*Msize + (y_cols +i)] = xtwx[(y_cols +i)*Msize + (y_cols +k)]+= wts[j*y_rows + (y_rows-1)];
			}
		}
	}

	/* the two other portions of the matrix */
	for (j =0; j < y_cols;j++){
		for (i= 0; i < y_rows-1;i++){
			xtwx[j*Msize + (y_cols + i)] = xtwx[(y_cols + i)*Msize + j] = wts[j*y_rows + i] - wts[j*y_rows + (y_rows-1)];
		}
	}
}

void XTWXinv(int_t y_rows, int_t y_cols, double *xtwx){
	int_t i,j,k;
	int_t Msize = y_cols +y_rows-1;
	double *P= (double *)calloc(y_cols,sizeof(double));
	double *RP = (double *)calloc(y_cols*(y_rows-1),sizeof(double));
	double *RPQ = (double *)calloc((y_rows-1)*(y_rows-1),sizeof(double));
	double *S = (double *)calloc((y_rows-1)*(y_rows-1),sizeof(double));
	double *work = (double *)calloc((y_rows-1)*(y_rows-1),sizeof(double));

	for (j=0;j < y_cols;j++){
		for (i=0; i < y_rows -1; i++){
			RP[j*(y_rows-1) + i] = xtwx[j*Msize + (y_cols + i)]*(1.0/xtwx[j*Msize+j]);
		}
	} 

	for (i=0; i < y_rows -1; i++){
		for (j=i;j <  y_rows -1; j++){
			for (k=0; k < y_cols;k++){
				RPQ[j*(y_rows-1) + i] +=  RP[k*(y_rows-1) + j]*xtwx[k*Msize + (y_cols + i)];
			}
			RPQ[i*(y_rows-1) + j] = RPQ[j*(y_rows-1) + i];
		}
	}

	for (j=0; j < y_rows-1;j++){
		for (i=j; i < y_rows-1;i++){
			RPQ[i*(y_rows-1) + j] = RPQ[j*(y_rows-1)+i] =  xtwx[(y_cols + j)*Msize + (y_cols + i)] - RPQ[j*(y_rows-1) + i]; 
		}
	} 


	/*for (i =0; i<  y_rows-1; i++){
	  for (j=0; j <  y_cols; j++){ 
	  printf("%4.4f ",RP[j*(y_rows-1) + i]);
	  }
	  printf("\n");
	  }

	  for (j=0;j <  y_rows -1; j++){
	  for (i=0; i < y_rows -1; i++){
	  printf("%4.4f ",RPQ[j*(y_rows-1) + i]);
	  }
	  printf("\n");
	  }


	  for (i=0; i < y_rows -1; i++){
	  for (j=0;j <  y_rows -1; j++){
	  printf("%4.4f ",S[j*(y_rows-1) + i]);
	  }
	  printf("\n");
	  }
	  */

	/* Lets start making the inverse */

	Choleski_inverse(RPQ, S, work, (int)(y_rows-1), 0);

	for (j=0; j< y_cols;j++){
		for (i=0; i < y_rows -1; i++){
			xtwx[j*Msize + (y_cols + i)] = 0.0;
			for (k=0; k < y_rows -1; k++){
				xtwx[j*Msize + (y_cols + i)]+= -1.0*(S[i*(y_rows-1) + k])*RP[j*(y_rows-1) + k];
			}
			xtwx[(y_cols + i)*Msize + j]=xtwx[j*Msize + (y_cols + i)];
		}
	}


	for (j=0;j < y_cols;j++){
		P[j] = 1.0/xtwx[j*Msize+j];
	} 


	for (j=0; j < y_cols; j++){
		for (i=j; i < y_cols;i++){
			xtwx[i*Msize + j]=0.0;
			for (k=0;k < y_rows-1; k++){
				xtwx[i*Msize + j]+= RP[i*(y_rows-1) + k]*xtwx[j*Msize + (y_cols + k)];
			}
			xtwx[i*Msize + j]*=-1.0;
			xtwx[j*Msize + i] =  xtwx[i*Msize + j];
		}
		xtwx[j*Msize + j]+=P[j];
	}


	for (j=0; j < y_rows-1;j++){
		for (i=0; i < y_rows-1;i++){
			xtwx[(y_cols + j)*Msize + (y_cols + i)] = S[j*(y_rows-1)+i];
		}
	}

	free(P);
	free(work);
	free(RP);
	free(RPQ);
	free(S);

}

