/*********************************************************************************
 **
 ** This file contains various functions used to calculate 
 **
 ** rlm_fit_anova and rlm_compute_se_anova
 **
 *********************************************************************************/

#include "PLMsummarize_utils.h"

double irls_delta(double* old, double* newv, int_t length) {
	double sum = 0.0;
	double sum2 = 0.0;
	double divisor = 1e-20;

	for (int_t i = 0; i < length; i++) {
		sum = sum + (old[i] - newv[i]) * (old[i] - newv[i]);
		sum2 = sum2 + old[i] * old[i];
	}

	if (sum2 >= divisor) {
		divisor = sum2;
	}

	return sqrt(sum / divisor);
}


void XTWY(int_t y_rows, int_t y_cols, double* wts, double* y, double* xtwy) {

	/*sweep columns (ie chip effects)*/
	for (int_t j = 0; j < y_cols; j++) {
		xtwy[j] = 0.0;
		for (int_t i = 0; i < y_rows; i++) {
			xtwy[j] += wts[j * y_rows + i] * y[j * y_rows + i];
		}
	}

	/*sweep rows  (ie probe effects)*/
	for (int_t i = y_rows - 1; i >= 0; i--) {
		xtwy[i + y_cols] = 0.0;
		for (int_t j = 0; j < y_cols; j++) {
			xtwy[i + y_cols] += wts[j * y_rows + i] * y[j * y_rows + i];
		}
		if (i < y_rows - 1) {
			xtwy[i + y_cols] = xtwy[i + y_cols] - xtwy[y_cols + y_rows - 1];
		}
	}
}

void XTWX(int_t y_rows, int_t y_cols, double* wts, double* xtwx) {

	int_t Msize = y_cols + y_rows - 1;

	/* diagonal elements of first part of matrix ie upper partition*/
	/*        for (int_t j = 0; j < y_cols; j++) {
		  for (int_t i = 0; i < y_rows-1; i++) {
		  xtwx[j * Msize + j] += wts[j * y_rows + i];
		  xtwx[(y_cols + i) * Msize + (y_cols + i)] += wts[j * y_rows + i];
		  for (int_t k = i; k < y_rows - 1; k++) {
		  xtwx[(y_cols + k) * Msize + (y_cols + i)] = xtwx[(y_cols + i) * Msize + (y_cols + k)] += wts[j * y_rows + (y_rows - 1)];
		  }
		  xtwx[j * Msize + (y_cols + i)] = xtwx[(y_cols + i) * Msize + j] = wts[j * y_rows + i] - wts[j * y_rows + (y_rows - 1)];
		  }
		  xtwx[j * Msize + j] += wts[j * y_rows + y_rows-1];
		  }
		  */
	double sum = 0;
	for (int_t j = 0; j < y_cols; j++) {
		for (int_t i = 0; i < y_rows - 1; i++) {
			xtwx[j * Msize + j] += wts[j * y_rows + i];
			xtwx[(y_cols + i) * Msize + (y_cols + i)] += wts[j * y_rows + i];
			xtwx[j * Msize + (y_cols + i)] = xtwx[(y_cols + i) * Msize + j] = wts[j * y_rows + i] - wts[j * y_rows + (y_rows - 1)];

		}
		sum += wts[j * y_rows + (y_rows - 1)];
		xtwx[j * Msize + j] += wts[j * y_rows + y_rows - 1];
	}

	for (int_t i = 0; i < y_rows - 1; i++) {
		for (int_t k = i; k < y_rows - 1; k++) {
			xtwx[(y_cols + k) * Msize + (y_cols + i)] = xtwx[(y_cols + i) * Msize + (y_cols + k)] += sum;
		}
	}

}

void XTWXinv(int_t y_rows, int_t y_cols, double* xtwx) {
	int_t Msize = y_cols + y_rows - 1;

	double* P = new double[y_cols]();
	double* RP = new double[y_cols * (y_rows - 1)]();
	double* RPQ = new double[(y_rows - 1) * (y_rows - 1)]();
	double* S = new double[(y_rows - 1) * (y_rows - 1)]();
	double* work = new double[(y_rows - 1) * (y_rows - 1)]();

	for (int_t j = 0; j < y_cols; j++) {
		for (int_t i = 0; i < y_rows - 1; i++) {
			RP[j * (y_rows - 1) + i] = xtwx[j * Msize + (y_cols + i)] * (1.0 / xtwx[j * Msize + j]);
		}
	}

	for (int_t i = 0; i < y_rows - 1; i++) {
		for (int_t j = i; j < y_rows - 1; j++) {
			for (int_t k = 0; k < y_cols; k++) {
				RPQ[j * (y_rows - 1) + i] += RP[k * (y_rows - 1) + j] * xtwx[k * Msize + (y_cols + i)];
			}
			RPQ[i * (y_rows - 1) + j] = RPQ[j * (y_rows - 1) + i];
			RPQ[j * (y_rows - 1) + i] = RPQ[i * (y_rows - 1) + j]
				= xtwx[(y_cols + i) * Msize + (y_cols + j)] - RPQ[i * (y_rows - 1) + j];
		}
	}

	/* Lets start making the inverse */
	//MatrixFunctions.Choleski_inverse(RPQ, S, work, y_rows - 1, false);
	Choleski_inverse(RPQ, S, work, (int)(y_rows-1), 0);

	for (int_t j = 0; j < y_cols; j++) {
		for (int_t i = 0; i < y_rows - 1; i++) {
			xtwx[j * Msize + (y_cols + i)] = 0.0;
			for (int_t k = 0; k < y_rows - 1; k++) {
				xtwx[j * Msize + (y_cols + i)] += -1.0 * (S[i * (y_rows - 1) + k]) * RP[j * (y_rows - 1) + k];
			}
			xtwx[(y_cols + i) * Msize + j] = xtwx[j * Msize + (y_cols + i)];
		}

		P[j] = 1.0 / xtwx[j * Msize + j];

		for (int_t i = j; i < y_cols; i++) {
			xtwx[i * Msize + j] = 0.0;
			for (int_t k = 0; k < y_rows - 1; k++) {
				xtwx[i * Msize + j] += RP[i * (y_rows - 1) + k] * xtwx[j * Msize + (y_cols + k)];
			}
			xtwx[i * Msize + j] *= -1.0;
			xtwx[j * Msize + i] = xtwx[i * Msize + j];
		}
		xtwx[j * Msize + j] += P[j];
	}


	for (int_t j = 0; j < y_rows - 1; j++) {
		for (int_t i = 0; i < y_rows - 1; i++) {
			xtwx[(y_cols + j) * Msize + (y_cols + i)] = S[j * (y_rows - 1) + i];
		}
	}

	delete[] P;
	delete[] RP;
	delete[] RPQ;
	delete[] S;
	delete[] work;

}

/*********************************************************************
 **
 ** double psi_huber(double u, double k,int deriv)
 **
 ** double u - data value
 ** double k - tuning constant
 ** int type - if 0 then return the evaluation of the weight function, if 1 returns the derivative
 **            other wise return psi itself
 **
 ** This function computes Hubers suggested PSI function.
 **
 *********************************************************************/

double psi_huber(double u, double k,int deriv){

	if (deriv == 0){
		if ( 1 < k/fabs(u)){
			return 1.0;
		} else {
			return  k/fabs(u);
		}
	} else if (deriv == 1){
		if (fabs(u) <= k){
			return 1.0;
		} else {
			return 0.0;
		}
	} else {
		if (fabs(u) <= k){
			return u;
		} else {
			if (u < 0){
				return -k;
			} else {
				return k;
			}
		}
	}
}


/**********************************************************
 **
 ** int sort_double(const void *a1,const void *a2)
 ** 
 ** a comparison function used when sorting doubles.
 **
 **********************************************************/

int sort_double(const double *a1,const double *a2){
	if (*a1 < *a2)
		return (-1);
	if (*a1 > *a2)
		return (1);
	return 0;
}

/**********************************************************************************
 **
 ** double med_abs(double *x, int length)
 **
 ** double *x - a vector of data
 ** int length - length of the vector.
 ** 
 ** returns the median of the absolute values.
 **
 ** computes the median of the absolute values of a given vector.
 **
 **********************************************************************************/

double med_abs(double *x, int length){
	int i;
	double med_abs;
	double *buffer = new double[length]();

	for (i = 0; i < length; i++)
		buffer[i] = fabs(x[i]);

	med_abs = median(buffer,length);

	delete[] buffer;
	return(med_abs);
}


/**************************************************************************
 **
 ** double median(double *x, int length)
 **
 ** double *x - vector
 ** int length - length of *x
 **
 ** returns the median of *x
 **
 *************************************************************************/

double  median(double *x, int length){
	/*  int i; */
	int half;
	double med;
	double *buffer = new double[length]();

	memcpy(buffer,x,length*sizeof(double));

	half = (length + 1)/2;

	qsort(buffer,length,sizeof(double), (int(*)(const void*, const void*))sort_double);  

	if (length % 2 == 1){
		med = buffer[half - 1];
	} 
	else {
		med = (buffer[half] + buffer[half-1])/2.0;
	}

	/*
	   rPsort(buffer, length, half-1);
	   med = buffer[half-1];
	   if (length % 2 == 0){
	   rPsort(buffer, length, half);
	   med = (med + buffer[half])/2.0;
	   }
	   */
	delete[] buffer;
	return med;
}


