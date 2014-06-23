#include "PLMsummarize_cpu.h"

void PLMsummarize_cpu(double* z, int_t numprobes, int_t numchips, double* results, double* SEresults /*double* affinities*/) {	

	double* weights = new double[numprobes * numchips]();
	double* resids = new double[numprobes * numchips]();
	double* y = new double[numprobes * numchips]();
	double* beta = new double[numprobes + numchips]();
	double* se = new double[numprobes + numchips]();

	double* varcov = NULL;

	for (int_t j = 0; j < numchips; j++) {
		for (int_t i = 0; i < numprobes; i++) {
			resids[j * numprobes + i] = y[j * numprobes + i] = z[i*numchips + j];
		}
	}

	/*rlm_fit_anova*/
	rlm_fit_anova(y, numprobes, numchips, beta, resids, weights, 1.345, 20);

	/*rlm_compute_se_anova*/
	rlm_compute_se_anova(numprobes, numchips, resids, weights, se, varcov, 1.345);

	for (int_t j = 0; j < numchips; j++) {
		results[j] = beta[j];
		SEresults[j] = se[j];
	}

	for (int_t j = 0; j < numchips; j++) {
		for (int_t i = 0; i < numprobes; i++) {
			z[i*numchips + j] = resids[j * numprobes + i];
		}
	}

	delete[] weights;
	delete[] resids;
	delete[] y;
	delete[] beta;
	delete[] se;
}

void rlm_fit_anova(double* y, int_t y_rows, int_t y_cols, double* out_beta, double* resids, double* weights, double psi_k, int_t max_iter) {


	//struct timeval start, end;
	//long utime;	

	double acc = 1e-4;
	double scale = 0.0;

	double* old_resids = new double[y_rows * y_cols]();
	double* rowmeans = new double[y_rows]();

	double* xtwx = new double[(y_rows + y_cols - 1) * (y_rows + y_cols - 1)]();
	double* xtwy = new double[y_rows + y_cols]();

	int_t rows = y_rows * y_cols;

	/* intially use equal weights */
	for (int_t i = 0; i < rows; i++) {
		weights[i] = 1.0;
	}


	/* starting matrix */
	/*			
				for (int_t i = 0; i < y_rows; i++) {
				for (int_t j = 0; j < y_cols; j++) {
				resids[j * y_rows + i] = y[j * y_rows + i];
				}
				}
				*/

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

		int_t xtwx_length = (y_rows + y_cols - 1) * (y_rows + y_cols - 1);
		for (int_t i = 0; i < xtwx_length; i++) {
			xtwx[i] = 0.0;
		}

		/***************** GPU offload start ***************/

		XTWX(y_rows, y_cols, weights, xtwx);
		XTWXinv(y_rows, y_cols, xtwx);
		XTWY(y_rows, y_cols, weights, y, xtwy);

		for (int_t i = 0; i < y_rows + y_cols - 1; i++) {
			out_beta[i] = 0.0;
			for (int_t j = 0; j < y_rows + y_cols - 1; j++) {
				out_beta[i] += xtwx[j * (y_rows + y_cols - 1) + i] * xtwy[j];
			}
		}

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
	delete[] xtwx;
	delete[] xtwy;
}

void rlm_compute_se_anova(int_t y_rows, int_t y_cols, double* resids, double* weights, double* se_estimates, double* varcov, double psi_k) {

	//struct timeval start, end;
	//long utime;	

	//double k1 = psi_k;   /*  was 1.345; */
	//double sumpsi2=0.0;  /* sum of psi(r_i)^2 */
	/*  double sumpsi=0.0; */
	//double sumderivpsi=0.0; /* sum of psi'(r_i) */
	//double Kappa=0.0;      /* A correction factor */
	//double scale=0.0;
	int_t n = y_rows * y_cols;
	int_t p = y_rows + y_cols - 1;
	double* XTX = new double[p * p]();

	//double [] W = new double [p*p];
	//double [] work = new double [p*p];
	double RMSEw = 0.0;
	//double vs=0.0,m,varderivpsi=0.0;
	//double [] W_tmp= new double [n];


	for (int_t i = 0; i < n; i++) {
		RMSEw += weights[i] * resids[i] * resids[i];
	}

	RMSEw = sqrt(RMSEw / (double) (n - p));
	//residSE[0] =  RMSEw;

	/***************** GPU offload part ***************/		
	XTWX(y_rows, y_cols, weights, XTX);
	if (y_rows > 1) {
		XTWXinv(y_rows, y_cols, XTX);
	} 
	else {
		for (int_t i = 0; i < p; i++) {
			XTX[i * p + i] = 1.0 / XTX[i * p + i];
		}
	}
	/***************** GPU offload end ***************/	

	for (int_t i = 0; i < p; i++) {
		se_estimates[i] = RMSEw * sqrt(XTX[i * p + i]);
	}


	if (varcov != NULL){
		for (int_t i = 0; i < p; i++)
			for (int_t j = i; j < p; j++)
				varcov[j * p + i] = RMSEw * RMSEw * XTX[j * p + i];
	}
	delete[] XTX;
}
