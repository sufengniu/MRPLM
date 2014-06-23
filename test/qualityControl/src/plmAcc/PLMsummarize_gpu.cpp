#include "PLMsummarize_gpu.h"

void PLMsummarize_gpu(double* z, int_t numprobes, int_t numchips, double* results, double* SEresults /*double* affinities*/) {	

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

	/*rlm_fit_anova_gpu*/
	rlm_fit_anova_gpu(y, numprobes, numchips, beta, resids, weights, 1.345, 20);

	/*rlm_compute_se_anova_gpu*/
	rlm_compute_se_anova_gpu(numprobes, numchips, resids, weights, se, varcov, 1.345);

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

