#include "define_t.h"

void PLMsummarize_gpu(double* z, int_t numprobes, int_t numchips, double* results, double* SEresults /*double* affinities*/);

void rlm_fit_anova_gpu(double* y, int_t y_rows, int_t y_cols, double* out_beta, double* resids, double* weights, double psi_k, int_t max_iter);
void rlm_compute_se_anova_gpu(int_t y_rows, int_t y_cols, double* resids, double* weights, double* se_estimates, double* varcov, double psi_k);

