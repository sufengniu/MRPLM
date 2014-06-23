#ifndef _PLMSUMMARIZE_UTILS_H_
#define _PLMSUMMARIZE_UTILS_H_

#include "define_t.h"
#include "PLMsummarize_utils.h"

void PLMsummarize_cpu(double* z, int_t numprobes, int_t numchips, double* results, double* SEresults /*double* affinities*/);
void rlm_fit_anova(double* y, int_t y_rows, int_t y_cols, double* out_beta, double* resids, double* weights, double psi_k, int_t max_iter);
void rlm_compute_se_anova(int_t y_rows, int_t y_cols, double* resids, double* weights, double* se_estimates, double* varcov, double psi_k);

#endif
