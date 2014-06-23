#include "define_t.h"

double irls_delta(double* old, double* newv, int_t length);
void XTWY(int_t y_rows, int_t y_cols, double* wts, double* y, double* xtwy);
void XTWX(int_t y_rows, int_t y_cols, double* wts, double* xtwx);
void XTWXinv(int_t y_rows, int_t y_cols, double* xtwx);
double psi_huber(double u, double k,int deriv);

int sort_double(const double *a1,const double *a2);
double med_abs(double *x, int length);
double  median(double *x, int length);

int Choleski_inverse(double *X, double *Xinv, double *work, int n, int upperonly);

