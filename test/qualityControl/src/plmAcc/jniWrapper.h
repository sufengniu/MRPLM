#ifndef JNIWRAPPER_H_
#define JNIWRAPPER_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "define_t.h"
#include "wls_acc.h"

extern void wls_acc(int_t y_cols, int_t y_rows, double* wts, double* y, double* out_beta_gpu);

JNIEXPORT void JNICALL Java_jniWrapper_wlsAcc(JNIEnv *env, jobject thisObj, jdoubleArray inJNIWeights, jdoubleArray outJNIY, jdoubleArray outJNIBeta, jint y_rows, jint y_cols);
JNIEXPORT jdoubleArray JNICALL Java_jniWrapper_seAcc(JNIEnv *env, jobject thisObj, jdoubleArray inJNIArray, jint y_rows, jint y_cols);

#endif
