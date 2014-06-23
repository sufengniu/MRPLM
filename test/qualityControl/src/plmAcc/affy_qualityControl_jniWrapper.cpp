#include <jni.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include "define_t.h"
#include "affy_qualityControl_jniWrapper.h"

//#include "PLMsummarize_cpu.h"

void rlm_fit_anova_gpu(double* y, int_t y_rows, int_t y_cols, double* out_beta, double* resids, double* weights, double psi_k, int_t max_iter);
void rlm_compute_se_anova_gpu(int_t y_rows, int_t y_cols, double* resids, double* weights, double* se_estimates, double* varcov, double psi_k);
 
JNIEXPORT void JNICALL Java_affy_qualityControl_jniWrapper_rlm_1fit_1anova
(JNIEnv *env, jobject _obj, jdoubleArray inJNIY, jint numprobes, jint numchips, jdoubleArray outJNIBeta, jdoubleArray outJNIResids, jdoubleArray outJNIWeights, jdouble psi_k, jint max_iter) {
		
	// Convert the incoming JNI jintarray to C's jint[]
	jdouble *weights = env->GetDoubleArrayElements(outJNIWeights, NULL);	
	jdouble *y = env->GetDoubleArrayElements(inJNIY, NULL);
	jdouble *beta = env->GetDoubleArrayElements(outJNIBeta, NULL);
	jdouble *resids = env->GetDoubleArrayElements(outJNIResids, NULL);

//	printf("starting CPU job...\n"); fflush(stdout);
//	rlm_fit_anova_cpu(y, numprobes, numchips, beta, resids, weights, psi_k, max_iter);
	
//	printf("starting GPU job...\n"); fflush(stdout);

	/* application code here */
	rlm_fit_anova_gpu(y, numprobes, numchips, beta, resids, weights, psi_k, max_iter);


//	printf("GPU job done\n"); fflush(stdout);
		
	env->ReleaseDoubleArrayElements(inJNIY, y, JNI_ABORT); // release resources

//	env->SetDoubleArrayRegion(outJNIResids, 0, numprobes*numchips, resids);
//	env->SetDoubleArrayRegion(outJNIBeta, 0, numprobes+numchips, beta);  // copy
//	env->SetDoubleArrayRegion( outJNIWeights, 0, numprobes*numchips, weights);
		
	env->ReleaseDoubleArrayElements(outJNIBeta, beta, 0); // release resources
	env->ReleaseDoubleArrayElements(outJNIResids, resids, 0); // release resources
	env->ReleaseDoubleArrayElements(outJNIWeights, weights, 0); // release resources


}

JNIEXPORT void JNICALL Java_affy_qualityControl_jniWrapper_rlm_1compute_1se_1anova
(JNIEnv *env, jobject _obj, jint numprobes, jint numchips, jdoubleArray inJNIResids, jdoubleArray inJNIWeights, jdoubleArray outJNISe, jdoubleArray outJNIVarcov, jdouble psi_k) {

	jdouble *weights = env->GetDoubleArrayElements( inJNIWeights, NULL);
	jdouble *resids = env->GetDoubleArrayElements( inJNIResids, NULL);
	jdouble *se = env->GetDoubleArrayElements( outJNISe, NULL);
	jdouble *varcov;
//	jdouble *varcov = env->GetDoubleArrayElements(outJNIVarcov, NULL);


	/* cuda code here */
	rlm_compute_se_anova_gpu(numprobes, numchips, resids, weights, se, varcov, 1.345);

	env->ReleaseDoubleArrayElements( inJNIResids, resids, JNI_ABORT); // release resources
	env->ReleaseDoubleArrayElements( inJNIWeights, weights, JNI_ABORT);

//	env->SetDoubleArrayRegion( outJNISe, 0, numprobes+numchips, se);
//	env->SetDoubleArrayRegion( outJNIVarcov, 0, numprobes+numchips, varcov);

	// Convert the C's Native jdouble[] to JNI jdoublearray, and return
	env->ReleaseDoubleArrayElements( outJNISe, se, 0); // release resources
//	env->ReleaseDoubleArrayElements( outJNIVarcov, varcov, 0);

}

