#include <jni.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include "affy_qualityControl_jniWrapper.h"
 
JNIEXPORT void JNICALL Java_affy_qualityControl_jniWrapper_wlsAcc
(JNIEnv *env, jobject _obj, jdoubleArray inJNIWeights, jdoubleArray outJNIY, jdoubleArray outJNIBeta, jint y_rows, jint y_cols) {
		
	// Convert the incoming JNI jintarray to C's jint[]
	jdouble *wts = (*env)->GetDoubleArrayElements(env, inJNIWeights, NULL);	
	jsize length = (*env)->GetArrayLength(env, inJNIWeights);

	jdouble *y = (*env)->GetDoubleArrayElements(env, outJNIY, NULL);
	jdouble *out_beta = (*env)->GetDoubleArrayElements(env, outJNIBeta, NULL);
	
/*	jclass _class = (*env)->GetObjectClass(env, _obj);

	// Get the field ID of instance array y
	jfieldID fidy = (*env)->GetFieldID(env, _class, "y", "[D");
	if (fidy == NULL) return;
	
	jdoubleArray fieldy = (jdoubleArray)(*env)->GetObjectField(env, _obj, fidy);		
	jdouble *y = (*env)->GetDoubleArrayElements(env, fieldy, NULL);

	// Get the field ID of instance array beta
	jfieldID fidbeta = (*env)->GetFieldID(env, _class, "beta", "[D");
	if (fidbeta == NULL) return;

	jdoubleArray fieldbeta = (jdoubleArray)(*env)->GetObjectField(env, _obj, fidbeta);
	jdouble *beta = (*env)->GetDoubleArrayElements(env, fieldbeta, NULL);
*/
	
//	printf("starting CPU job...\n"); fflush(stdout);
//	wls_cpu(y_cols, y_rows, wts, y, out_beta);
	
//	printf("starting GPU job...\n"); fflush(stdout);

	/* application code here */
	wls_gpu(y_cols, y_rows, wts, y, out_beta);

//	printf("GPU job done\n"); fflush(stdout);
		
	(*env)->ReleaseDoubleArrayElements(env, inJNIWeights, wts, 0); // release resources

	//check for out_beta
/*	double out_beta_err = 0.0;
	int j;
	for(j=0; j<(y_rows+y_cols-1); j++){
		out_beta_err+=fabs((out_beta[j] - out_beta_cpu[j])/out_beta[j]);
	}
	out_beta_err/=(y_rows+y_cols-1);
	printf("\nError Out_beta_error = %e\n",out_beta_err); 	
*/

//	(*env)->SetDoubleArrayRegion(env, outJNIY, 0, y_rows*y_cols, y);
	(*env)->SetDoubleArrayRegion(env, outJNIBeta, 0, y_rows+y_cols, out_beta);  // copy

	(*env)->ReleaseDoubleArrayElements(env, outJNIY, y, 0); // release resources
	(*env)->ReleaseDoubleArrayElements(env, outJNIBeta, out_beta, 0); // release resources
	
}

JNIEXPORT void JNICALL Java_affy_qualityControl_jniWrapper_seAcc
(JNIEnv *env, jobject _obj, jdoubleArray inJNIWeights, jdoubleArray outJNIXTX, jint y_rows, jint y_cols) {
	// Convert the incoming JNI jintarray to C's jint[]
	jdouble *weights = (*env)->GetDoubleArrayElements(env, inJNIWeights, NULL);
	if (NULL == inCArray) return NULL;	
	
	jdouble *xtx = (*env)->GetDoubleArrayElements(env, outJNIXTX, NULL);
	
	/* cuda code here */
	se_gpu(y_cols, y_rows, weights, xtx);

	(*env)->ReleaseDoubleArrayElements(env, inJNIWeights, weights, 0); // release resources
	(*env)->ReleaseDoubleArrayElements(env, outJNIXTX, xtx, 0);

}

