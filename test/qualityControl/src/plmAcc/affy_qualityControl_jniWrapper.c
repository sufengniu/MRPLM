#include <jni.h>
#include <stdio.h>
#include <stdlib.h>

#include "affy_qualityControl_jniWrapper.h"
 
JNIEXPORT void JNICALL Java_affy_qualityControl_jniWrapper_wlsAcc
(JNIEnv *env, jobject _obj, jdoubleArray inJNIWeights, jdoubleArray outJNIY, jdoubleArray outJNIBeta, jint y_rows, jint y_cols) {
		
	// Convert the incoming JNI jintarray to C's jint[]
	jdouble *wts = (*env)->GetDoubleArrayElements(env, inJNIWeights, NULL);	
	jsize length = (*env)->GetArrayLength(env, inJNIWeights);

//	jdouble *y = (*env)->GetDoubleArrayElements(env, outJNIY, NULL);
//	jdouble *out_beta = (*env)->GetDoubleArrayElements(env, outJNIBeta, NULL);
	
	jdouble *y = (jdouble *)malloc(y_rows*y_cols*sizeof(jdouble));
	jdouble *out_beta = (jdouble *)malloc((y_rows+y_cols)*sizeof(jdouble));

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
	
	printf("starting GPU job...\n");

	/* application code here */
	wls_gpu(y_cols, y_rows, wts, y, out_beta);

	(*env)->ReleaseDoubleArrayElements(env, inJNIWeights, wts, 0); // release resources

	(*env)->SetDoubleArrayRegion(env, outJNIY, 0, y_rows*y_cols, y);
	(*env)->SetDoubleArrayRegion(env, outJNIBeta, 0, y_rows+y_cols, out_beta);  // copy
	
}


JNIEXPORT jdoubleArray JNICALL Java_affy_qualityControl_jniWrapper_seAcc
(JNIEnv *env, jobject _obj, jdoubleArray inJNIArray, jint y_rows, jint y_cols) {
	// Convert the incoming JNI jintarray to C's jint[]
	jdouble *inCArray = (*env)->GetDoubleArrayElements(env, inJNIArray, NULL);
	if (NULL == inCArray) return NULL;
	jsize length = (*env)->GetArrayLength(env, inJNIArray);

	/* cuda code here */
	double outArray[y_rows*y_cols];

	(*env)->ReleaseDoubleArrayElements(env, inJNIArray, inCArray, 0); // release resources

	// Convert the C's Native jdouble[] to JNI jdoublearray, and return
	jdoubleArray outJNIArray = (*env)->NewDoubleArray(env, y_rows+y_cols-1);  // allocate
	if (NULL == outJNIArray) return NULL;
	(*env)->SetDoubleArrayRegion(env, outJNIArray, 0 , y_rows+y_cols-1, outArray);  // copy
	return outJNIArray;

}

