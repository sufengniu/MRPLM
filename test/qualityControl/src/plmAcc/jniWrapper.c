#include <jni.h>
#include <stdio.h>

#include "jniWrapper.h"

int main()
{

	return 0;
}


// output: out_beta
JNIEXPORT void JNICALL Java_jniWrapper_wlsAcc
(JNIEnv *env, jobject thisObj, jdoubleArray inJNIWeights, jdoubleArray outJNIY, jdoubleArray outJNIBeta, jint y_rows, jint y_cols) {

	// Convert the incoming JNI jintarray to C's jint[]
	jdouble *wts = (*env)->GetDoubleArrayElements(env, inJNIWeights, NULL);
	
	jsize length = (*env)->GetArrayLength(env, inJNIWeights);

	jdouble *y = (*env)->GetDoubleArrayElements(env, outJNIY, NULL);
	jdouble *out_beta = (*env)->GetDoubleArrayElements(env, outJNIBeta, NULL);
	

	/* application code here */
	wls_gpu(y_cols, y_rows, wts, y, out_beta);

	//(*env)->ReleaseDoubleArrayElements(env, inJNIArray, inCArray, 0); // release resources
	
	// Convert the C's Native jdouble[] to JNI jdoublearray, and return
/*	jdoubleArray outJNIArray = (*env)->NewDoubleArray(env, y_rows+y_cols);  // allocate
	if (NULL == outJNIArray) return NULL;
	(*env)->SetDoubleArrayRegion(env, outJNIArray, 0 , y_rows+y_cols, outArray);  // copy
	return outJNIArray;
*/

}


JNIEXPORT jdoubleArray JNICALL Java_jniWrapper_seAcc
(JNIEnv *env, jobject thisObj, jdoubleArray inJNIArray, jint y_rows, jint y_cols) {
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

